"""
Regression tests for background dataset warm-up.

The deploy failure this prevents: the container ran a blocking prewarm
before uvicorn bound, so the platform healthcheck could not pass until
every download finished. The auxiliary layers pushed measured cold start
to ~848s against a 900s healthcheck budget -- 94% consumed -- from a home
connection, with the MAF loader retrying up to 4x per flaky file.

The invariant worth pinning is not "warm-up works" but "nothing that
answers a request waits on it".
"""
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from app.services import warmup as warmup_module
from app.services.warmup import WarmupTracker, warm_datasets


class _FakeDataset:
    def __init__(self, n_samples: int):
        import pandas as pd

        self.matrix = pd.DataFrame({f"S{i}": [1.0] for i in range(n_samples)}, index=["G"])


@pytest.fixture(autouse=True)
def _fresh_tracker(monkeypatch):
    """Each test gets its own tracker -- the module-level one is global
    state and would leak between tests."""
    monkeypatch.setattr(warmup_module, "tracker", WarmupTracker())


def test_states_advance_pending_to_ready():
    warmup_module.tracker.register(["a", "b"])
    assert {d["state"] for d in warmup_module.tracker.snapshot()["datasets"]} == {"pending"}

    warm_datasets(lambda ds: _FakeDataset(3), ["a", "b"])

    snap = warmup_module.tracker.snapshot()
    assert snap["finished"] is True
    assert snap["n_ready"] == 2
    assert snap["n_failed"] == 0
    assert all(d["n_samples"] == 3 for d in snap["datasets"])


def test_one_failure_does_not_stop_the_others():
    """A single unreachable source must not strand every other dataset --
    that would turn one flaky upstream into a fully dead deployment."""

    def loader(dataset_id):
        if dataset_id == "broken":
            raise RuntimeError("GDC unreachable")
        return _FakeDataset(2)

    warm_datasets(loader, ["broken", "fine"])

    snap = warmup_module.tracker.snapshot()
    states = {d["dataset_id"]: d["state"] for d in snap["datasets"]}
    assert states == {"broken": "failed", "fine": "ready"}
    assert snap["n_ready"] == 1 and snap["n_failed"] == 1
    broken = next(d for d in snap["datasets"] if d["dataset_id"] == "broken")
    assert "GDC unreachable" in broken["error"]


def test_begin_is_single_shot():
    """A double-start would run every download twice."""
    assert warmup_module.tracker.begin() is True
    assert warmup_module.tracker.begin() is False


def test_start_background_warmup_does_not_block_the_caller():
    """The whole point: the caller (app startup) must return immediately,
    long before a slow loader finishes."""
    release = threading.Event()

    def slow_loader(dataset_id):
        release.wait(timeout=10)
        return _FakeDataset(1)

    start = time.monotonic()
    thread = warmup_module.start_background_warmup(slow_loader, ["slow"])
    elapsed = time.monotonic() - start

    try:
        # Returning at all while the loader is still blocked is the
        # assertion; the generous bound just avoids CI flakiness.
        assert elapsed < 1.0
        assert warmup_module.tracker.snapshot()["finished"] is False
    finally:
        release.set()
        if thread:
            thread.join(timeout=10)

    assert warmup_module.tracker.snapshot()["n_ready"] == 1


def test_snapshot_is_readable_while_a_load_is_in_flight():
    """Request handlers read this map from another thread while the
    warm-up thread writes it."""
    release = threading.Event()
    seen = []

    def slow_loader(dataset_id):
        seen.append(warmup_module.tracker.state_of(dataset_id))
        release.wait(timeout=10)
        return _FakeDataset(1)

    thread = warmup_module.start_background_warmup(slow_loader, ["x"])
    try:
        deadline = time.monotonic() + 5
        while warmup_module.tracker.state_of("x") != "loading" and time.monotonic() < deadline:
            time.sleep(0.01)
        # Readable, and honest about being mid-flight.
        assert warmup_module.tracker.state_of("x") == "loading"
        assert warmup_module.tracker.snapshot()["n_ready"] == 0
    finally:
        release.set()
        if thread:
            thread.join(timeout=10)

    assert seen == ["loading"]
