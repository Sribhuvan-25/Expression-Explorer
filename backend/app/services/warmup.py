"""
Background dataset warm-up, and the readiness state the API reports.

**Why this is not done before the server binds.** It used to be: the
container ran `scripts/prewarm_cache.py` to completion and only then
started uvicorn, so `/health` stayed dead until every dataset was on
disk. That made deploy success a race against a download:

- The first calibration of that race set Railway's `healthcheckTimeout`
  to 30s, which killed every deploy (Build and Deploy succeeded,
  Healthcheck failed at exactly 30s). It was raised to 900s.
- Then the auxiliary layers landed (METHODS.md §8), adding a 739-file
  GDC MAF pull plus a 198MB CRISPR matrix to the same startup path.
  Measured cold cost is now ~848s against that 900s budget -- 94%
  consumed, ~52s of headroom, from a *home* connection. GDC is
  typically slower to datacenter IPs, and the MAF loader retries up to
  4x per flaky file.

Raising the timeout again would just restart the same race, and the next
dataset would lose it. So the race is removed instead: uvicorn binds
immediately, `/health` answers at once, and the data loads on a daemon
thread behind it. Deploy success no longer depends on download speed at
all.

The cost is that a dataset can be *absent but coming* -- a state the API
must state plainly rather than paper over, which is what `status()` and
the `warming` flag on `/datasets` are for. Answering "no such dataset"
while it is still downloading would be a lie of exactly the kind
METHODS.md §6.1 exists to prevent.
"""
from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass, field
from typing import Literal

logger = logging.getLogger(__name__)

State = Literal["pending", "loading", "ready", "failed"]


@dataclass
class DatasetWarmup:
    dataset_id: str
    state: State = "pending"
    seconds: float | None = None
    error: str | None = None
    n_samples: int | None = None


@dataclass
class WarmupTracker:
    """Thread-safe readiness map. Written by the warm-up thread, read by
    request handlers, so every mutation takes the lock."""

    _lock: threading.Lock = field(default_factory=threading.Lock)
    _datasets: dict[str, DatasetWarmup] = field(default_factory=dict)
    _started: bool = False
    _finished: bool = False

    def register(self, dataset_ids: list[str]) -> None:
        with self._lock:
            for dataset_id in dataset_ids:
                self._datasets.setdefault(dataset_id, DatasetWarmup(dataset_id))

    def mark(self, dataset_id: str, state: State, **fields) -> None:
        with self._lock:
            entry = self._datasets.setdefault(dataset_id, DatasetWarmup(dataset_id))
            entry.state = state
            for key, value in fields.items():
                setattr(entry, key, value)

    def state_of(self, dataset_id: str) -> State | None:
        with self._lock:
            entry = self._datasets.get(dataset_id)
            return entry.state if entry else None

    def snapshot(self) -> dict:
        with self._lock:
            datasets = [
                {
                    "dataset_id": d.dataset_id,
                    "state": d.state,
                    "seconds": round(d.seconds, 1) if d.seconds is not None else None,
                    "error": d.error,
                    "n_samples": d.n_samples,
                }
                for d in self._datasets.values()
            ]
            finished = self._finished
        ready = sum(1 for d in datasets if d["state"] == "ready")
        failed = sum(1 for d in datasets if d["state"] == "failed")
        return {
            "finished": finished,
            "n_ready": ready,
            "n_failed": failed,
            "n_total": len(datasets),
            "datasets": sorted(datasets, key=lambda d: d["dataset_id"]),
        }

    def begin(self) -> bool:
        """True if this call started the run; False if one is already going.
        Guards against a double-start (e.g. a reload-mode worker importing
        the module twice), which would duplicate every download."""
        with self._lock:
            if self._started:
                return False
            self._started = True
            return True

    def finish(self) -> None:
        with self._lock:
            self._finished = True


tracker = WarmupTracker()


def warm_datasets(loader_for, dataset_ids: list[str]) -> None:
    """Load every dataset in turn, recording state as it goes.

    A failure is logged and skipped rather than raised: one unreachable
    source must not stop the others from warming, and the lazy
    load-on-first-request path still exists as a fallback.
    """
    tracker.register(dataset_ids)
    for dataset_id in dataset_ids:
        tracker.mark(dataset_id, "loading")
        start = time.monotonic()
        try:
            dataset = loader_for(dataset_id)
        except Exception as exc:  # noqa: BLE001 - deliberately broad, see docstring
            elapsed = time.monotonic() - start
            logger.warning(
                "[warmup] %s FAILED after %.0fs (%s) -- will lazy-load on first request",
                dataset_id,
                elapsed,
                exc,
            )
            tracker.mark(dataset_id, "failed", seconds=elapsed, error=str(exc)[:300])
            continue
        elapsed = time.monotonic() - start
        logger.info("[warmup] %s ready (%d samples, %.0fs)", dataset_id, dataset.matrix.shape[1], elapsed)
        tracker.mark(dataset_id, "ready", seconds=elapsed, n_samples=int(dataset.matrix.shape[1]))
    tracker.finish()


def start_background_warmup(loader_for, dataset_ids: list[str]) -> threading.Thread | None:
    """Kick off `warm_datasets` on a daemon thread. Returns None if a
    warm-up is already running."""
    tracker.register(dataset_ids)
    if not tracker.begin():
        return None
    thread = threading.Thread(
        target=warm_datasets,
        args=(loader_for, dataset_ids),
        name="dataset-warmup",
        daemon=True,
    )
    thread.start()
    return thread
