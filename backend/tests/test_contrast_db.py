"""
Design contrasts survive ingest -> database -> reload unchanged.

The contrast statistics are computed once, at ingest, and every later
read comes from the stored parquet. A lossy round trip (dropped untested
genes, a re-ordered index, a spec field lost) would change numbers a
user reads with no error anywhere, so the comparison is exact.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import registry
from app.db import session as db_session
from app.db.models import Base, ContrastRow, DatasetRow
from app.ingest import gse324984
from tests.test_gse324984 import write_dataset  # noqa: F401  (fixture)

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "ingest_dataset.py"


@pytest.fixture()
def ingest_script(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_URL", f"sqlite:///{tmp_path / 'contrast.db'}")
    db_session.reset_engine()
    Base.metadata.create_all(db_session.get_engine())
    spec = importlib.util.spec_from_file_location("ingest_dataset_under_test", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "_matrix_dir", lambda dataset_id: tmp_path / "datasets" / dataset_id)
    yield module
    db_session.reset_engine()


def test_contrasts_round_trip_exactly(write_dataset, ingest_script):  # noqa: F811
    write_dataset()
    original = gse324984.load()
    ingest_script.ingest("gse324984")

    with db_session.session_scope() as s:
        row = s.get(DatasetRow, "gse324984")
        assert {c.contrast_id for c in row.contrasts} == {"c42b_mycl_oe", "h660_mycl_kd"}

    reloaded = registry.get_descriptor("gse324984").loader()
    assert set(reloaded.contrasts) == set(original.contrasts)
    for cid, before in original.contrasts.items():
        after = reloaded.contrasts[cid]
        pd.testing.assert_frame_equal(before.table, after.table, check_exact=True)
        assert after.spec == before.spec
        assert after.reference_samples == before.reference_samples
        assert after.treatment_samples == before.treatment_samples
        assert after.details == before.details


def test_reingest_replaces_contrasts_rather_than_duplicating(write_dataset, ingest_script):  # noqa: F811
    write_dataset()
    ingest_script.ingest("gse324984")
    ingest_script.ingest("gse324984")
    with db_session.session_scope() as s:
        n = s.query(ContrastRow).filter_by(dataset_id="gse324984").count()
    assert n == 2


def test_datasets_without_contrasts_load_with_none(write_dataset, ingest_script, monkeypatch):  # noqa: F811
    write_dataset()
    real_load = gse324984.load

    def load_without_contrasts():
        ds = real_load()
        ds.contrasts = {}
        return ds

    # DatasetDescriptor is frozen; replace the registry entry instead.
    descriptor = registry._REGISTRY["gse324984"]
    monkeypatch.setitem(
        registry._REGISTRY,
        "gse324984",
        registry.DatasetDescriptor(**{**descriptor.__dict__, "loader": load_without_contrasts}),
    )
    ingest_script.ingest("gse324984")
    assert registry.get_descriptor("gse324984").loader().contrasts == {}
