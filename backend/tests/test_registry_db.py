"""
Tests for the database-backed registry.

This is the cutover that makes "add a dataset without a redeploy" real:
`list_descriptors()` now merges rows from the `datasets` table over the
code-registered ones. The properties worth pinning are the ones that
make it safe -- an unreachable database must degrade to the previous
behaviour rather than emptying the app, and a database-loaded dataset
must be indistinguishable from a code-loaded one.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from app import registry
from app.db import session as db_session
from app.db.models import Base, DatasetRow, FeatureRow, SampleRow


@pytest.fixture()
def db(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_URL", f"sqlite:///{tmp_path / 'reg.db'}")
    db_session.reset_engine()
    Base.metadata.create_all(db_session.get_engine())
    yield tmp_path
    db_session.reset_engine()


def _seed(tmp_path, dataset_id="dbonly", status="ready"):
    """A dataset that exists ONLY in the database -- no ingest module."""
    matrix = pd.DataFrame(
        {"S1": [1.0, 4.0], "S2": [2.0, 5.0], "S3": [3.0, 6.0]}, index=["GENEA", "GENEB"]
    )
    path = tmp_path / f"{dataset_id}.parquet"
    matrix.to_parquet(path)

    with db_session.session_scope() as s:
        row = DatasetRow(
            dataset_id=dataset_id,
            display_name="DB Only Cohort",
            accession="ACC",
            repository="test",
            assay_type="rna_seq",
            expression_unit="tpm",
            n_samples=3,
            supports_survival=False,
            group_columns=["etp_status"],
            matrix_uri=path.resolve().as_uri(),
            status=status,
            notes="a provenance caveat",
        )
        row.samples = [
            SampleRow(dataset_id=dataset_id, sample_id=f"S{i}",
                      attributes={"etp_status": "ETP" if i == 1 else "non-ETP"})
            for i in (1, 2, 3)
        ]
        row.features = [
            FeatureRow(dataset_id=dataset_id, feature_id=g, symbol=g, aliases=[])
            for g in ("GENEA", "GENEB")
        ]
        s.add(row)
    return path


def test_a_database_only_dataset_is_listed_and_loadable(db):
    """The whole point: this dataset has no ingest module, no register()
    call, and no code change -- only a row and a parquet file."""
    _seed(db)
    ids = {d.dataset_id for d in registry.list_descriptors()}
    assert "dbonly" in ids

    ds = registry.get_descriptor("dbonly").loader()
    assert ds.matrix.shape == (2, 3)
    assert ds.source.accession == "ACC"
    assert ds.source.notes == "a provenance caveat"


def test_group_values_survive_the_round_trip(db):
    """Ingest flattens clinical fields into one JSON blob; the analysis
    layer reads them as columns AND from a nested group_columns dict.
    Both shapes must be restored or grouping silently breaks."""
    _seed(db)
    ds = registry.get_descriptor("dbonly").loader()

    assert list(ds.samples["etp_status"]) == ["ETP", "non-ETP", "non-ETP"]
    assert ds.samples["group_columns"].iloc[0]["etp_status"] == "ETP"


def test_sample_order_follows_the_matrix(db):
    """Sample rows come back from the database in arbitrary order; the
    analysis layer aligns samples against matrix columns, so a drifting
    order would misalign values with their labels."""
    _seed(db)
    ds = registry.get_descriptor("dbonly").loader()
    assert list(ds.samples["sample_id"]) == list(ds.matrix.columns)


def test_only_ready_datasets_are_served(db):
    """A half-ingested dataset must be invisible rather than broken."""
    _seed(db, dataset_id="pending", status="draft")
    assert "pending" not in {d.dataset_id for d in registry.list_descriptors()}


def test_unreachable_database_falls_back_to_code(monkeypatch):
    """The property that makes this cutover safe: the previous behaviour
    is the floor, never the casualty."""
    monkeypatch.setenv("DATABASE_URL", "postgresql+psycopg://nobody@127.0.0.1:9/none")
    db_session.reset_engine()
    try:
        registry.ensure_loaded()
        ids = {d.dataset_id for d in registry.list_descriptors()}
        # The code-registered datasets are still there.
        assert {"depmap", "target_all_p2", "gds4299"} <= ids
    finally:
        db_session.reset_engine()


def test_database_wins_over_code_for_the_same_id(db):
    """Where both exist, the database row is the one an operator most
    recently wrote, so it takes precedence."""
    _seed(db, dataset_id="gds4299")
    descriptor = registry.get_descriptor("gds4299")
    assert descriptor.display_name == "DB Only Cohort"


def test_missing_matrix_file_names_the_path(db):
    """A `ready` row pointing at a deleted file is a real inconsistency;
    the error must say which file so it is fixable."""
    path = _seed(db)
    path.unlink()
    with pytest.raises(FileNotFoundError, match="ingest_dataset"):
        registry.get_descriptor("dbonly").loader()
