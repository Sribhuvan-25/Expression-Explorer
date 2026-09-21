"""
Tests for the dataset metadata schema and the ingest command.

The point of this schema is that adding a dataset stops being a code
change. These pin the properties that make that true: a dataset's
description round-trips through the database, cohorts with genuinely
different clinical fields coexist, and a failed ingest stays visibly
failed rather than silently absent.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest
from sqlalchemy import func, select

from app.db import session as db_session
from app.db.models import AuxLayerRow, Base, DatasetRow, FeatureRow, SampleRow


@pytest.fixture()
def db(tmp_path, monkeypatch):
    """A throwaway SQLite database per test."""
    monkeypatch.setenv("DATABASE_URL", f"sqlite:///{tmp_path / 'test.db'}")
    db_session.reset_engine()
    Base.metadata.create_all(db_session.get_engine())
    yield db_session
    db_session.reset_engine()


def _dataset_row(dataset_id="ds1", **overrides):
    defaults = dict(
        dataset_id=dataset_id,
        display_name="Test Cohort",
        accession="ACC1",
        repository="test",
        assay_type="rna_seq",
        expression_unit="tpm",
        n_samples=2,
        supports_survival=True,
        group_columns=["etp_status"],
        matrix_uri="file:///tmp/x.parquet",
        status="ready",
    )
    defaults.update(overrides)
    return DatasetRow(**defaults)


def test_dataset_round_trips(db):
    with db.session_scope() as s:
        s.add(_dataset_row(notes="a provenance caveat"))

    with db.session_scope() as s:
        row = s.get(DatasetRow, "ds1")
        assert row.display_name == "Test Cohort"
        assert row.group_columns == ["etp_status"]
        assert row.notes == "a provenance caveat"
        assert row.matrix_uri.endswith(".parquet")


def test_cohorts_with_different_clinical_fields_coexist(db):
    """The reason `attributes` is JSON rather than typed columns: TARGET
    carries vital_status/days_to_death, DepMap carries lineage/subtype,
    GDS4299 carries only etp_status. Typed columns would mean a schema
    migration per cohort, defeating the point of the whole change."""
    with db.session_scope() as s:
        target = _dataset_row("target")
        target.samples = [
            SampleRow(
                dataset_id="target",
                sample_id="P1",
                attributes={"etp_status": "ETP", "vital_status": "Dead", "days_to_death": 412},
            )
        ]
        depmap = _dataset_row("depmap", expression_unit="log2_tpm")
        depmap.samples = [
            SampleRow(
                dataset_id="depmap",
                sample_id="ACH-1",
                attributes={"lineage": "Lymphoid", "subtype": "T-ALL"},
            )
        ]
        s.add_all([target, depmap])

    with db.session_scope() as s:
        p1 = s.scalar(select(SampleRow).where(SampleRow.sample_id == "P1"))
        ach = s.scalar(select(SampleRow).where(SampleRow.sample_id == "ACH-1"))
        assert p1.attributes["days_to_death"] == 412
        assert ach.attributes["lineage"] == "Lymphoid"
        # Neither cohort needed the other's columns to exist.
        assert "lineage" not in p1.attributes
        assert "vital_status" not in ach.attributes


def test_can_query_inside_attributes(db):
    """Filtering by a clinical field must still work despite the blob --
    otherwise JSON would have traded away the reason to use a database."""
    with db.session_scope() as s:
        row = _dataset_row("ds1")
        row.samples = [
            SampleRow(dataset_id="ds1", sample_id=f"S{i}",
                      attributes={"etp_status": "ETP" if i < 3 else "non-ETP"})
            for i in range(10)
        ]
        s.add(row)

    with db.session_scope() as s:
        n = s.scalar(
            select(func.count()).select_from(SampleRow).where(
                SampleRow.attributes["etp_status"].as_string() == "ETP"
            )
        )
        assert n == 3


def test_feature_lookup_spans_datasets_without_loading_matrices(db):
    """'Which datasets have MEF2C' was previously unanswerable without
    reading every parquet file into memory."""
    with db.session_scope() as s:
        a = _dataset_row("a")
        a.features = [FeatureRow(dataset_id="a", feature_id="ENSG001", symbol="MEF2C", aliases=[])]
        b = _dataset_row("b")
        b.features = [FeatureRow(dataset_id="b", feature_id="MEF2C", symbol="MEF2C", aliases=[])]
        c = _dataset_row("c")
        c.features = [FeatureRow(dataset_id="c", feature_id="OTHER", symbol="OTHER", aliases=[])]
        s.add_all([a, b, c])

    with db.session_scope() as s:
        hits = s.execute(
            select(FeatureRow.dataset_id, FeatureRow.feature_id).where(FeatureRow.symbol == "MEF2C")
        ).all()
        assert {h[0] for h in hits} == {"a", "b"}
        # The dataset's own id space is preserved -- the matrix is indexed
        # by this exact value, so normalising it would break the join.
        assert dict(hits)["a"] == "ENSG001"


def test_deleting_a_dataset_removes_its_children(db):
    """Re-ingest replaces wholesale; orphaned sample rows from an older
    release must not survive to mix with a newer matrix."""
    with db.session_scope() as s:
        row = _dataset_row("ds1")
        row.samples = [SampleRow(dataset_id="ds1", sample_id="S1", attributes={})]
        row.features = [FeatureRow(dataset_id="ds1", feature_id="F1", symbol="F1", aliases=[])]
        row.aux_layers = [
            AuxLayerRow(
                dataset_id="ds1", layer="crispr_gene_effect", value_label="x",
                value_description="x", source_note="x", n_features=1,
                n_usable_features=1, n_samples_covered=1, matrix_uri="file:///tmp/a.parquet",
            )
        ]
        s.add(row)

    with db.session_scope() as s:
        s.delete(s.get(DatasetRow, "ds1"))

    with db.session_scope() as s:
        assert s.scalar(select(func.count()).select_from(SampleRow)) == 0
        assert s.scalar(select(func.count()).select_from(FeatureRow)) == 0
        assert s.scalar(select(func.count()).select_from(AuxLayerRow)) == 0


def test_a_failed_ingest_stays_visible(db):
    """A dataset whose ingest failed should be recorded as failed, not
    silently missing -- an absent dataset looks like it was never
    requested, which hides the problem."""
    with db.session_scope() as s:
        s.add(_dataset_row("broken", status="failed", status_detail="GDC unreachable"))

    with db.session_scope() as s:
        row = s.get(DatasetRow, "broken")
        assert row.status == "failed"
        assert "GDC unreachable" in row.status_detail


def test_sample_attributes_flattening_merges_both_cohort_shapes():
    """TARGET keeps clinical fields as real DataFrame columns; DepMap
    nests them in a per-sample `group_columns` dict. Ingest must capture
    both, and drop NaN rather than writing it (not JSON-serialisable, and
    it means 'absent')."""
    from scripts.ingest_dataset import _sample_attributes

    class _DS:
        samples = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "dataset_id": ["d", "d"],
                "vital_status": ["Dead", None],
                "days_to_death": [412.0, float("nan")],
                "group_columns": [{"etp_status": "ETP"}, {"etp_status": "non-ETP"}],
            }
        )

    rows = _sample_attributes(_DS())
    assert rows[0]["attributes"] == {
        "vital_status": "Dead",
        "days_to_death": 412.0,
        "etp_status": "ETP",
    }
    # S2 had no vital_status and a NaN follow-up -- both absent, not null.
    assert rows[1]["attributes"] == {"etp_status": "non-ETP"}
