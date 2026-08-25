import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from fastapi.encoders import jsonable_encoder
from starlette.responses import JSONResponse

from app.analysis.ranking import _assemble_ranking, rank_by_gene, rank_by_signature


def test_assemble_ranking_missing_metadata_key_serializes_as_null():
    """Reproduces a production 500: a sample missing a group_columns key
    (e.g. TARGET-ALL-P2 samples without an etp_status classification)
    left NaN in the ranking output, and Starlette's default JSONResponse
    rejects NaN outright (allow_nan=False), turning every ranking request
    that touched an incompletely-classified sample into an unhandled 500."""
    samples = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "group_columns": [{"etp_status": "ETP"}, {}, {"etp_status": "near-ETP", "mrd_status": "MRD_neg"}],
        }
    )
    values = pd.Series({"S1": 3.0, "S2": 1.0, "S3": 2.0})

    result = _assemble_ranking(samples, values)

    assert result.loc[result["sample_id"] == "S2", "etp_status"].iloc[0] is None
    assert result.loc[result["sample_id"] == "S1", "etp_status"].iloc[0] == "ETP"

    # The actual regression: this must not raise.
    encoded = jsonable_encoder(result.to_dict(orient="records"))
    JSONResponse(content=encoded)


def test_assemble_ranking_sorts_descending_and_drops_unmatched_samples():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C", "D"],
            "group_columns": [{}, {}, {}, {}],
        }
    )
    # D has no value at all -- must be dropped, not turned into NaN in the output.
    values = pd.Series({"A": 1.0, "B": 3.0, "C": 2.0})

    result = _assemble_ranking(samples, values)

    assert list(result["sample_id"]) == ["B", "C", "A"]
    assert list(result["value"]) == [3.0, 2.0, 1.0]


def test_rank_by_gene_and_signature_end_to_end_with_missing_metadata():
    genes = ["G1", "G2", "G3"]
    matrix = pd.DataFrame(
        {"S1": [5.0, 1.0, 2.0], "S2": [3.0, 4.0, 1.0], "S3": [1.0, 1.0, 1.0]},
        index=genes,
    )
    samples = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "group_columns": [{"etp_status": "ETP"}, {}, {"mrd_status": "MRD_pos"}],
        }
    )

    gene_result = rank_by_gene(matrix, samples, "G1")
    assert list(gene_result["sample_id"]) == ["S1", "S2", "S3"]
    jsonable_encoder(gene_result.to_dict(orient="records"))  # must not raise

    sig_result = rank_by_signature(matrix, samples, ["G1", "G2"], method="log2_mean")
    assert len(sig_result) == 3
    jsonable_encoder(sig_result.to_dict(orient="records"))  # must not raise


def test_assemble_ranking_includes_top_level_group_columns():
    """A dataset may expose a group column as a top-level samples column
    rather than a key inside group_columns -- TARGET-ALL-P2's vital_status
    is a declared field on SampleMetadata, and the dataset advertises it in
    group_columns. Ranking previously flattened only the nested dict, so
    vital_status was silently absent from every ranking row even though
    /compare grouped by it fine."""
    samples = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3"],
            "dataset_id": ["d", "d", "d"],
            "vital_status": ["Alive", "Dead", None],
            "group_columns": [{"etp_status": "ETP"}, {}, {"etp_status": "near-ETP"}],
        }
    )
    values = pd.Series({"S1": 3.0, "S2": 1.0, "S3": 2.0})

    result = _assemble_ranking(samples, values)

    assert "vital_status" in result.columns
    assert result.loc[result["sample_id"] == "S1", "vital_status"].iloc[0] == "Alive"
    assert result.loc[result["sample_id"] == "S2", "vital_status"].iloc[0] == "Dead"
    # Missing values still serialize as null, not NaN.
    assert result.loc[result["sample_id"] == "S3", "vital_status"].iloc[0] is None
    # Nested keys still work alongside top-level ones.
    assert result.loc[result["sample_id"] == "S1", "etp_status"].iloc[0] == "ETP"
    # Bookkeeping columns are not echoed back as metadata.
    assert "dataset_id" not in result.columns
    assert "group_columns" not in result.columns


def test_assemble_ranking_omits_survival_input_columns():
    """days_to_death / days_to_last_follow_up are inputs to the survival
    analysis, not annotation for a ranking table -- including them just
    widens the table with numbers that are meaningless out of that
    context."""
    samples = pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "vital_status": ["Alive", "Dead"],
            "days_to_death": [None, 400.0],
            "days_to_last_follow_up": [1200.0, None],
            "group_columns": [{}, {}],
        }
    )
    result = _assemble_ranking(samples, pd.Series({"S1": 2.0, "S2": 1.0}))

    assert "vital_status" in result.columns
    assert "days_to_death" not in result.columns
    assert "days_to_last_follow_up" not in result.columns

    JSONResponse(content=jsonable_encoder(result.to_dict(orient="records")))


def test_rank_by_gene_missing_gene_raises():
    matrix = pd.DataFrame({"S1": [1.0]}, index=["G1"])
    samples = pd.DataFrame({"sample_id": ["S1"], "group_columns": [{}]})
    with pytest.raises(KeyError):
        rank_by_gene(matrix, samples, "NOTAREALGENE")
