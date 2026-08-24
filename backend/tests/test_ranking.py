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


def test_rank_by_gene_missing_gene_raises():
    matrix = pd.DataFrame({"S1": [1.0]}, index=["G1"])
    samples = pd.DataFrame({"sample_id": ["S1"], "group_columns": [{}]})
    with pytest.raises(KeyError):
        rank_by_gene(matrix, samples, "NOTAREALGENE")
