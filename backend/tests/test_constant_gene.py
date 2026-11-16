"""
A gene with the same value in every sample (unexpressed, all-zero TPM)
has no variance. Correlation and expression-vs-aux returned NaN, which is
not valid JSON, so the endpoints answered 500; co-expression returned an
empty network with no reason. 30,809 of GSE324984's genes and 2,445 of
TARGET's are constant.
"""
import json

import numpy as np
import pandas as pd
import pytest

from app.analysis.aux_analysis import expression_vs_aux_correlation
from app.analysis.correlation import co_expression_network, gene_correlation
from tests.test_api import DATASET_ID, client  # noqa: F401  (fixture)
from app.api import main as api_main

MATRIX = pd.DataFrame(
    {f"S{i}": [0.0, float(i), float(i * i % 7)] for i in range(8)}, index=["FLAT", "UP", "OTHER"]
)


@pytest.mark.parametrize("method", ["pearson", "spearman", "kendall"])
def test_correlation_with_a_constant_gene_is_a_clear_error(method):
    with pytest.raises(ValueError, match="same value .* in every sample"):
        gene_correlation(MATRIX, "FLAT", "UP", method=method)
    with pytest.raises(ValueError, match="FLAT has the same value"):
        gene_correlation(MATRIX, "UP", "FLAT", method=method)
    json.dumps(gene_correlation(MATRIX, "UP", "OTHER", method=method), allow_nan=False)


def test_co_expression_of_a_constant_gene_says_why():
    with pytest.raises(ValueError, match="no variance"):
        co_expression_network(MATRIX, "FLAT")


def test_expression_vs_aux_with_a_constant_gene_is_a_422(client):  # noqa: F811
    ds = api_main._get_dataset(DATASET_ID)
    original = ds.matrix.loc["MYCN"].copy()
    try:
        ds.matrix.loc["MYCN"] = 3.0
        r = client.get(
            f"/datasets/{DATASET_ID}/expression-vs-aux",
            params={"gene": "MYCN", "layer": "crispr_gene_effect", "aux_feature": "DENSE"},
        )
        assert r.status_code == 422, r.text
        assert "no variance" in r.json()["detail"]
    finally:
        ds.matrix.loc["MYCN"] = original
