"""
One gene must count once, and one symbol must resolve to one sensible row.

- TARGET has 110 symbols over several Ensembl rows; the first-listed row
  was used, and for RGS5 / DGCR5 that is a near-silent or all-zero row.
- A gene given twice (a symbol and its feature id) entered PCA k^2 times
  and was weighted double in the log2-mean score.
"""
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from app import registry
from app.analysis.dimensionality import pca
from app.analysis.signature import log2_mean_signature_score
from app.api import main as api_main
from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

DATASET_ID = "duplicate_symbol_test"
N = 12


def _build():
    rng = np.random.default_rng(2)
    ids = ["ENSG_DUP_LOW", "ENSG_DUP_HIGH", "ENSG_A", "ENSG_B", "ENSG_C"]
    symbols = ["DUP", "DUP", "AAA", "BBB", "CCC"]
    matrix = pd.DataFrame(rng.uniform(1, 10, size=(5, N)), index=ids, columns=[f"S{i}" for i in range(N)])
    matrix.loc["ENSG_DUP_LOW"] = rng.uniform(0, 0.2, N)  # listed first, barely expressed
    matrix.loc["ENSG_DUP_HIGH"] = rng.uniform(20, 40, N)
    samples = pd.DataFrame(
        {"sample_id": matrix.columns, "dataset_id": DATASET_ID, "group_columns": [{} for _ in range(N)]}
    )
    features = pd.DataFrame({"feature_id": ids, "symbol": symbols, "aliases": [[] for _ in ids]}, index=ids)
    source = DatasetSource(
        dataset_id=DATASET_ID, display_name="dup", accession="none", repository="test",
        assay_type=AssayType.RNA_SEQ, expression_unit=ExpressionUnit.TPM, n_samples=N,
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


@pytest.fixture(scope="module")
def client():
    if DATASET_ID not in {d.dataset_id for d in registry.list_descriptors()}:
        registry.register(registry.DatasetDescriptor(dataset_id=DATASET_ID, display_name="dup", loader=_build))
    api_main._get_dataset.cache_clear()
    return TestClient(api_main.app)


def test_a_duplicated_symbol_resolves_to_its_most_expressed_row(client):
    ds = api_main._get_dataset(DATASET_ID)
    assert api_main._resolve_gene(ds, "DUP") == "ENSG_DUP_HIGH"
    body = client.get(f"/datasets/{DATASET_ID}/correlation", params={"gene_a": "DUP", "gene_b": "AAA"}).json()
    expected = np.log2(ds.matrix.loc["ENSG_DUP_HIGH"].to_numpy() + 1)  # TPM, on the log2(x+1) scale
    assert sorted(p["x"] for p in body["points"]) == pytest.approx(sorted(expected))


def test_pca_counts_a_repeated_gene_once(client):
    once = client.post(f"/datasets/{DATASET_ID}/pca", json={"genes": ["AAA", "BBB", "CCC"]}).json()
    twice = client.post(f"/datasets/{DATASET_ID}/pca", json={"genes": ["AAA", "ENSG_A", "aaa", "BBB", "CCC"]}).json()
    assert twice["variance_explained"] == pytest.approx(once["variance_explained"])
    assert twice["genes_used"] == ["AAA", "BBB", "CCC"]  # symbols, each once


def test_analysis_functions_dedupe_too():
    ds = _build()
    m = ds.matrix
    assert pca(m, ["ENSG_A", "ENSG_A", "ENSG_B"], ds.samples)["variance_explained"] == pytest.approx(
        pca(m, ["ENSG_A", "ENSG_B"], ds.samples)["variance_explained"]
    )
    pd.testing.assert_series_equal(
        log2_mean_signature_score(m, ["ENSG_A", "ENSG_A", "ENSG_B"]),
        log2_mean_signature_score(m, ["ENSG_A", "ENSG_B"]),
    )
