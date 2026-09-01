import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from app.analysis.dimensionality import pca


def _matrix(n_genes=20, n_samples=30, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)]
    data = rng.normal(size=(n_genes, n_samples))
    return pd.DataFrame(data, index=genes, columns=samples)


def _samples_df(sample_ids):
    return pd.DataFrame({"sample_id": sample_ids, "dataset_id": "test", "group_columns": [{} for _ in sample_ids]})


def test_pca_variance_explained_sums_to_at_most_one():
    matrix = _matrix()
    samples = _samples_df(matrix.columns.tolist())
    result = pca(matrix, list(matrix.index), samples, n_components=2)
    assert sum(result["variance_explained"]) <= 1.0 + 1e-9
    assert all(v >= 0 for v in result["variance_explained"])


def test_pca_two_perfectly_separated_clusters_separate_on_pc1():
    """Two groups of samples with a clear expression shift on a handful of
    genes should separate cleanly on the first principal component -- the
    same "do ETP and non-ETP actually structure apart" question this
    feature exists to answer."""
    genes = [f"GENE{i}" for i in range(10)]
    group_a = [f"A{i}" for i in range(15)]
    group_b = [f"B{i}" for i in range(15)]
    rng = np.random.default_rng(1)

    data = {}
    for s in group_a:
        data[s] = rng.normal(loc=0, scale=0.1, size=len(genes))
    for s in group_b:
        data[s] = rng.normal(loc=5, scale=0.1, size=len(genes))
    matrix = pd.DataFrame(data, index=genes)
    samples = _samples_df(matrix.columns.tolist())

    result = pca(matrix, genes, samples, n_components=2)
    pc1_by_sample = {p["sample_id"]: p["pc1"] for p in result["points"]}

    a_scores = [pc1_by_sample[s] for s in group_a]
    b_scores = [pc1_by_sample[s] for s in group_b]
    # The two groups must land on opposite sides of some threshold on PC1
    # -- whichever sign convention this SVD happens to produce, the two
    # groups' PC1 ranges must not overlap.
    assert max(a_scores) < min(b_scores) or max(b_scores) < min(a_scores)
    # With such a strong planted signal, PC1 should explain the vast
    # majority of variance.
    assert result["variance_explained"][0] > 0.8


def test_pca_matches_sklearn_variance_explained():
    """Cross-check against sklearn's PCA (available in this dev/test
    environment, not a runtime dependency of the app) on the same
    standardized data -- confirms the hand-rolled SVD approach isn't
    subtly wrong, not just internally self-consistent."""
    sklearn_decomposition = pytest.importorskip("sklearn.decomposition")

    matrix = _matrix(n_genes=15, n_samples=25, seed=2)
    samples = _samples_df(matrix.columns.tolist())
    genes = list(matrix.index)

    result = pca(matrix, genes, samples, n_components=3)

    sub = matrix.loc[genes].T
    standardized = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=1)
    sk_pca = sklearn_decomposition.PCA(n_components=3)
    sk_pca.fit_transform(standardized.to_numpy())

    np.testing.assert_allclose(
        sorted(result["variance_explained"], reverse=True),
        sorted(sk_pca.explained_variance_ratio_.tolist(), reverse=True),
        atol=1e-8,
    )


def test_pca_drops_missing_genes_rather_than_erroring():
    matrix = _matrix(n_genes=10)
    samples = _samples_df(matrix.columns.tolist())
    genes = list(matrix.index[:5]) + ["NOTAGENE1", "NOTAGENE2"]

    result = pca(matrix, genes, samples, n_components=2)

    assert result["genes_missing"] == ["NOTAGENE1", "NOTAGENE2"]
    assert set(result["genes_used"]) == set(matrix.index[:5])


def test_pca_drops_zero_variance_genes():
    matrix = _matrix(n_genes=10, n_samples=20)
    matrix.loc["CONSTANT"] = 5.0  # identical value for every sample
    samples = _samples_df(matrix.columns.tolist())

    result = pca(matrix, list(matrix.index), samples, n_components=2)

    assert "CONSTANT" in result["genes_zero_variance"]
    assert "CONSTANT" not in result["genes_used"]


def test_pca_too_few_present_genes_raises():
    matrix = _matrix(n_genes=10)
    samples = _samples_df(matrix.columns.tolist())
    with pytest.raises(ValueError, match="at least 2"):
        pca(matrix, ["GENE0", "NOTAGENE"], samples)


def test_pca_n_components_capped_at_number_of_genes():
    matrix = _matrix(n_genes=10, n_samples=20)
    samples = _samples_df(matrix.columns.tolist())
    result = pca(matrix, ["GENE0", "GENE1"], samples, n_components=5)
    assert result["n_components"] == 2
    assert all("pc3" not in p for p in result["points"])


def test_pca_points_carry_group_columns_for_ui_coloring():
    matrix = _matrix(n_genes=10, n_samples=5)
    samples = pd.DataFrame(
        {
            "sample_id": matrix.columns.tolist(),
            "dataset_id": "test",
            "group_columns": [{"lineage": "Lymphoid"} for _ in matrix.columns],
        }
    )
    result = pca(matrix, list(matrix.index), samples, n_components=2)
    assert all(p["group_columns"] == {"lineage": "Lymphoid"} for p in result["points"])
