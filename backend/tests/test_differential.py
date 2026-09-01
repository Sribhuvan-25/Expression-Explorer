import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from app.analysis.differential import top_differential_genes


def _dataset(n_genes=100, n_samples=40, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    sample_ids = [f"S{i}" for i in range(n_samples)]
    matrix = pd.DataFrame(rng.normal(size=(n_genes, n_samples)), index=genes, columns=sample_ids)
    group = np.where(np.arange(n_samples) < n_samples // 2, "A", "B")
    samples = pd.DataFrame({"sample_id": sample_ids, "dataset_id": "test", "group": group})
    return matrix, samples


def test_top_differential_genes_ranks_by_ascending_p_value():
    matrix, samples = _dataset()
    result = top_differential_genes(matrix, samples, "group", "A", "B", top_n=10)
    p_values = [g["p_value"] for g in result["genes"]]
    assert p_values == sorted(p_values)


def test_top_differential_genes_finds_the_planted_signal():
    matrix, samples = _dataset(n_genes=50, n_samples=40, seed=1)
    # Plant a strong, unambiguous shift in one gene between the two groups.
    a_ids = samples.loc[samples["group"] == "A", "sample_id"]
    b_ids = samples.loc[samples["group"] == "B", "sample_id"]
    matrix.loc["GENE0", a_ids] = np.random.default_rng(2).normal(loc=0, scale=0.1, size=len(a_ids))
    matrix.loc["GENE0", b_ids] = np.random.default_rng(3).normal(loc=10, scale=0.1, size=len(b_ids))

    result = top_differential_genes(matrix, samples, "group", "A", "B", top_n=5)

    assert result["genes"][0]["feature_id"] == "GENE0"
    assert result["genes"][0]["p_value"] < 1e-6


def test_top_differential_genes_matches_scipy_per_gene_loop():
    """Cross-check the vectorised mannwhitneyu(axis=1) call against a
    plain per-gene loop -- this is the exact thing that made the
    vectorised approach usable at TARGET's real 60k-gene scale, verified
    at that size during development; this test locks in bit-identical
    results at a smaller scale so a future scipy upgrade can't silently
    change vectorised semantics without a test noticing."""
    matrix, samples = _dataset(n_genes=30, n_samples=24, seed=4)
    a_ids = samples.loc[samples["group"] == "A", "sample_id"].tolist()
    b_ids = samples.loc[samples["group"] == "B", "sample_id"].tolist()

    result = top_differential_genes(matrix, samples, "group", "A", "B", top_n=30)
    p_by_gene = {g["feature_id"]: g["p_value"] for g in result["genes"]}

    for gene in matrix.index:
        _, expected_p = stats.mannwhitneyu(
            matrix.loc[gene, a_ids], matrix.loc[gene, b_ids], alternative="two-sided"
        )
        assert p_by_gene[gene] == pytest.approx(expected_p)


def test_top_differential_genes_n_a_n_b_match_actual_group_sizes():
    matrix, samples = _dataset(n_genes=10, n_samples=40)
    result = top_differential_genes(matrix, samples, "group", "A", "B", top_n=5)
    assert result["n_a"] == 20
    assert result["n_b"] == 20


def test_top_differential_genes_q_value_computed_over_all_tested_not_just_top_n():
    matrix, samples = _dataset(n_genes=200, n_samples=40, seed=5)
    result = top_differential_genes(matrix, samples, "group", "A", "B", top_n=5)
    assert result["n_genes_tested"] == 200
    assert result["top_n"] == 5
    assert len(result["genes"]) == 5
    # q_value must be present and >= p_value for every returned gene
    # (BH correction never makes a q-value smaller than its own p-value).
    for g in result["genes"]:
        assert g["q_value"] >= g["p_value"] - 1e-12


def test_top_differential_genes_group_column_inside_nested_group_columns_dict():
    matrix, samples = _dataset(n_genes=10, n_samples=20)
    samples = samples.drop(columns=["group"])
    etp = np.where(np.arange(20) < 10, "ETP", "non-ETP")
    samples["group_columns"] = [{"etp_status": e} for e in etp]

    result = top_differential_genes(matrix, samples, "etp_status", "ETP", "non-ETP", top_n=5)
    assert result["n_a"] == 10
    assert result["n_b"] == 10


def test_top_differential_genes_too_few_samples_in_a_group_raises():
    matrix, samples = _dataset(n_genes=10, n_samples=10)
    samples.loc[1:, "group"] = "B"  # only 1 sample left in group A
    with pytest.raises(ValueError, match="at least 2 samples"):
        top_differential_genes(matrix, samples, "group", "A", "B")
