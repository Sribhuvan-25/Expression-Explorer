import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from app.analysis.correlation import co_expression_network, gene_correlation


def _matrix(n_genes=50, n_samples=40, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)]
    data = rng.normal(size=(n_genes, n_samples))
    return pd.DataFrame(data, index=genes, columns=samples)


def test_gene_correlation_perfect_positive():
    matrix = _matrix()
    matrix.loc["GENE1"] = matrix.loc["GENE0"]  # identical -> r = 1

    result = gene_correlation(matrix, "GENE0", "GENE1", method="pearson")

    assert result["coefficient"] == pytest.approx(1.0)
    assert result["n"] == matrix.shape[1]
    assert len(result["points"]) == matrix.shape[1]


def test_gene_correlation_perfect_negative():
    matrix = _matrix()
    matrix.loc["GENE1"] = -matrix.loc["GENE0"]

    result = gene_correlation(matrix, "GENE0", "GENE1", method="pearson")

    assert result["coefficient"] == pytest.approx(-1.0)


def test_gene_correlation_spearman_matches_scipy_directly():
    from scipy import stats

    matrix = _matrix(seed=1)
    result = gene_correlation(matrix, "GENE2", "GENE5", method="spearman")

    expected_stat, expected_p = stats.spearmanr(matrix.loc["GENE2"], matrix.loc["GENE5"])
    assert result["coefficient"] == pytest.approx(expected_stat)
    assert result["p_value"] == pytest.approx(expected_p)


def test_gene_correlation_unknown_gene_raises_keyerror():
    matrix = _matrix()
    with pytest.raises(KeyError):
        gene_correlation(matrix, "NOTAGENE", "GENE0")


def test_gene_correlation_same_gene_twice_raises_value_error_not_crash():
    """Regression test: gene_a == gene_b used to crash with an unhandled
    exception inside scipy, not a clean error -- matrix.loc[[g, g]].T
    produces two identically-named columns, so paired[gene_a] resolved to
    a 2-column DataFrame instead of a Series, and pearsonr/spearmanr/
    kendalltau all threw on that shape rather than raising a message a
    caller could act on. Found in QA (an easy UI misclick: leaving gene_b
    at its default while typing a new gene_a)."""
    matrix = _matrix()
    for method in ("pearson", "spearman", "kendall"):
        with pytest.raises(ValueError, match="must be different"):
            gene_correlation(matrix, "GENE0", "GENE0", method=method)


def test_gene_correlation_unknown_method_raises_value_error():
    matrix = _matrix()
    with pytest.raises(ValueError, match="Unknown correlation method"):
        gene_correlation(matrix, "GENE0", "GENE1", method="bogus")


def test_gene_correlation_drops_samples_missing_either_gene():
    matrix = _matrix(n_samples=10)
    matrix.loc["GENE0", "S3"] = np.nan
    matrix.loc["GENE1", "S7"] = np.nan

    result = gene_correlation(matrix, "GENE0", "GENE1")

    assert result["n"] == 8  # 10 - 2 samples with a missing value


def test_gene_correlation_too_few_paired_samples_raises():
    matrix = _matrix(n_samples=10)
    matrix.loc["GENE0"] = np.nan
    matrix.loc["GENE0", ["S0", "S1"]] = [1.0, 2.0]  # only 2 paired samples

    with pytest.raises(ValueError, match="at least 3"):
        gene_correlation(matrix, "GENE0", "GENE1")


def test_co_expression_network_finds_the_planted_correlate():
    matrix = _matrix(n_genes=30, n_samples=40, seed=2)
    matrix.loc["GENE29"] = matrix.loc["GENE0"] * 2 + 0.001  # near-perfect positive

    top = co_expression_network(matrix, "GENE0", top_n=5, direction="positive")

    assert top[0]["gene"] == "GENE29"
    assert top[0]["coefficient"] == pytest.approx(1.0, abs=1e-6)
    assert len(top) == 5


def test_co_expression_network_negative_direction_returns_most_anticorrelated():
    matrix = _matrix(n_genes=30, n_samples=40, seed=3)
    matrix.loc["GENE29"] = -matrix.loc["GENE0"] * 2

    top = co_expression_network(matrix, "GENE0", top_n=5, direction="negative")

    assert top[0]["gene"] == "GENE29"
    assert top[0]["coefficient"] == pytest.approx(-1.0, abs=1e-6)
    # Results must actually be sorted ascending (most negative first), not
    # just contain the right gene somewhere in the list.
    coeffs = [t["coefficient"] for t in top]
    assert coeffs == sorted(coeffs)


def test_co_expression_network_excludes_the_query_gene_itself():
    matrix = _matrix(n_genes=20, n_samples=30, seed=4)
    top = co_expression_network(matrix, "GENE0", top_n=25)  # ask for more than exist
    assert all(t["gene"] != "GENE0" for t in top)
    assert len(top) == 19  # every other gene, GENE0 excluded


def test_co_expression_network_rejects_kendall():
    matrix = _matrix()
    with pytest.raises(ValueError, match="Unsupported method"):
        co_expression_network(matrix, "GENE0", method="kendall")


def test_co_expression_network_unknown_gene_raises_keyerror():
    matrix = _matrix()
    with pytest.raises(KeyError):
        co_expression_network(matrix, "NOTAGENE")
