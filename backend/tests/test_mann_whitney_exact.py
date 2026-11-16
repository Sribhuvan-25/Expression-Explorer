"""Small-group Mann-Whitney p-values (METHODS.md 1.4a).

The defect: scipy's vectorised `mannwhitneyu(axis=1)` switches every row to
the normal approximation if any row has ties, so at 2 vs 2 a perfectly
separated gene was served p = 0.245 or 0.22 next to a stated floor of
0.333, and GSE324984's genome-wide ranking was ordered by ties.
"""
from math import comb

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from app.analysis.compare import mann_whitney_pvalues, min_p_achievable, pairwise_tests
from app.analysis.differential import top_differential_genes

SIZES = [(2, 2), (2, 3), (3, 5), (5, 3), (4, 6), (5, 2)]


@pytest.mark.parametrize("n_a,n_b", SIZES)
def test_tie_free_rows_equal_scipy_exact(n_a, n_b):
    rng = np.random.default_rng(n_a * 10 + n_b)
    a, b = rng.normal(size=(200, n_a)), rng.normal(0.7, size=(200, n_b))
    _, ours = mann_whitney_pvalues(a, b)
    expected = [stats.mannwhitneyu(x, y, method="exact").pvalue for x, y in zip(a, b)]
    np.testing.assert_allclose(ours, expected, rtol=1e-12)


@pytest.mark.parametrize("n_a,n_b", SIZES)
def test_tied_rows_equal_exact_permutation_and_respect_the_floor(n_a, n_b):
    rng = np.random.default_rng(99 + n_a + n_b)
    # Heavy ties: values from {0, 1, 2}, like zero-inflated expression.
    a, b = rng.integers(0, 3, size=(150, n_a)).astype(float), rng.integers(0, 3, size=(150, n_b)).astype(float)
    _, ours = mann_whitney_pvalues(a, b)

    def u_stat(x, y, axis):
        return stats.mannwhitneyu(x, y, axis=axis, method="asymptotic").statistic

    for row, p in zip(range(len(a)), ours):
        ref = stats.permutation_test(
            (a[row], b[row]), u_stat, permutation_type="independent", vectorized=True,
            n_resamples=np.inf, alternative="two-sided",
        ).pvalue
        assert p == pytest.approx(ref, rel=1e-9, abs=1e-12)
    assert (ours >= min_p_achievable(n_a, n_b) - 1e-12).all()


def test_one_tied_row_does_not_change_the_others():
    """The reviewer's reproducer: a clean, perfectly separated 2 vs 2 row
    must get the exact 1/3 whatever the other rows look like."""
    a = np.array([[1.0, 2.0], [0.0, 0.0]])
    b = np.array([[5.0, 6.0], [0.05, 0.04]])
    _, p = mann_whitney_pvalues(a, b)
    assert p[0] == pytest.approx(1 / 3)
    assert p[1] == pytest.approx(1 / 3)  # one arm all-zero, still exact


def test_large_groups_use_the_stated_approximation():
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=(20, 19)), rng.normal(size=(20, 146))
    _, ours = mann_whitney_pvalues(a, b)
    expected = stats.mannwhitneyu(a, b, axis=1, method="asymptotic").pvalue
    np.testing.assert_allclose(ours, expected, rtol=1e-12)
    assert min_p_achievable(19, 146) is None


def test_differential_never_reports_below_its_floor():
    rng = np.random.default_rng(3)
    genes = [f"G{i}" for i in range(300)]
    values = rng.gamma(2.0, 5.0, size=(300, 4))
    values[:150, :2] = 0.0  # zero arm in half the genes: the tie pattern that broke scipy
    values[0] = [1, 2, 50, 60]  # perfectly separated
    matrix = pd.DataFrame(values, index=genes, columns=["A1", "A2", "B1", "B2"])
    samples = pd.DataFrame({"sample_id": matrix.columns, "group_columns": [{"g": "a"}] * 2 + [{"g": "b"}] * 2})
    res = top_differential_genes(matrix, samples, "g", "a", "b", top_n=300)
    assert res["min_p_achievable"] == pytest.approx(1 / 3)
    assert min(g["p_value"] for g in res["genes"]) >= 1 / 3 - 1e-12


def test_pairwise_tests_use_the_exact_path():
    df = pd.DataFrame({"group": ["x", "x", "y", "y"], "value": [0.0, 0.0, 0.05, 0.04], "sample_id": list("abcd")})
    row = pairwise_tests(df).iloc[0]
    assert row["p_value"] == pytest.approx(1 / 3)
    assert row["p_value"] >= row["min_p_achievable"] - 1e-12


def test_ties_at_the_floor_are_ordered_by_effect_size():
    """At 2 vs 2 every perfectly separated gene has p = 1/3; the documented
    tie-break (|median difference|) must decide the order."""
    genes = ["SMALL", "HUGE", "MID"]
    matrix = pd.DataFrame([[1, 2, 3, 4], [1, 2, 900, 1000], [1, 2, 30, 40]], index=genes,
                          columns=["A1", "A2", "B1", "B2"], dtype=float)
    samples = pd.DataFrame({"sample_id": matrix.columns, "group_columns": [{"g": "a"}] * 2 + [{"g": "b"}] * 2})
    res = top_differential_genes(matrix, samples, "g", "a", "b", top_n=3)
    assert [g["feature_id"] for g in res["genes"]] == ["HUGE", "MID", "SMALL"]
    assert len({g["q_value"] for g in res["genes"]}) == 1  # tied p -> tied q


def test_mixed_exact_and_approximate_pairs_serialise_as_strict_json():
    """DepMap by subtype mixes tiny and large groups in one compare; the
    approximate pairs' None floor became NaN and made /compare 500."""
    import json

    rng = np.random.default_rng(5)
    df = pd.DataFrame({
        "group": ["tiny"] * 2 + ["big1"] * 40 + ["big2"] * 40,  # 2 vs 40: 861 splits, exact
        "value": rng.normal(size=82),
        "sample_id": [f"s{i}" for i in range(82)],
    })
    records = pairwise_tests(df).to_dict("records")
    json.dumps(records, allow_nan=False)
    floors = {(r["group_a"], r["group_b"]): r["min_p_achievable"] for r in records}
    assert floors[("big1", "big2")] is None
    assert floors[("big1", "tiny")] is not None
