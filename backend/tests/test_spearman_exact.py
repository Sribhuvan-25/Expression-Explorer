"""
Spearman p must match R's cor.test (AS 89: exact for n <= 9, Edgeworth
for n < 1290 without ties). scipy's t-approximation was used at every n;
on GSE324984 (n = 8) MYCL vs YBX2 served p = 0.0465, exact 0.0576.
Values below were computed with R 4.2 cor.test(method = "spearman").
"""
import numpy as np
import pandas as pd
import pytest

from app.analysis.correlation import gene_correlation
from app.analysis.spearman import spearman_test

X8 = [3.1, 0.4, 2.2, 5.9, 4.4, 1.7, 7.3, 6.0]
Y8 = [2.0, 1.1, 0.3, 4.8, 3.9, 2.6, 5.5, 7.7]
X12 = [0.5, 2.9, 1.4, 7.7, 3.3, 6.1, 4.8, 9.2, 8.4, 2.2, 5.6, 10.1]
Y12 = [1.0, 0.7, 3.4, 5.1, 2.2, 8.8, 3.9, 6.6, 10.3, 4.4, 7.2, 9.5]


def test_exact_small_n_matches_r():
    rho, p = spearman_test(X8, Y8)
    assert rho == pytest.approx(R_8[0], abs=1e-12)
    assert p == pytest.approx(R_8[1], rel=1e-12)


def test_edgeworth_matches_r():
    assert spearman_test(X12, Y12)[1] == pytest.approx(R_12, rel=1e-12)


def test_ties_fall_back_to_the_t_approximation_as_r_does():
    x = [1, 2, 2, 3, 4, 5, 6, 7]
    from scipy import stats

    assert spearman_test(x, Y8)[1] == pytest.approx(stats.spearmanr(x, Y8).pvalue)


def test_correlation_endpoint_function_uses_it():
    m = pd.DataFrame([X8, Y8], index=["A", "B"], columns=[f"S{i}" for i in range(8)])
    assert gene_correlation(m, "A", "B", method="spearman")["p_value"] == pytest.approx(R_8[1], rel=1e-12)


R_8 = (0.85714285714285732, 0.010714285714285714)
R_12 = 0.0023689303978236161


def test_no_p_of_exactly_zero_where_the_edgeworth_series_breaks_down():
    """AS 89's Edgeworth series clips to 0 far in the tail (R prints
    p < 2.2e-16); n = 52, rho ~ -0.78 has an asymptotic p near 1e-11."""
    from scipy import stats

    rng = np.random.default_rng(4)
    x = rng.normal(size=52)
    y = -x + rng.normal(scale=0.55, size=52)
    rho, p = spearman_test(x, y)
    assert rho < -0.7
    assert p > 0
    assert p == pytest.approx(stats.spearmanr(x, y).pvalue)
