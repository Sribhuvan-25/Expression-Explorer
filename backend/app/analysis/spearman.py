"""
Spearman correlation p-value as R's `cor.test(method = "spearman")`
computes it (METHODS.md 7.6).

scipy's `spearmanr` uses the t-approximation at every n. R, without ties
and for n < 1290, uses AS 89 (Best & Roberts, Appl. Statist. 1975):
exact enumeration for n <= 9 and an Edgeworth series above. The
difference is not cosmetic at this app's sample sizes: on GSE324984
(n = 8) 39 of 200 tie-free genes against MYCL moved across p = 0.05
(MYCL vs YBX2: t-approximation 0.0465, exact 0.0576).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import permutations
from math import exp, sqrt

import numpy as np
from scipy import stats

_EXACT_MAX_N = 9  # AS 89 / R: exact enumeration up to here
_EDGEWORTH = (0.2274, 0.2531, 0.1745, 0.0758, 0.1033, 0.3932, 0.0879, 0.0151, 0.0072, 0.0831, 0.0131, 4.6e-4)


@lru_cache(maxsize=None)
def _null_s(n: int) -> np.ndarray:
    """Sorted S = sum (i - pi_i)^2 over all n! rank permutations."""
    ranks = np.arange(1, n + 1)
    return np.sort(np.fromiter(
        (int(((ranks - np.array(p)) ** 2).sum()) for p in permutations(range(1, n + 1))),
        dtype=np.int64,
    ))


def _prho(n: int, s: float, lower_tail: bool) -> float:
    """AS 89: P[S >= s], or P[S < s] if lower_tail."""
    pv = 0.0 if lower_tail else 1.0
    if s <= 0:
        return pv
    n3 = n * (n * n - 1) / 3
    if s > n3:
        return 1 - pv
    if n <= _EXACT_MAX_N:
        null = _null_s(n)
        n_ge = 1 if s == n3 else int(len(null) - np.searchsorted(null, s, side="left"))
        return (len(null) - n_ge if lower_tail else n_ge) / len(null)
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = _EDGEWORTH
    b = 1 / n
    x = (6 * (s - 1) * b / (n * n - 1) - 1) * sqrt(n - 1)
    y = x * x
    u = x * b * (c1 + b * (c2 + c3 * b) + y * (-c4 + b * (c5 + c6 * b) - y * b * (
        c7 + c8 * b - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))))
    y = u / exp(y / 2)
    pv = (-y if lower_tail else y) + (stats.norm.cdf(x) if lower_tail else stats.norm.sf(x))
    return min(max(pv, 0.0), 1.0)


def spearman_test(x, y) -> tuple[float, float]:
    """(rho, two-sided p), matching R's cor.test default (exact = NULL)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    rho = float(stats.spearmanr(x, y).statistic)
    ties = len(np.unique(x)) < n or len(np.unique(y)) < n
    if n >= 1290 or ties:
        # R falls back to the t-approximation here too (with ties it warns).
        return rho, float(stats.spearmanr(x, y).pvalue)
    q = (n ** 3 - n) * (1 - rho) / 6
    if q > (n ** 3 - n) / 6:
        p = _prho(n, round(q), lower_tail=False)
    else:
        p = _prho(n, round(q) + 2, lower_tail=True)
    return rho, min(2 * p, 1.0)
