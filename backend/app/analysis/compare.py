"""
Group comparison — TNMplot's core operation, generalized to any sample
metadata column. Answers PDF comment 1 (MYCN/PP2A across TARGET by ETP
status) and any future "gene X across groups Y" question without new code.
"""
from __future__ import annotations

from itertools import combinations
from math import comb

import pandas as pd
import numpy as np
from scipy import stats


def expression_by_group(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    gene_id: str,
    group_column: str,
) -> pd.DataFrame:
    """Long-format frame of (sample_id, group, value) for one gene, ready
    to hand to a boxplot. `group_column` must be a key present in each
    sample's `group_columns` dict, or a top-level samples column."""
    if gene_id not in matrix.index:
        raise KeyError(f"{gene_id} not found in this dataset's matrix.")

    values = matrix.loc[gene_id]
    if group_column in samples.columns:
        groups = samples[group_column]
    else:
        groups = samples["group_columns"].apply(lambda d: d.get(group_column))

    df = pd.DataFrame({"sample_id": samples["sample_id"], "group": groups, "value": values.reindex(samples["sample_id"]).to_numpy()})
    return df.dropna(subset=["group", "value"])


# Largest number of group splits enumerated for an exact test. 5 vs 3 is
# 56, 7 vs 7 is 3,432; 8 vs 8 (12,870) and anything bigger use the normal
# approximation, which is accurate at those sizes.
EXACT_SPLIT_LIMIT = 5000


def min_p_achievable(n_a: int, n_b: int) -> float | None:
    """Smallest two-sided p-value the test used for these group sizes can
    return: 2 / C(n_a + n_b, n_a) when the exact test is used (complete
    separation), None when the normal approximation is used (no fixed
    floor worth stating at those sizes).

    Reported next to every p-value because at small n it is the whole
    story: 2 vs 2 samples cannot get below 0.333, 2 vs 5 not below 0.095,
    so "not significant" there says nothing about the data (METHODS.md
    1.4a). The bound genuinely holds for every served p: small groups get
    the exact permutation test below, ties included."""
    if comb(n_a + n_b, n_a) > EXACT_SPLIT_LIMIT:
        return None
    return float(min(1.0, 2.0 / comb(n_a + n_b, n_a)))


def mann_whitney_pvalues(a_vals, b_vals):
    """Two-sided Mann-Whitney p-values, one per row of `a_vals`/`b_vals`
    (genes x samples). Returns (u_statistic, p_value) arrays.

    Small groups (<= EXACT_SPLIT_LIMIT possible splits): the exact
    permutation distribution of U, enumerating every split of the pooled
    samples, with average ranks for ties -- p = min(1, 2 x the smaller
    tail), scipy's exact-test convention. For tie-free rows this is
    identical to `mannwhitneyu(method="exact")`; for tied rows it is the
    exact conditional test, where scipy would silently switch to the
    normal approximation.

    Why not scipy: with `axis=1`, one tied row anywhere in the matrix
    switches *every* row to the normal approximation, which at 2 vs 2
    served p = 0.22 for a gene whose exact p cannot go below 0.333, and
    ranked GSE324984's genome-wide list by ties rather than effect
    (found in review, METHODS.md 1.4a).

    Larger groups: the normal approximation with tie and continuity
    correction (`mannwhitneyu(method="asymptotic")`), stated explicitly
    rather than left to `method="auto"`.
    """
    a_vals = np.atleast_2d(np.asarray(a_vals, dtype=float))
    b_vals = np.atleast_2d(np.asarray(b_vals, dtype=float))
    n_a, n_b = a_vals.shape[1], b_vals.shape[1]
    if comb(n_a + n_b, n_a) > EXACT_SPLIT_LIMIT:
        u, p = stats.mannwhitneyu(a_vals, b_vals, axis=1, alternative="two-sided", method="asymptotic")
        return np.asarray(u, dtype=float), np.asarray(p, dtype=float)
    pooled = np.concatenate([a_vals, b_vals], axis=1)
    ranks = stats.rankdata(pooled, axis=1)  # average ranks for ties
    offset = n_a * (n_a + 1) / 2
    u_obs = ranks[:, :n_a].sum(axis=1) - offset
    n_le = np.zeros(len(pooled))
    n_ge = np.zeros(len(pooled))
    n_splits = 0
    eps = 1e-9
    for idx in combinations(range(n_a + n_b), n_a):
        u = ranks[:, list(idx)].sum(axis=1) - offset
        n_le += u <= u_obs + eps
        n_ge += u >= u_obs - eps
        n_splits += 1
    p = np.minimum(1.0, 2 * np.minimum(n_le, n_ge) / n_splits)
    return u_obs, p


def pairwise_tests(df: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U between every pair of groups, matching the paper's
    stated default (two-sided Wilcoxon/Mann-Whitney unless noted otherwise
    in a figure legend)."""
    groups = sorted(df["group"].unique())
    rows = []
    for g1, g2 in combinations(groups, 2):
        v1 = df.loc[df["group"] == g1, "value"]
        v2 = df.loc[df["group"] == g2, "value"]
        if len(v1) < 2 or len(v2) < 2:
            continue
        u, pv = mann_whitney_pvalues(v1.to_numpy()[None, :], v2.to_numpy()[None, :])
        stat, p = float(u[0]), float(pv[0])
        rows.append({
            "group_a": g1, "group_b": g2, "n_a": len(v1), "n_b": len(v2), "u_stat": stat, "p_value": p,
            "min_p_achievable": min_p_achievable(len(v1), len(v2)),
        })
    result = pd.DataFrame(rows)
    if not result.empty and len(result) > 1:
        # Benjamini-Hochberg FDR across the pairwise tests shown together
        result = result.sort_values("p_value").reset_index(drop=True)
        m = len(result)
        result["q_value"] = (result["p_value"] * m / (result.index + 1)).clip(upper=1)
        result["q_value"] = result["q_value"][::-1].cummin()[::-1]
    if not result.empty:
        # A mix of exact pairs (a number) and approximate pairs (None) in one
        # frame makes pandas store None as NaN, which is not valid JSON and
        # 500'd /compare on DepMap by subtype. Keep None as None.
        floor = result["min_p_achievable"].astype(object)
        result["min_p_achievable"] = floor.where(floor.notna(), None)
    return result


def kruskal_wallis(df: pd.DataFrame) -> dict:
    """Omnibus test across all groups at once, used when there are more
    than two groups (e.g. ETP / near-ETP / non-ETP)."""
    groups = [g["value"].to_numpy() for _, g in df.groupby("group") if len(g) >= 2]
    if len(groups) < 2:
        return {"h_stat": None, "p_value": None}
    # Every value identical (an unexpressed gene: all-zero TPM): H is 0/0,
    # scipy raises and /compare returned 500 for 47% of GSE324984's genes
    # and 3,169 of TARGET's. R's kruskal.test gives NA; so do we, saying why.
    if len(np.unique(np.concatenate(groups))) == 1:
        return {
            "h_stat": None,
            "p_value": None,
            "note": "Every value is identical, so there is nothing to rank and no test is possible.",
        }
    h, p = stats.kruskal(*groups)
    return {"h_stat": h, "p_value": p}
