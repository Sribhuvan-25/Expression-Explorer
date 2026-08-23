"""
Group comparison — TNMplot's core operation, generalized to any sample
metadata column. Answers PDF comment 1 (MYCN/PP2A across TARGET by ETP
status) and any future "gene X across groups Y" question without new code.
"""
from __future__ import annotations

from itertools import combinations

import pandas as pd
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
        stat, p = stats.mannwhitneyu(v1, v2, alternative="two-sided")
        rows.append({"group_a": g1, "group_b": g2, "n_a": len(v1), "n_b": len(v2), "u_stat": stat, "p_value": p})
    result = pd.DataFrame(rows)
    if not result.empty and len(result) > 1:
        # Benjamini-Hochberg FDR across the pairwise tests shown together
        result = result.sort_values("p_value").reset_index(drop=True)
        m = len(result)
        result["q_value"] = (result["p_value"] * m / (result.index + 1)).clip(upper=1)
        result["q_value"] = result["q_value"][::-1].cummin()[::-1]
    return result


def kruskal_wallis(df: pd.DataFrame) -> dict:
    """Omnibus test across all groups at once, used when there are more
    than two groups (e.g. ETP / near-ETP / non-ETP)."""
    groups = [g["value"].to_numpy() for _, g in df.groupby("group") if len(g) >= 2]
    if len(groups) < 2:
        return {"h_stat": None, "p_value": None}
    h, p = stats.kruskal(*groups)
    return {"h_stat": h, "p_value": p}
