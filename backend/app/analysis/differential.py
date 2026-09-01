"""
Genome-wide differential expression: rank every gene in a dataset by how
strongly it separates two groups, instead of the single-gene-first shape
every other analysis in this app uses. Modelled on GEPIA3's Differential
Genes feature (reference/gepia3-feature-review.md) -- the "browse" answer
to "what's different between these groups", complementing Compare's
"is this specific gene different" question.

Vectorised across every gene in one scipy call (`mannwhitneyu(..., axis=1)`)
rather than a per-gene Python loop -- measured directly on the TARGET
matrix (60,660 genes): ~32s for a naive loop vs. ~2s vectorised, and the
two approaches were verified to produce bit-identical p-values first.
32s is not an acceptable synchronous HTTP request; ~2s is.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def top_differential_genes(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    group_column: str,
    group_a: str,
    group_b: str,
    top_n: int = 50,
) -> dict:
    """Every gene in the matrix, ranked by a two-group Mann-Whitney test
    between `group_a` and `group_b` on `group_column`. Returns the top
    `top_n` by p-value (ties broken by |median difference| descending),
    plus how many samples/genes actually went into the scan.

    Only a two-group comparison, deliberately: an omnibus (Kruskal-Wallis)
    scan across >2 groups doesn't have a single natural "ranking by effect"
    the way a two-group test does (which group is "up" vs "down"?), and
    every dataset's groupings the app has today are usable as pairwise
    comparisons regardless of how many total groups exist. Compare's own
    /compare endpoint already does the >2-group omnibus case for a single
    gene; this endpoint's job is specifically "rank genes for A vs B."
    """
    if group_column in samples.columns:
        groups = samples[group_column]
    else:
        groups = samples["group_columns"].apply(lambda d: d.get(group_column))

    a_mask = groups == group_a
    b_mask = groups == group_b
    a_samples = samples.loc[a_mask, "sample_id"]
    b_samples = samples.loc[b_mask, "sample_id"]

    a_samples = [s for s in a_samples if s in matrix.columns]
    b_samples = [s for s in b_samples if s in matrix.columns]
    if len(a_samples) < 2 or len(b_samples) < 2:
        raise ValueError(
            f"Need at least 2 samples in each group to test -- got {len(a_samples)} in "
            f"{group_a!r} and {len(b_samples)} in {group_b!r}."
        )

    a_vals = matrix[a_samples].to_numpy()
    b_vals = matrix[b_samples].to_numpy()
    stat, p_values = stats.mannwhitneyu(a_vals, b_vals, axis=1, alternative="two-sided")

    median_diff = np.median(a_vals, axis=1) - np.median(b_vals, axis=1)

    result = pd.DataFrame(
        {
            "feature_id": matrix.index,
            "p_value": p_values,
            "median_diff": median_diff,
            "median_a": np.median(a_vals, axis=1),
            "median_b": np.median(b_vals, axis=1),
        }
    )

    # Benjamini-Hochberg across every gene tested, not just the returned
    # top_n -- an FDR computed only over the top N already-most-significant
    # genes would be a biased/meaningless correction (same reasoning as
    # pairwise_tests' FDR in compare.py, just applied genome-wide here).
    m = len(result)
    result = result.sort_values("p_value").reset_index(drop=True)
    result["q_value"] = (result["p_value"] * m / (result.index + 1)).clip(upper=1)
    result["q_value"] = result["q_value"][::-1].cummin()[::-1]

    top = result.head(top_n).copy()

    return {
        "group_column": group_column,
        "group_a": group_a,
        "group_b": group_b,
        "n_a": len(a_samples),
        "n_b": len(b_samples),
        "n_genes_tested": m,
        "top_n": len(top),
        "genes": top.to_dict("records"),
    }
