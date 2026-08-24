"""
Gene-signature scoring, matching the paper's AUCell-based method
(rank genes within each sample, compute area-under-the-curve over the top
fraction of expressed genes). This is what makes gene sets a first-class
input rather than hardcoded — the ETP-TF5, ZMIZ1-5, and any future
MYCN/PP2A signature all flow through this one function.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# np.trapezoid only exists from numpy 2.0 onward (np.trapz, its
# predecessor, is deprecated there but still the only option pre-2.0).
# requirements.txt pins numpy>=2.0, but that doesn't retroactively fix an
# already-deployed environment still running an older numpy underneath a
# transitive pin from pandas/scipy -- confirmed this crashed with
# AttributeError in production (numpy 1.26.4) even with the pin in place
# until the next full rebuild. Falling back at import time means this
# keeps working regardless of which numpy is actually installed anywhere.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def auc_signature_score(
    matrix: pd.DataFrame,
    gene_ids: list[str],
    top_fraction: float = 0.25,
) -> pd.Series:
    """AUCell-equivalent score per sample for an arbitrary gene set.

    For each sample, genes are ranked by expression (rank 1 = highest).
    We consider only the top `top_fraction` of expressed genes (aucMaxRank),
    matching the paper's methods (aucMaxRank = 0.25 for bulk RNA-seq).
    The score is the normalized area under the recovery curve of the
    signature genes within that ranked list — higher means the signature's
    genes are concentrated near the top of this sample's expression ranks.

    Returns one score per sample, NaN for samples where none of the
    requested genes are present in the matrix.
    """
    present = [g for g in gene_ids if g in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns)

    max_rank = max(1, int(round(top_fraction * matrix.shape[0])))
    ranks = matrix.rank(axis=0, ascending=False, method="average")

    scores = {}
    for sample in matrix.columns:
        sample_ranks = ranks.loc[present, sample].sort_values()
        # step function: cumulative count of signature genes found by rank r,
        # truncated at max_rank; AUC = area under that step function
        within = sample_ranks[sample_ranks <= max_rank]
        if within.empty:
            scores[sample] = 0.0
            continue
        hit_ranks = within.sort_values().to_numpy()
        cum_hits = np.arange(1, len(hit_ranks) + 1)
        # trapezoid from rank 0 to max_rank, step increments at each hit
        x = np.concatenate([[0], hit_ranks, [max_rank]])
        y = np.concatenate([[0], cum_hits, [cum_hits[-1]]])
        auc = _trapezoid(y, x)
        max_possible = max_rank * len(present)
        scores[sample] = auc / max_possible if max_possible > 0 else 0.0

    return pd.Series(scores, name="signature_score")


def log2_mean_signature_score(matrix: pd.DataFrame, gene_ids: list[str]) -> pd.Series:
    """Simpler alternative used in a few of the paper's figures (e.g. Fig.
    10A/B): mean of log2-transformed expression across signature genes.
    Offered alongside AUC scoring since the paper switches between both
    depending on the figure."""
    present = [g for g in gene_ids if g in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns)
    log_vals = np.log2(matrix.loc[present] + 1)
    return log_vals.mean(axis=0).rename("signature_score")
