"""
Gene-signature scoring, matching the paper's AUCell-based method
(rank genes within each sample, compute AUCell's area under the recovery
curve over the top fraction of expressed genes). This is what makes gene
sets a first-class input rather than hardcoded — the ETP-TF5, ZMIZ1-5, and any future
MYCN/PP2A signature all flow through this one function.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def auc_signature_score(
    matrix: pd.DataFrame,
    gene_ids: list[str],
    top_fraction: float = 0.25,
) -> pd.Series:
    """AUCell score per sample for an arbitrary gene set (METHODS.md 2.1-2.3).

    Genes are ranked within each sample (rank 1 = highest expression) and
    the score is AUCell's normalised area under the recovery curve over
    ranks below `aucThreshold = round(top_fraction * n_genes)` -- the
    paper's aucMaxRank of 25% for bulk RNA-seq. As in AUCell's `.auc`:

        AUC    = sum(diff(c(x, aucThreshold)) * seq_along(x)),
                 x = sorted ranks of signature genes, x < aucThreshold
        maxAUC = the same area with the signature genes at ranks 1..n
        score  = AUC / maxAUC                          (normAUC = TRUE)

    The step-function area telescopes to sum(aucThreshold - x), so each
    signature gene ranked r contributes max(aucThreshold - r, 0).

    Ties: AUCell breaks them at random, so its score varies between runs.
    A gene tied over ranks a..b gets the mean of its contribution over
    that block -- the expected AUCell score, exactly, and deterministic.
    With no ties this is AUCell's score.

    Returns one score per sample in [0, 1]; NaN for every sample when
    none of the requested genes are in the matrix. Duplicated genes count
    once, as in AUCell.
    """
    present = list(dict.fromkeys(g for g in gene_ids if g in matrix.index))
    if not present:
        return pd.Series(np.nan, index=matrix.columns)

    threshold = max(1, int(round(top_fraction * matrix.shape[0])))
    # A tie block's ranks: a = 1 + (# strictly higher), b = # higher or equal.
    neg_sorted = np.sort(-matrix.to_numpy(dtype=float), axis=0)
    neg_sig = -matrix.loc[present].to_numpy(dtype=float)
    a = np.empty_like(neg_sig)
    b = np.empty_like(neg_sig)
    for j in range(matrix.shape[1]):
        a[:, j] = np.searchsorted(neg_sorted[:, j], neg_sig[:, j], side="left") + 1
        b[:, j] = np.searchsorted(neg_sorted[:, j], neg_sig[:, j], side="right")
    # Mean over r in a..b of max(threshold - r, 0): sum the ranks below
    # the threshold, divide by the block size.
    hi = np.minimum(b, threshold - 1)
    count = np.clip(hi - a + 1, 0, None)
    contribution = count * threshold - count * (a + hi) / 2
    contribution = np.where(count > 0, contribution / (b - a + 1), 0.0)
    contribution[np.isnan(neg_sig)] = 0.0
    auc = contribution.sum(axis=0)

    best = np.arange(1, min(len(present), threshold - 1) + 1)
    max_auc = float(np.sum(threshold - best))
    scores = auc / max_auc if max_auc > 0 else np.zeros(matrix.shape[1])
    return pd.Series(scores, index=matrix.columns, name="signature_score")


LOG2_UNITS = {"log2_tpm", "log2_intensity"}


def is_log2_unit(unit) -> bool:
    """True for expression units already on a log2 scale."""
    return getattr(unit, "value", unit) in LOG2_UNITS


def log2_mean_signature_score(
    matrix: pd.DataFrame, gene_ids: list[str], already_log2: bool = False
) -> pd.Series:
    """Simpler alternative used in a few of the paper's figures (e.g. Fig.
    10A/B): mean of log2-transformed expression across signature genes.
    Offered alongside AUC scoring since the paper switches between both
    depending on the figure.

    `already_log2`: the matrix is already log2 (DepMap log2(TPM+1),
    microarray log2 intensity) and is averaged as is. Transforming it
    again gave the mean of log2(log2(x)+1) on three of five datasets
    (METHODS.md 2.4)."""
    present = list(dict.fromkeys(g for g in gene_ids if g in matrix.index))  # each gene once, as AUC
    if not present:
        return pd.Series(np.nan, index=matrix.columns)
    log_vals = matrix.loc[present] if already_log2 else np.log2(matrix.loc[present] + 1)
    return log_vals.mean(axis=0).rename("signature_score")
