"""
PCA over a user-supplied gene set, for an unbiased "do these samples
actually separate" check -- e.g. does ETP vs non-ETP structure fall out
of a signature's genes on more than one axis, not just the single-gene
box plots Compare already shows. Modelled on GEPIA3's Dimensionality
Reduction (see reference/gepia3-feature-review.md), identified there as
a cheap win since it's computable from a dataset's existing matrix.

No new dependency: implemented via SVD directly (numpy, already a
requirement) rather than scikit-learn. PCA-via-SVD is a handful of lines,
not a fragile reimplementation of something genuinely complex -- adding a
dependency as heavy as scikit-learn for one function isn't proportionate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def pca(
    matrix: pd.DataFrame,
    gene_ids: list[str],
    samples: pd.DataFrame,
    n_components: int = 2,
) -> dict:
    """PCA of samples over the given gene set. Genes not present in this
    dataset's matrix are silently dropped (same convention as
    auc_signature_score) rather than erroring -- a gene set is often
    reused across datasets on different platforms, and a handful of
    absent genes shouldn't block the other N-k that are present.

    Standardises each gene to zero mean, unit variance before the SVD
    (i.e. this is PCA on the correlation matrix, not the covariance
    matrix) -- without it, a handful of high-variance/high-expression
    genes would dominate every component regardless of whether they
    actually separate samples, which defeats the point of using a curated
    gene set in the first place.
    """
    present = [g for g in gene_ids if g in matrix.index]
    missing = [g for g in gene_ids if g not in matrix.index]
    if len(present) < 2:
        raise ValueError(
            f"Only {len(present)} of {len(gene_ids)} requested genes are present in this dataset -- "
            "need at least 2 to run PCA."
        )
    n_components = min(n_components, len(present))

    sub = matrix.loc[present].T  # rows=samples, cols=genes
    sub = sub.dropna(how="any")
    if len(sub) < 3:
        raise ValueError(f"Only {len(sub)} samples have complete data across the selected genes -- need at least 3.")

    std = sub.std(axis=0, ddof=1)
    zero_variance = std[std == 0].index.tolist()
    # A gene with zero variance across every sample (e.g. undetected
    # throughout, or a constant) contributes nothing to separating
    # samples and, worse, divides by zero under standardisation -- drop
    # it rather than propagating NaN through the whole computation.
    keep = [g for g in present if g not in zero_variance]
    if len(keep) < 2:
        raise ValueError(
            f"Only {len(keep)} of {len(present)} present genes have nonzero variance across these samples -- "
            "need at least 2 to run PCA."
        )
    sub = sub[keep]
    standardized = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=1)

    # PCA via SVD: for mean-centered, unit-variance data X (samples x
    # genes), the principal component scores are U*S and the loadings are
    # V. Dividing S by sqrt(n-1) turns the singular values into the
    # standard per-component scale used for variance-explained.
    u, s, _vt = np.linalg.svd(standardized.to_numpy(), full_matrices=False)
    n = len(standardized)
    explained_variance = (s**2) / (n - 1)
    total_variance = explained_variance.sum()
    variance_ratio = explained_variance / total_variance if total_variance > 0 else explained_variance

    scores = u[:, :n_components] * s[:n_components]

    indexed_samples = samples.set_index("sample_id")
    points = []
    for i, sample_id in enumerate(standardized.index):
        point = {"sample_id": sample_id}
        for c in range(n_components):
            point[f"pc{c + 1}"] = float(scores[i, c])
        if sample_id in indexed_samples.index and "group_columns" in indexed_samples.columns:
            point["group_columns"] = indexed_samples.loc[sample_id, "group_columns"]
        points.append(point)

    return {
        "n_components": n_components,
        "n_samples": n,
        "genes_used": keep,
        "genes_missing": missing,
        "genes_zero_variance": zero_variance,
        "variance_explained": [float(v) for v in variance_ratio[:n_components]],
        "points": points,
    }
