"""
Gene-vs-gene correlation, and the top-N co-expression network derived from
it -- both computable from a dataset's existing expression matrix, no new
data source needed. Modelled on GEPIA3's Correlation Analysis / Expression
Network features (see reference/gepia3-feature-review.md); identified
there as a "cheap win" specifically because nothing here needs an
ingestion pipeline change.
"""
from __future__ import annotations

import pandas as pd
from scipy import stats

_METHODS = {
    "pearson": stats.pearsonr,
    "spearman": stats.spearmanr,
    "kendall": stats.kendalltau,
}


def gene_correlation(
    matrix: pd.DataFrame,
    gene_a: str,
    gene_b: str,
    method: str = "pearson",
) -> dict:
    """Correlation between two genes' expression across every sample in
    this dataset. Samples missing a value for either gene are dropped
    (a gene absent from a dataset's matrix entirely is a caller error --
    checked before this is called, same as every other gene lookup in
    this app -- but a per-sample NaN, if a future dataset ever has one,
    is handled here rather than crashing scipy)."""
    if method not in _METHODS:
        raise ValueError(f"Unknown correlation method: {method!r}")
    if gene_a not in matrix.index:
        raise KeyError(f"{gene_a} not found in this dataset's matrix.")
    if gene_b not in matrix.index:
        raise KeyError(f"{gene_b} not found in this dataset's matrix.")
    if gene_a == gene_b:
        # A gene against itself isn't a meaningful correlation query (r is
        # trivially 1.0) -- reject explicitly rather than let it through:
        # `matrix.loc[[gene_a, gene_b]]` with equal ids produces a frame
        # with two identically-named columns, so `.T` makes `paired[gene_a]`
        # resolve to a 2-column DataFrame instead of a Series, and scipy's
        # correlation functions throw an unhandled exception on that shape
        # (confirmed in QA: an easy UI misclick -- leaving gene_b at its
        # default while typing gene_a -- crashed the endpoint with a raw
        # 500 instead of a clean error).
        raise ValueError(f"Gene A and Gene B must be different (both were {gene_a!r}).")

    paired = matrix.loc[[gene_a, gene_b]].T.dropna()
    n = len(paired)
    if n < 3:
        raise ValueError(
            f"Only {n} samples have values for both genes -- need at least 3 to compute a correlation."
        )

    stat, p_value = _METHODS[method](paired[gene_a], paired[gene_b])
    return {
        "gene_a": gene_a,
        "gene_b": gene_b,
        "method": method,
        "n": n,
        "coefficient": float(stat),
        "p_value": float(p_value),
        "points": [
            {"sample_id": sid, "x": float(row[gene_a]), "y": float(row[gene_b])}
            for sid, row in paired.iterrows()
        ],
    }


def co_expression_network(
    matrix: pd.DataFrame,
    gene: str,
    top_n: int = 20,
    method: str = "pearson",
    direction: str = "positive",
) -> list[dict]:
    """The `top_n` genes most correlated with `gene` across this dataset's
    samples -- GEPIA3's "Expression Network". `direction` picks whether
    "most correlated" means strongest positive or strongest negative
    (anti-correlated) coefficient; it does not mean "top by |r|", since a
    user asking for anti-correlation wants the most negative end
    specifically, not whichever sign happens to be larger in magnitude.

    Computed via matrix-wide correlation (pandas' own `.corrwith`, backed
    by numpy) rather than a per-gene Python loop -- this dataset's matrix
    can be tens of thousands of rows, and a Python-level loop over every
    other gene for a single query would be needlessly slow next to a
    single vectorised call.
    """
    if method not in ("pearson", "spearman"):
        # kendall's O(n log n) implementation still means an O(genes) sweep
        # is one kendall-tau per gene -- prohibitively slow at this matrix
        # size, unlike pearson/spearman which pandas computes vectorised.
        raise ValueError(f"Unsupported method for a network sweep: {method!r} (use 'pearson' or 'spearman')")
    if gene not in matrix.index:
        raise KeyError(f"{gene} not found in this dataset's matrix.")
    if direction not in ("positive", "negative"):
        raise ValueError(f"Unknown direction: {direction!r}")

    others = matrix.drop(index=gene)
    query = matrix.loc[gene]
    corr = others.T.corrwith(query, method=method).dropna()

    ranked = corr.sort_values(ascending=(direction == "negative"))
    top = ranked.head(top_n)
    return [{"gene": g, "coefficient": float(c)} for g, c in top.items()]
