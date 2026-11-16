"""
Two-group differential expression for small replicated perturbation
experiments (e.g. GSE324984: 2 control vs 2 MYCL-overexpressing C4-2B
samples), computed on raw counts with the edgeR + limma-trend pipeline.

Why not the app's usual Mann-Whitney (METHODS.md 1.1)? With n=2 per
group the smallest two-sided Mann-Whitney p-value that *can exist* is
2/C(4,2) = 0.333 -- every gene in the genome is "not significant" by
construction, whatever the data say. A replicated knockdown/overexpression
experiment needs a test that borrows variance information across genes,
which is exactly what limma's empirical-Bayes moderated t does. See
METHODS.md 11.2 for the decision and its alternatives.

This is a faithful re-implementation, not an approximation: each step
mirrors the named R function with that function's default arguments, and
`tests/test_contrast.py` checks the output against values produced by R
(edgeR 3.40 / limma 3.54) on the real GSE324984 counts.

    filterByExpr(y, group)                  -> filter_by_expr()
    calcNormFactors(y, method = "TMM")      -> tmm_norm_factors()
    cpm(y, log = TRUE, prior.count = 3)     -> log_cpm()
    lmFit(logCPM, ~group); eBayes(trend=T)  -> moderated_t()
    topTable(..., adjust = "BH")            -> bh_adjust()

Implemented with numpy/scipy rather than a new dependency so that the
dataset can be (re-)ingested from the project's existing requirements.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import special, stats

PRIOR_COUNT = 3  # limma user's guide, limma-trend section: cpm(..., prior.count=3)
NORM_CPM_PREFIX = "norm_cpm:"


# --------------------------------------------------------------------------
# edgeR pieces
# --------------------------------------------------------------------------

def filter_by_expr(
    counts: np.ndarray,
    group: np.ndarray,
    min_count: float = 10,
    min_total_count: float = 15,
    large_n: int = 10,
    min_prop: float = 0.7,
) -> np.ndarray:
    """edgeR::filterByExpr.default with a `group` factor.

    Keeps a gene if it reaches `min_count` reads (expressed as CPM at the
    median library size) in at least as many samples as the smallest
    group, and has at least `min_total_count` reads in total. The rule
    uses group *sizes* only, never which sample is in which group, so it
    cannot favour genes that differ between the groups being compared.
    """
    lib_size = counts.sum(axis=0)
    _, sizes = np.unique(group, return_counts=True)
    min_sample_size = float(sizes[sizes > 0].min())
    if min_sample_size > large_n:
        min_sample_size = large_n + (min_sample_size - large_n) * min_prop
    cpm_cutoff = min_count / np.median(lib_size) * 1e6
    cpm = counts / lib_size * 1e6
    tol = 1e-14
    keep_cpm = (cpm >= cpm_cutoff).sum(axis=1) >= (min_sample_size - tol)
    keep_total = counts.sum(axis=1) >= (min_total_count - tol)
    return keep_cpm & keep_total


def _r_quantile(x: np.ndarray, p: float) -> float:
    # R's default quantile (type 7) is numpy's default 'linear' method.
    return float(np.quantile(x, p))


def _rank(x: np.ndarray) -> np.ndarray:
    # R's rank() default is ties.method = "average".
    return stats.rankdata(x, method="average")


def _tmm_factor(obs, ref, n_obs, n_ref, logratio_trim=0.3, sum_trim=0.05, a_cutoff=-1e10) -> float:
    """edgeR:::.calcFactorTMM with doWeighting = TRUE."""
    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((obs / n_obs) / (ref / n_ref))
        abs_e = (np.log2(obs / n_obs) + np.log2(ref / n_ref)) / 2
        v = (n_obs - obs) / n_obs / obs + (n_ref - ref) / n_ref / ref
    fin = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > a_cutoff)
    log_r, abs_e, v = log_r[fin], abs_e[fin], v[fin]
    # R: max() of an empty vector is -Inf, so an all-non-finite sample
    # returns factor 1 rather than erroring.
    if len(log_r) == 0 or np.max(np.abs(log_r)) < 1e-6:
        return 1.0
    n = len(log_r)
    lo_l = np.floor(n * logratio_trim) + 1
    hi_l = n + 1 - lo_l
    lo_s = np.floor(n * sum_trim) + 1
    hi_s = n + 1 - lo_s
    rl, ra = _rank(log_r), _rank(abs_e)
    keep = (rl >= lo_l) & (rl <= hi_l) & (ra >= lo_s) & (ra <= hi_s)
    f = np.nansum(log_r[keep] / v[keep]) / np.nansum(1 / v[keep])
    if np.isnan(f):
        f = 0.0
    return float(2**f)


def tmm_norm_factors(counts: np.ndarray) -> np.ndarray:
    """edgeR::calcNormFactors(method = "TMM") with default arguments.

    The reference library is the one whose upper-quartile-scaled counts
    sit closest to the mean upper quartile; factors are rescaled to a
    geometric mean of 1.
    """
    counts = counts[(counts > 0).sum(axis=1) > 0]  # edgeR drops all-zero rows first
    lib_size = counts.sum(axis=0)
    scaled = counts / lib_size
    f75 = np.array([_r_quantile(scaled[:, j], 0.75) for j in range(counts.shape[1])])
    if np.median(f75) < 1e-20:
        ref_col = int(np.argmax(np.sqrt(counts).sum(axis=0)))
    else:
        ref_col = int(np.argmin(np.abs(f75 - f75.mean())))
    factors = np.array(
        [
            _tmm_factor(counts[:, j], counts[:, ref_col], lib_size[j], lib_size[ref_col])
            for j in range(counts.shape[1])
        ]
    )
    return factors / np.exp(np.mean(np.log(factors)))


def log_cpm(counts: np.ndarray, lib_size: np.ndarray, prior_count: float = PRIOR_COUNT) -> np.ndarray:
    """edgeR::cpm(log = TRUE): a prior count scaled to each library's size
    is added to every count, and twice that to the library size."""
    prior_scaled = lib_size / lib_size.mean() * prior_count
    adj_lib = lib_size + 2 * prior_scaled
    return np.log2((counts + prior_scaled) / adj_lib * 1e6)


# --------------------------------------------------------------------------
# limma pieces
# --------------------------------------------------------------------------

def _trigamma_inverse(x: float) -> float:
    """limma::trigammaInverse (Newton iteration, Smyth 2004)."""
    if x > 1e7:
        return 1 / np.sqrt(x)
    if x < 1e-6:
        return 1 / x
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if -dif / y < 1e-8:
            break
    return float(y)


def _natural_spline_basis(x: np.ndarray, df: int) -> np.ndarray:
    """A basis spanning the same space as R's splines::ns(x, df, intercept=TRUE).

    ns() places df-2 interior knots at quantiles of x and boundary knots
    at its range. Any basis for natural cubic splines on those knots gives
    identical least-squares fitted values, which is all fitFDist uses, so
    the truncated-power form (Hastie, Tibshirani & Friedman, ESL eq. 5.4-5.5)
    is used here rather than ns()'s B-spline form.
    """
    n_interior = df - 2
    probs = np.linspace(0, 1, n_interior + 2)[1:-1]
    knots = np.concatenate([[x.min()], np.quantile(x, probs), [x.max()]])
    # Rescale for conditioning; an affine change of x maps natural splines
    # to natural splines, so the spanned space is unchanged.
    centre, scale = x.mean(), x.std() or 1.0
    xs, ks = (x - centre) / scale, (knots - centre) / scale
    big_k = len(ks)

    def d(k):
        return (np.clip(xs - ks[k], 0, None) ** 3 - np.clip(xs - ks[-1], 0, None) ** 3) / (ks[-1] - ks[k])

    cols = [np.ones_like(xs), xs] + [d(k) - d(big_k - 2) for k in range(big_k - 2)]
    return np.column_stack(cols)


def fit_f_dist(var: np.ndarray, df1: float, covariate: np.ndarray | None) -> tuple[np.ndarray, float]:
    """limma::fitFDist -- the prior (scale s0^2, df d0) for the gene-wise
    variances. Returns (s0^2 per gene, d0)."""
    x = np.maximum(var, 0)
    m = np.median(x)
    if m == 0:
        m = 1.0
    x = np.maximum(x, 1e-5 * m)
    z = np.log(x)
    e = z - special.digamma(df1 / 2) + np.log(df1 / 2)
    nok = len(e)
    if covariate is None:
        emean = np.full(nok, e.mean())
        evar = np.sum((e - e.mean()) ** 2) / (nok - 1)
    else:
        spline_df = 1 + (nok >= 3) + (nok >= 6) + (nok >= 30)
        spline_df = min(spline_df, len(np.unique(covariate)))
        if spline_df < 2:
            emean = np.full(nok, e.mean())
            evar = np.sum((e - e.mean()) ** 2) / (nok - 1)
        else:
            design = _natural_spline_basis(covariate, spline_df)
            coef, _, rank, _ = np.linalg.lstsq(design, e, rcond=None)
            emean = design @ coef
            evar = np.sum((e - emean) ** 2) / (nok - rank)
    evar = evar - special.polygamma(1, df1 / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s20 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf
        # limma: with no covariate the MLE of the scale is the plain pooled
        # variance (of the floored variances); with a trend, exp(emean).
        # Using exp(emean) in both cases was a deviation (found in review,
        # checked against R on tests/fixtures/contrast/poisson_*).
        s20 = np.full(nok, x.mean()) if covariate is None else np.exp(emean)
    return s20, float(d0)


def moderated_t(log_expr: np.ndarray, is_treatment: np.ndarray, trend: bool = True) -> pd.DataFrame:
    """lmFit(log_expr, ~ group) followed by eBayes(trend = trend), for the
    treatment-vs-reference coefficient. One row per gene, in input order.
    """
    is_treatment = np.asarray(is_treatment, dtype=bool)
    a, b = log_expr[:, ~is_treatment], log_expr[:, is_treatment]
    n_a, n_b = a.shape[1], b.shape[1]
    coef = b.mean(axis=1) - a.mean(axis=1)
    df_resid = n_a + n_b - 2
    rss = ((a - a.mean(axis=1, keepdims=True)) ** 2).sum(axis=1) + ((b - b.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    sigma2 = rss / df_resid
    stdev_unscaled = np.sqrt(1 / n_a + 1 / n_b)
    amean = log_expr.mean(axis=1)

    s20, d0 = fit_f_dist(sigma2, df_resid, amean if trend else None)
    if np.isinf(d0):
        s2_post = s20
    else:
        s2_post = (d0 * s20 + df_resid * sigma2) / (d0 + df_resid)
    df_total = np.minimum(df_resid + d0, df_resid * len(sigma2))
    t = coef / stdev_unscaled / np.sqrt(s2_post)
    p = 2 * stats.t.sf(np.abs(t), df_total)
    return pd.DataFrame(
        {
            "log2_fc": coef,
            "ave_log_cpm": amean,
            "t": t,
            "p_value": p,
            "df_total": np.broadcast_to(df_total, coef.shape),
        }
    ).assign(df_prior=d0)


def _finite_or_none(x: float):
    return x if np.isfinite(x) else None


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg, identical to R's p.adjust(method = "BH")."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    if not ok.all():
        # R's p.adjust drops NA and adjusts over the rest (n = their count);
        # without this one NaN would turn every q into NaN.
        out = np.full(len(p), np.nan)
        out[ok] = bh_adjust(p[ok])
        return out
    n = len(p)
    order = np.argsort(p)[::-1]
    ranked = p[order] * n / np.arange(n, 0, -1)
    q = np.minimum.accumulate(ranked)
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


# --------------------------------------------------------------------------
# The whole pipeline for one contrast
# --------------------------------------------------------------------------

METHOD_LABEL = "limma-trend moderated t on TMM log-CPM (edgeR filterByExpr; BH across tested genes)"


def run_contrast(counts: pd.DataFrame, reference: list[str], treatment: list[str]) -> pd.DataFrame:
    """Differential expression of `treatment` vs `reference` samples.

    `counts` is genes x samples (raw, possibly fractional, counts). Only
    the listed samples enter the fit -- filtering, normalisation and the
    variance prior are all estimated within this experiment alone, so a
    second experiment in the same dataset (another cell line) has no
    influence on this one's numbers.

    Returns every input gene. Genes that fail the expression filter are
    kept with `tested=False` and NaN statistics, so "too lowly expressed
    to test" is reported as such rather than looking like an absent gene.
    log2_fc is treatment minus reference: positive = higher in treatment.
    """
    overlap = set(reference) & set(treatment)
    if overlap:
        raise ValueError(f"Samples in both groups: {sorted(overlap)}")
    if len(reference) < 2 or len(treatment) < 2:
        raise ValueError("Each group needs at least 2 replicates for a variance estimate.")
    samples = list(reference) + list(treatment)
    y = counts[samples].to_numpy(dtype=float)
    is_treatment = np.array([False] * len(reference) + [True] * len(treatment))

    keep = filter_by_expr(y, is_treatment)
    y_kept = y[keep]
    # keep.lib.sizes = FALSE: library sizes are recomputed after filtering,
    # as in the edgeR/limma user's guide workflow.
    lib_size = y_kept.sum(axis=0)
    norm = tmm_norm_factors(y_kept)
    eff_lib = lib_size * norm
    lcpm = log_cpm(y_kept, eff_lib)
    fit = moderated_t(lcpm, is_treatment, trend=True)
    fit["q_value"] = bh_adjust(fit["p_value"].to_numpy())

    out = pd.DataFrame(index=counts.index)
    out["tested"] = keep
    for col in ("log2_fc", "ave_log_cpm", "t", "p_value", "q_value"):
        out[col] = np.nan
        out.loc[keep, col] = fit[col].to_numpy()
    # The normalised counts the test actually ran on (TMM-scaled CPM, no
    # prior), one column per sample. Stored because a dataset's displayed
    # expression unit need not be derived from these counts -- GSE324984's
    # TPM comes from a different, isoform-aware quantification, and for
    # some genes its direction disagrees with this fold-change
    # (METHODS.md 11.3). The reader needs to see the values behind the p.
    norm_cpm = y_kept / eff_lib * 1e6
    for j, sample in enumerate(samples):
        col = f"{NORM_CPM_PREFIX}{sample}"
        out[col] = np.nan
        out.loc[keep, col] = norm_cpm[:, j]
    out.attrs.update(
        {
            "n_genes_input": int(len(counts)),
            "n_genes_tested": int(keep.sum()),
            # JSON (and Postgres JSON) cannot hold Infinity; an infinite
            # prior df is ordinary limma output for very regular data.
            "df_prior": _finite_or_none(float(fit["df_prior"].iloc[0])),
            "df_prior_infinite": bool(np.isinf(fit["df_prior"].iloc[0])),
            "df_residual": int(len(samples) - 2),
            "norm_factors": dict(zip(samples, norm.round(6).tolist())),
        }
    )
    return out


# --------------------------------------------------------------------------
# Reading stored results (used by the API; no statistics computed here)
# --------------------------------------------------------------------------

NOT_TESTED_REASON = (
    "Not tested: too few reads in this experiment's samples to estimate a "
    "fold-change reliably (edgeR filterByExpr)."
)


def _num(x):
    """JSON-safe float: NaN/inf -> None (a missing statistic, not zero)."""
    if x is None:
        return None
    x = float(x)
    return x if np.isfinite(x) else None


def _arm(matrix: pd.DataFrame, feature_id: str, label: str, samples: list[str], table: pd.DataFrame | None = None) -> dict:
    out = {"label": label, "mean_expression": None, "values": [], "mean_norm_cpm": None, "norm_cpm": []}
    if feature_id in matrix.index:
        values = matrix.loc[feature_id, samples]
        out["mean_expression"] = _num(values.mean())
        out["values"] = [{"sample_id": s, "value": _num(values[s])} for s in samples]
    cols = [f"{NORM_CPM_PREFIX}{s}" for s in samples]
    if table is not None and feature_id in table.index and all(c in table.columns for c in cols):
        cpm = table.loc[feature_id, cols].astype(float)
        if cpm.notna().all():
            out["mean_norm_cpm"] = _num(cpm.mean())
            out["norm_cpm"] = [{"sample_id": s, "value": _num(v)} for s, v in zip(samples, cpm)]
    return out


def _display_direction_agrees(stats: dict) -> bool | None:
    """Does the displayed unit's mean change point the same way as the
    count-based log2FC? None when it cannot be judged (untested, or a
    zero/absent mean)."""
    fc = stats.get("log2_fc")
    ref, trt = stats["reference"]["mean_expression"], stats["treatment"]["mean_expression"]
    # A zero arm still has a direction (47.5 -> 0 TPM is "down"); only no
    # change at all, a missing value, or a zero fold-change is unjudgeable.
    if fc is None or ref is None or trt is None or fc == 0 or ref == trt:
        return None
    return (trt > ref) == (fc > 0)


def contrast_gene_stats(matrix: pd.DataFrame, contrast, feature_id: str) -> dict:
    """One gene's stored statistics for one contrast, plus the per-replicate
    expression values behind them."""
    spec = contrast.spec
    out = {
        "contrast_id": spec.contrast_id,
        "in_contrast_table": feature_id in contrast.table.index,
        "reference": _arm(matrix, feature_id, spec.reference, contrast.reference_samples, contrast.table),
        "treatment": _arm(matrix, feature_id, spec.treatment, contrast.treatment_samples, contrast.table),
    }
    if not out["in_contrast_table"]:
        out.update(tested=False, not_tested_reason="Gene absent from this experiment's counts.",
                   display_direction_agrees=None)
        return out
    row = contrast.table.loc[feature_id]
    tested = bool(row["tested"])
    out.update(
        tested=tested,
        not_tested_reason=None if tested else NOT_TESTED_REASON,
        log2_fc=_num(row["log2_fc"]),
        ave_log_cpm=_num(row["ave_log_cpm"]),
        t=_num(row["t"]),
        p_value=_num(row["p_value"]),
        q_value=_num(row["q_value"]),
    )
    out["display_direction_agrees"] = _display_direction_agrees(out)
    return out


def top_contrast_genes(matrix: pd.DataFrame, contrast, n: int = 50, direction: str = "both") -> dict:
    """The contrast's most significant tested genes, by p-value (ties by
    |log2FC|). q-values are the stored genome-wide BH values, so this list
    and the gene panel always show the same number for the same gene."""
    table = contrast.table[contrast.table["tested"].astype(bool)].copy()
    if direction == "up":
        table = table[table["log2_fc"] > 0]
    elif direction == "down":
        table = table[table["log2_fc"] < 0]
    table["_abs"] = table["log2_fc"].abs()
    table = table.sort_values(["p_value", "_abs"], ascending=[True, False]).head(n)
    rows = []
    for fid, r in table.iterrows():
        ref = _arm(matrix, fid, contrast.spec.reference, contrast.reference_samples, contrast.table)
        trt = _arm(matrix, fid, contrast.spec.treatment, contrast.treatment_samples, contrast.table)
        agrees = _display_direction_agrees(
            {"log2_fc": _num(r["log2_fc"]), "reference": ref, "treatment": trt}
        )
        rows.append(
            {
                "feature_id": fid,
                "symbol": r["symbol"],
                "gene_name": r.get("gene_name", r["symbol"]),
                "n_loci_with_name": int(r["n_loci_with_name"]) if "n_loci_with_name" in r.index else 1,
                "log2_fc": _num(r["log2_fc"]),
                "ave_log_cpm": _num(r["ave_log_cpm"]),
                "p_value": _num(r["p_value"]),
                "q_value": _num(r["q_value"]),
                "in_expression_matrix": bool(r.get("in_expression_matrix", fid in matrix.index)),
                "mean_reference": ref["mean_expression"],
                "mean_treatment": trt["mean_expression"],
                "mean_norm_cpm_reference": ref["mean_norm_cpm"],
                "mean_norm_cpm_treatment": trt["mean_norm_cpm"],
                "display_direction_agrees": agrees,
            }
        )
    tested = contrast.table["tested"].astype(bool)
    return {
        "contrast": contrast.summary(),
        "direction": direction,
        # Counted in the same direction as the list it sits beside.
        "n_significant_q05": int(
            (
                (contrast.table.loc[tested, "q_value"] < 0.05)
                & {
                    "up": contrast.table.loc[tested, "log2_fc"] > 0,
                    "down": contrast.table.loc[tested, "log2_fc"] < 0,
                }.get(direction, True)
            ).sum()
        ),
        "genes": rows,
    }
