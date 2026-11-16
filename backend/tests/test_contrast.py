"""
The count-based contrast pipeline (app/analysis/contrast.py) must
reproduce R's edgeR + limma-trend exactly -- it is a re-implementation
whose whole claim to correctness is "same numbers as the reference
software".

Fixture: 2,720 genes of the real GSE324984 counts (random expressed
genes, 200 all-zero genes, and every PP2A-panel gene), with R's output on
that same subset vendored beside it. Produced by
`limma_ref.R` (edgeR 3.40.2 / limma 3.54.2, R 4.2):

    d <- DGEList(y, group); keep <- filterByExpr(d, group = group)
    d <- d[keep, , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
    lcpm <- cpm(d, log = TRUE, prior.count = 3)
    fit <- eBayes(lmFit(lcpm, model.matrix(~group)), trend = TRUE)
    topTable(fit, coef = 2, number = Inf, sort.by = "none")

On the full 78,691-gene file the same comparison agrees to <= 2e-12
absolute on every statistic (METHODS.md 11.2); the subset keeps the test
offline and fast.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from app.analysis.contrast import bh_adjust, run_contrast

FIXTURES = Path(__file__).parent / "fixtures" / "contrast"

DESIGNS = {
    "c42b": (["C4-2B_Ctrl_1", "C4-2B_Ctrl_2"], ["C4-2B_MYCL_1", "C4-2B_MYCL_2"]),
    "h660": (["shCtrl_H660_1", "shCtrl_H660_2"], ["shMYCL_H660_1", "shMYCL_H660_2"]),
}


@pytest.fixture(scope="module")
def counts():
    frame = pd.read_csv(FIXTURES / "gse324984_counts_subset.tsv.gz", sep="\t")
    frame.index = frame["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    return frame.drop(columns=["gene_id", "gene_name"])


@pytest.mark.parametrize("tag", sorted(DESIGNS))
def test_matches_r_limma_trend(counts, tag):
    reference, treatment = DESIGNS[tag]
    ours = run_contrast(counts, reference, treatment)
    r = pd.read_csv(FIXTURES / f"{tag}_limma.csv.gz").set_index("feature_id")

    # Same genes pass filterByExpr.
    assert set(ours.index[ours["tested"]]) == set(r.index)
    tested = ours.loc[r.index]
    for ours_col, r_col in [
        ("log2_fc", "logFC"),
        ("ave_log_cpm", "AveExpr"),
        ("t", "t"),
        ("p_value", "P"),
        ("q_value", "adjP"),
    ]:
        np.testing.assert_allclose(tested[ours_col], r[r_col], rtol=1e-8, atol=1e-12, err_msg=ours_col)


@pytest.mark.parametrize("tag", sorted(DESIGNS))
def test_tmm_factors_match_r(counts, tag):
    reference, treatment = DESIGNS[tag]
    ours = run_contrast(counts, reference, treatment).attrs["norm_factors"]
    r = pd.read_csv(FIXTURES / f"{tag}_nf.csv.gz").set_index("sample")["nf"]
    np.testing.assert_allclose([ours[s] for s in r.index], r.to_numpy(), atol=1e-6)


def test_bh_matches_r_p_adjust():
    p = np.array([0.0005, 0.004, 0.004, 0.03, 0.2, 0.21, 0.6, 0.99])
    # R: p.adjust(p, "BH")
    expected = [0.004, 0.010666666666666666, 0.010666666666666666, 0.06, 0.28, 0.28, 0.6857142857142856, 0.99]
    np.testing.assert_allclose(bh_adjust(p), expected, rtol=1e-12)
    # Order of input must not matter.
    shuffled = np.array([3, 0, 7, 5, 1, 6, 2, 4])
    np.testing.assert_allclose(bh_adjust(p[shuffled]), np.array(expected)[shuffled], rtol=1e-12)


def test_other_experiments_samples_do_not_influence_a_contrast(counts):
    """The H660 columns sit in the same counts frame as the C4-2B ones.
    A C4-2B contrast must be computed from the C4-2B samples alone --
    filtering, TMM and the variance prior included."""
    reference, treatment = DESIGNS["c42b"]
    alone = run_contrast(counts[reference + treatment], reference, treatment)
    alongside = run_contrast(counts, reference, treatment)
    pd.testing.assert_frame_equal(alone, alongside, check_exact=True)


def test_untested_genes_are_reported_not_dropped(counts):
    reference, treatment = DESIGNS["c42b"]
    result = run_contrast(counts, reference, treatment)
    assert len(result) == len(counts)
    untested = result[~result["tested"]]
    assert len(untested) > 0
    assert untested[["log2_fc", "p_value", "q_value"]].isna().all().all()
    assert result.loc[result["tested"], "p_value"].notna().all()


def test_direction_is_treatment_minus_reference(counts):
    reference, treatment = DESIGNS["c42b"]
    forward = run_contrast(counts, reference, treatment)
    reverse = run_contrast(counts, treatment, reference)
    tested = forward["tested"]
    np.testing.assert_allclose(forward.loc[tested, "log2_fc"], -reverse.loc[tested, "log2_fc"], atol=1e-12)
    # MYCL (ENSG00000116990) is the overexpressed gene: must be strongly positive.
    assert forward.loc["ENSG00000116990", "log2_fc"] > 5


def test_rejects_overlapping_or_unreplicated_groups(counts):
    with pytest.raises(ValueError, match="both groups"):
        run_contrast(counts, ["C4-2B_Ctrl_1", "C4-2B_Ctrl_2"], ["C4-2B_Ctrl_2", "C4-2B_MYCL_1"])
    with pytest.raises(ValueError, match="at least 2"):
        run_contrast(counts, ["C4-2B_Ctrl_1"], ["C4-2B_MYCL_1", "C4-2B_MYCL_2"])


@pytest.mark.parametrize("tag", sorted(DESIGNS))
def test_stored_normalised_cpm_matches_r_library_scaling(counts, tag):
    """norm_cpm:<sample> = count / (filtered lib size x TMM factor) x 1e6,
    with lib sizes and factors taken from R's DGEList -- an independent
    source, not this module's own intermediates."""
    reference, treatment = DESIGNS[tag]
    ours = run_contrast(counts, reference, treatment)
    r = pd.read_csv(FIXTURES / f"{tag}_nf.csv.gz").set_index("sample")
    tested = ours.index[ours["tested"]]
    for sample in reference + treatment:
        expected = counts.loc[tested, sample] / (r.loc[sample, "lib"] * r.loc[sample, "nf"]) * 1e6
        np.testing.assert_allclose(ours.loc[tested, f"norm_cpm:{sample}"], expected, rtol=1e-6)
    assert ours.loc[~ours["tested"], [f"norm_cpm:{s}" for s in reference]].isna().all().all()


def test_display_direction_agreement_cases():
    from app.analysis.contrast import _display_direction_agrees as agrees

    def stats(fc, ref, trt):
        return {"log2_fc": fc, "reference": {"mean_expression": ref}, "treatment": {"mean_expression": trt}}

    assert agrees(stats(1.0, 10, 20)) is True
    assert agrees(stats(-1.0, 10, 20)) is False  # HSPA12A-style disagreement
    assert agrees(stats(None, 10, 20)) is None  # untested
    assert agrees(stats(1.0, 0, 20)) is True  # zero arm still has a direction
    assert agrees(stats(1.0, 20, 0)) is False
    assert agrees(stats(1.0, 5, 5)) is None  # no change to judge
    assert agrees(stats(1.0, None, 20)) is None  # gene not in TPM matrix


def test_bh_ignores_nan_like_r():
    # R: p.adjust(c(0.01, NA, 0.5, 0.02), "BH") -> 0.03, NA, 0.5, 0.03
    q = bh_adjust(np.array([0.01, np.nan, 0.5, 0.02]))
    assert np.isnan(q[1])
    np.testing.assert_allclose(q[[0, 2, 3]], [0.03, 0.5, 0.03], rtol=1e-12)


def test_bh_ignores_nan_like_r_second_case():
    # R: p.adjust(c(0.01, NA, 0.03, 0.04), "BH") -> 0.03 NA 0.04 0.04
    out = bh_adjust(np.array([0.01, np.nan, 0.03, 0.04]))
    np.testing.assert_allclose(out[[0, 2, 3]], [0.03, 0.04, 0.04])
    assert np.isnan(out[1])


@pytest.mark.parametrize("trend", [True, False])
def test_infinite_prior_df_matches_r(trend):
    """Very regular (Poisson) replicates: limma's prior df is infinite and
    the posterior variance is the prior. R reference from limma_edge.R."""
    from app.analysis.contrast import log_cpm, moderated_t, tmm_norm_factors, filter_by_expr

    frame = pd.read_csv(FIXTURES / "poisson_infinite_prior.tsv.gz", sep="\t")
    frame.index = frame["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    y = frame[["R1", "R2", "T1", "T2"]].to_numpy(float)
    grp = np.array([False, False, True, True])
    keep = filter_by_expr(y, grp)
    yk = y[keep]
    lib = yk.sum(axis=0) * tmm_norm_factors(yk)
    fit = moderated_t(log_cpm(yk, lib), grp, trend=trend)
    assert np.isinf(fit["df_prior"].iloc[0])
    r = pd.read_csv(FIXTURES / f"poisson_trend_{str(trend).upper()}.csv.gz").set_index("feature_id")
    np.testing.assert_allclose(fit["t"].to_numpy(), r.loc[frame.index[keep], "t"].to_numpy(), rtol=1e-8)
    np.testing.assert_allclose(fit["p_value"].to_numpy(), r.loc[frame.index[keep], "P"].to_numpy(), rtol=1e-8)


def test_infinite_prior_df_is_stored_json_safe():
    import json

    frame = pd.read_csv(FIXTURES / "poisson_infinite_prior.tsv.gz", sep="\t").set_index("gene_id")
    result = run_contrast(frame[["R1", "R2", "T1", "T2"]], ["R1", "R2"], ["T1", "T2"])
    assert result.attrs["df_prior"] is None and result.attrs["df_prior_infinite"] is True
    json.dumps(result.attrs, allow_nan=False)


def test_tmm_tolerates_a_sample_with_no_finite_ratios():
    from app.analysis.contrast import _tmm_factor

    obs = np.zeros(5)
    ref = np.array([1.0, 2, 3, 4, 5])
    assert _tmm_factor(obs, ref, 1.0, ref.sum()) == 1.0
