import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from lifelines.exceptions import ConvergenceError, StatError

from app.analysis.survival import (
    binarize_by_cutoff,
    binarize_by_median,
    build_survival_frame,
    cox_model,
    kaplan_meier_curves,
)


def _samples(n=30, with_group_columns=True):
    rng = np.random.default_rng(0)
    vital = np.where(rng.random(n) < 0.5, "Dead", "Alive")
    df = pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(n)],
            "vital_status": vital,
            "days_to_death": np.where(vital == "Dead", rng.uniform(50, 2000, n), np.nan),
            "days_to_last_follow_up": np.where(vital == "Alive", rng.uniform(50, 2000, n), np.nan),
        }
    )
    if with_group_columns:
        etp = rng.choice(["ETP", "near-ETP", "non-ETP"], n)
        df["group_columns"] = [{"etp_status": e} for e in etp]
    return df


def test_build_survival_frame_flattens_nested_group_columns_onto_the_frame():
    """etp_status is a group column TARGET-ALL-P2 exposes, but it lives
    inside each sample's nested group_columns dict rather than as a
    top-level field. cox_model selects covariates by column name on the
    survival frame, so requesting 'etp_status' as a covariate raised a raw
    KeyError ("['etp_status'] not in index") -- an unhandled 500 -- purely
    because the column didn't exist on the frame at all, regardless of
    whether the fit itself would succeed."""
    samples = _samples(40)
    scores = pd.Series(np.random.default_rng(1).random(40), index=samples["sample_id"])

    frame = build_survival_frame(samples, scores)

    assert "etp_status" in frame.columns


def test_cox_model_raw_categorical_covariate_raises_rather_than_500s():
    """cox_model does not one-hot encode covariates (no covariate picker
    exists in the UI today, so this path is API-only) -- a raw string
    column like etp_status fails inside lifelines' fit with a plain
    ValueError ('could not convert string to float'). The API route
    catches ValueError/ConvergenceError/StatError around cox_model and
    degrades to cox: None rather than propagating a 500; this test locks
    in that the failure mode stays within that catchable set instead of
    silently becoming some other exception type on a future lifelines
    upgrade."""
    samples = _samples(40)
    scores = pd.Series(np.random.default_rng(1).random(40), index=samples["sample_id"])
    frame = build_survival_frame(samples, scores)

    with pytest.raises((ValueError, ConvergenceError, StatError)):
        cox_model(frame, covariates=["etp_status"])


def test_cox_model_unknown_covariate_raises_value_error_not_keyerror():
    samples = _samples(30, with_group_columns=False)
    scores = pd.Series(np.random.default_rng(2).random(30), index=samples["sample_id"])
    frame = build_survival_frame(samples, scores)

    with pytest.raises(ValueError, match="Unknown covariate"):
        cox_model(frame, covariates=["not_a_real_column"])


def test_cox_model_includes_hazard_ratio_confidence_interval():
    """A hazard ratio without its CI gives no sense of how precisely it's
    estimated -- exactly the thing a reader checks it for."""
    samples = _samples(60, with_group_columns=False)
    rng = np.random.default_rng(3)
    scores = pd.Series(rng.random(60), index=samples["sample_id"])
    frame = build_survival_frame(samples, scores)

    result = cox_model(frame)
    coef = result["coefficients"]["signature_score"]
    assert "exp(coef) lower 95%" in coef
    assert "exp(coef) upper 95%" in coef
    assert coef["exp(coef) lower 95%"] <= coef["exp(coef)"] <= coef["exp(coef) upper 95%"]


def test_cox_model_too_few_cases_raises_value_error():
    samples = _samples(5, with_group_columns=False)
    scores = pd.Series(np.random.default_rng(4).random(5), index=samples["sample_id"])
    frame = build_survival_frame(samples, scores)

    with pytest.raises(ValueError, match="Too few complete cases"):
        cox_model(frame)


def _scored_frame(n=100, seed=5):
    samples = _samples(n, with_group_columns=False)
    scores = pd.Series(np.random.default_rng(seed).random(n), index=samples["sample_id"])
    return build_survival_frame(samples, scores)


def test_binarize_by_cutoff_median_matches_binarize_by_median():
    frame = _scored_frame()
    assert binarize_by_cutoff(frame, method="median").equals(binarize_by_median(frame))


def test_binarize_by_cutoff_quartile_keeps_only_top_and_bottom_quarter():
    frame = _scored_frame(n=100)
    group = binarize_by_cutoff(frame, method="quartile")

    assert (group == "HIGH").sum() == 25
    assert (group == "LOW").sum() == 25
    assert group.isna().sum() == 50

    # The HIGH group really is the top scores, not an arbitrary 25%.
    high_scores = frame.loc[group == "HIGH", "signature_score"]
    low_scores = frame.loc[group == "LOW", "signature_score"]
    assert high_scores.min() >= low_scores.max()


def test_binarize_by_cutoff_custom_asymmetric_percentiles():
    frame = _scored_frame(n=100)
    group = binarize_by_cutoff(frame, method="custom", custom_high_pct=30, custom_low_pct=10)

    assert (group == "HIGH").sum() == 30
    assert (group == "LOW").sum() == 10
    assert group.isna().sum() == 60


def test_binarize_by_cutoff_rejects_overlapping_percentiles():
    frame = _scored_frame(n=50)
    with pytest.raises(ValueError, match="exceed 100%"):
        binarize_by_cutoff(frame, method="custom", custom_high_pct=60, custom_low_pct=60)


def test_binarize_by_cutoff_rejects_out_of_range_percentile():
    frame = _scored_frame(n=50)
    with pytest.raises(ValueError, match="between 0 and 100"):
        binarize_by_cutoff(frame, method="custom", custom_high_pct=0, custom_low_pct=50)


def test_binarize_by_cutoff_unknown_method_raises():
    frame = _scored_frame(n=50)
    with pytest.raises(ValueError, match="Unknown cutoff method"):
        binarize_by_cutoff(frame, method="bogus")


def test_kaplan_meier_curves_drops_excluded_middle_band():
    """A quartile/custom cutoff's excluded middle band must not surface as
    a spurious third curve (e.g. a 'None'/'<NA>' group) -- only HIGH/LOW
    should ever appear, and only the kept samples should count toward n."""
    frame = _scored_frame(n=100)
    group = binarize_by_cutoff(frame, method="quartile")

    result = kaplan_meier_curves(frame, group)

    assert set(result["curves"].keys()) == {"HIGH", "LOW"}
    assert result["curves"]["HIGH"]["n"] == 25
    assert result["curves"]["LOW"]["n"] == 25
    assert "logrank_p_value" in result  # still exactly 2 groups -> log-rank runs
