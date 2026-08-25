"""
Kaplan-Meier + Cox proportional hazards on a signature score, matching the
paper's method (binarize by median signature score, log-rank test, Cox
model via survfit-equivalent). Answers PDF comment 3 (PP2A survival) and
rebuilds the mechanics behind Figure 10 C/D/N/O/P for any signature.
"""
from __future__ import annotations

import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test


def build_survival_frame(
    samples: pd.DataFrame,
    signature_scores: pd.Series,
) -> pd.DataFrame:
    """Combine TARGET-style clinical fields with a signature score into the
    (duration, event, score) shape survival analysis needs.

    duration = days_to_death if the patient died, else days_to_last_follow_up
    event = 1 if vital_status == 'Dead', else 0 (right-censored)
    """
    df = samples.copy()
    df["signature_score"] = signature_scores.reindex(df["sample_id"]).to_numpy()
    df["event"] = (df["vital_status"] == "Dead").astype(int)
    df["duration"] = df["days_to_death"].where(df["event"] == 1, df["days_to_last_follow_up"])
    # A covariate a caller wants to control for (e.g. etp_status, MRD
    # status) may live inside each sample's nested group_columns dict
    # rather than as a top-level field -- flatten it here so cox_model can
    # select it by name either way, instead of raising KeyError.
    if "group_columns" in df.columns:
        meta = pd.json_normalize(df["group_columns"])
        meta.index = df.index
        for col in meta.columns:
            if col not in df.columns:
                df[col] = meta[col]
    return df.dropna(subset=["duration", "signature_score"])


def binarize_by_median(df: pd.DataFrame, score_col: str = "signature_score") -> pd.Series:
    """LOW/HIGH split at the median, matching the paper's binarization
    approach for signature-based survival groups (e.g. ZMIZ1-5^LOW/HIGH)."""
    median = df[score_col].median()
    return (df[score_col] > median).map({True: "HIGH", False: "LOW"})


def kaplan_meier_curves(df: pd.DataFrame, group: pd.Series) -> dict[str, dict]:
    """One KM curve per group; returns step-function points ready to plot,
    plus a log-rank p-value between the two groups when there are exactly
    two (matching the paper's two-group LOW/HIGH comparisons)."""
    curves = {}
    for label in group.unique():
        mask = group == label
        kmf = KaplanMeierFitter()
        kmf.fit(df.loc[mask, "duration"], event_observed=df.loc[mask, "event"], label=str(label))
        sf = kmf.survival_function_.reset_index()
        sf.columns = ["days", "survival_probability"]
        curves[str(label)] = {
            "n": int(mask.sum()),
            "points": sf.to_dict("records"),
        }

    result: dict = {"curves": curves}
    labels = list(group.unique())
    if len(labels) == 2:
        m1, m2 = group == labels[0], group == labels[1]
        test = logrank_test(
            df.loc[m1, "duration"], df.loc[m2, "duration"],
            event_observed_A=df.loc[m1, "event"], event_observed_B=df.loc[m2, "event"],
        )
        result["logrank_p_value"] = test.p_value
    return result


def cox_model(df: pd.DataFrame, covariates: list[str] | None = None) -> dict:
    """Cox proportional hazards with the signature score as the primary
    predictor, plus optional covariates (paper controls for day-29 MRD,
    CNS status, diagnostic age, WBC where available)."""
    covariates = covariates or []
    unknown = [c for c in covariates if c not in df.columns]
    if unknown:
        raise ValueError(f"Unknown covariate(s): {', '.join(unknown)}")
    cols = ["duration", "event", "signature_score"] + covariates
    cph_df = df[cols].dropna()
    if len(cph_df) < 10:
        raise ValueError("Too few complete cases for a Cox model (need >= 10).")
    cph = CoxPHFitter()
    cph.fit(cph_df, duration_col="duration", event_col="event")
    summary = cph.summary
    # Include the HR confidence interval lifelines already computes --
    # without it a hazard ratio has no indication of how precisely it's
    # estimated, which matters most for exactly the case a reader would
    # check it: is this HR actually distinguishable from 1.
    cols = ["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]
    return {
        "n": len(cph_df),
        "coefficients": summary[cols].to_dict("index"),
        "log_likelihood_p_value": cph.log_likelihood_ratio_test().p_value,
    }
