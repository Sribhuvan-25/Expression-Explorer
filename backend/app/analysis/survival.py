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


def binarize_by_cutoff(
    df: pd.DataFrame,
    method: str = "median",
    custom_high_pct: float = 50.0,
    custom_low_pct: float = 50.0,
    score_col: str = "signature_score",
    tie_policy: str = "trim",
) -> pd.Series:
    """Group samples by signature score, with a null (excluded) entry for
    anyone in neither group -- generalizes `binarize_by_median` to match
    GEPIA3's three cutoff modes (see reference/gepia3-feature-review.md):

    - "median": every sample is LOW or HIGH, split at the median. Same
      result as `binarize_by_median`, just returned through this shared
      interface so callers only need one function.
    - "quartile": only the top 25% (HIGH) and bottom 25% (LOW) by score
      are kept -- the middle 50% get `None` and are excluded from
      everything downstream (KM curves, log-rank, Cox). This is a
      materially different comparison from median-split, not just a
      different threshold: it trades n (half the cohort, by design) for
      a cleaner separation between "clearly high" and "clearly low"
      scores, which is why it needs its own explicit mode rather than
      being folded into a "percentile" parameter on the median function.
    - "custom": same top/bottom idea as quartile, but with independently
      chosen percentiles for each side (e.g. top 30% vs bottom 40%) --
      GEPIA3 exposes these as separate Cutoff-High/Cutoff-Low fields
      rather than one symmetric percentile, so a real use case is an
      asymmetric split, not just "quartile with a different number".

    `tie_policy` decides what happens when a percentile threshold lands
    exactly on a block of equal scores (see the body for why this is
    common here, not an edge case). Only meaningful for quartile/custom:
    - "trim" (default) -- keep the arm at its requested size, dropping
      surplus tied samples by a deterministic rank order.
    - "exclude" -- drop the whole tied block from both arms.
    - "inclusive" -- the original behaviour: every tied sample joins the
      arm, so arms can be much larger than the requested percentile.
      Kept so an expert can reproduce prior numbers deliberately.

    Callers must apply the same exclusion-transparency pattern used
    elsewhere in this app (METHODS.md 6.1) whenever method != "median":
    the returned Series has real `None` entries, and how many of the
    input samples got excluded is exactly `.isna().sum()`.
    """
    if method == "median":
        return binarize_by_median(df, score_col=score_col)

    if method == "quartile":
        high_pct, low_pct = 25.0, 25.0
    elif method == "custom":
        high_pct, low_pct = custom_high_pct, custom_low_pct
    else:
        raise ValueError(f"Unknown cutoff method: {method!r}")

    if tie_policy not in ("trim", "exclude", "inclusive"):
        raise ValueError(
            f"Unknown tie_policy: {tie_policy!r} (expected 'trim', 'exclude', or 'inclusive')"
        )

    for name, pct in (("custom_high_pct", high_pct), ("custom_low_pct", low_pct)):
        if not (0 < pct < 100):
            raise ValueError(f"{name} must be between 0 and 100 (exclusive), got {pct}")
    if high_pct + low_pct > 100:
        raise ValueError(
            f"High ({high_pct}%) + low ({low_pct}%) cutoffs exceed 100% of the cohort -- "
            "the two groups would overlap."
        )

    scores = df[score_col]
    # Quantile convention: cutoff_high_pct=25 means "top 25% by score", so
    # the threshold is the 75th percentile (scores ABOVE it are the top
    # 25%) -- not the 25th percentile, which would be the bottom quarter's
    # own threshold. Symmetric for the low side.
    high_threshold = scores.quantile(1 - high_pct / 100)
    low_threshold = scores.quantile(low_pct / 100)

    group = pd.Series(pd.NA, index=df.index, dtype="object")
    group[scores >= high_threshold] = "HIGH"
    group[scores <= low_threshold] = "LOW"

    if tie_policy == "inclusive":
        return group

    # A threshold can land exactly ON a tied block of scores, and then the
    # inclusive comparison above sweeps the whole block into one arm. This
    # is not hypothetical: `auc_signature_score` returns exactly 0.0 for
    # every sample where the gene never reaches the top-25% rank cut, so
    # for a single gene expressed in under ~75% of the cohort the 25th
    # percentile IS 0.0. MEF2C on TARGET gave LOW n=217 (46% of the
    # cohort) against HIGH n=118, while the UI still called it a quartile
    # split -- the log-rank p and Cox HR were computed on arms that were
    # not what the caption claimed (caught in QA).
    #
    # There is no single correct answer here, which is why the caller
    # chooses and the result reports what happened:
    #   "trim"    -- keep the requested arm size, breaking the tie by rank
    #                order. Arms match the stated definition, but which
    #                equals land in the arm is arbitrary (deterministic).
    #   "exclude" -- drop the entire tied block from both arms. Never puts
    #                equal scores in different arms, at the cost of n.
    for label, threshold, ascending in (
        ("HIGH", high_threshold, False),
        ("LOW", low_threshold, True),
    ):
        in_arm = group == label
        tied = in_arm & (scores == threshold)
        n_tied = int(tied.sum())
        if n_tied <= 1:
            continue
        # Requested size for this arm, from the percentile the caller asked
        # for -- not the arm's current (tie-inflated) size.
        pct = high_pct if label == "HIGH" else low_pct
        target = int(round(len(scores) * pct / 100))
        untied = int((in_arm & ~tied).sum())
        if untied >= target:
            # The untied members already fill the arm; every tied sample is
            # surplus under either policy.
            group[tied] = pd.NA
            continue
        if tie_policy == "exclude":
            group[tied] = pd.NA
            continue
        keep = target - untied
        if keep >= n_tied:
            continue
        # Deterministic: order the tied block by sample id so the same
        # query always yields the same arms.
        surplus = sorted(scores.index[tied])[keep:]
        group[surplus] = pd.NA

    return group


def kaplan_meier_curves(df: pd.DataFrame, group: pd.Series) -> dict[str, dict]:
    """One KM curve per group; returns step-function points ready to plot,
    plus a log-rank p-value between the two groups when there are exactly
    two (matching the paper's two-group LOW/HIGH comparisons).

    `group` may contain nulls (quartile/custom cutoffs deliberately exclude
    a middle band, see `binarize_by_cutoff`) -- those samples don't belong
    to any curve and must be dropped before iterating labels, not plotted
    as a spurious "None"/"<NA>" group.
    """
    labeled = group.dropna()
    df = df.loc[labeled.index]
    group = labeled

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
