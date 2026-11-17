"""Fields the Kaplan-Meier figure needs beyond the estimate itself: the
endpoint the curves measure, censoring marks, and numbers at risk.

A clinician reviewer could not publish the figure without them: it did not
say it was overall survival, drew no censoring ticks (the API returned no
censoring times), and had no at-risk table. Every expected value below is
worked by hand from a six-patient example, not read back from lifelines.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from app.analysis.survival import ENDPOINT, kaplan_meier_curves

# Six patients in one arm (event 1 = died, 0 = censored):
#   day 1 died, day 3 died, day 3 censored, day 5 died, day 8 censored,
#   day 10 died.
# KM: S(1) = 5/6; S(3) = 5/6 * 4/5 = 2/3 (one death of 5 at risk; the
# patient censored at 3 is still at risk at 3); S(5) = 2/3 * 2/3 = 4/9;
# S(8) = 4/9 (censoring only); S(10) = 0.
# At risk (follow-up >= t): t=0 -> 6, 1 -> 6, 3 -> 5, 5 -> 3, 8 -> 2, 10 -> 1.
DAYS = [1, 3, 3, 5, 8, 10]
EVENTS = [1, 1, 0, 1, 0, 1]


def _curves():
    df = pd.DataFrame(
        {"duration": DAYS + [2, 4, 6, 7], "event": EVENTS + [1, 0, 1, 0]},
        index=[f"S{i}" for i in range(10)],
    )
    group = pd.Series(["A"] * 6 + ["B"] * 4, index=df.index)
    return kaplan_meier_curves(df, group)["curves"]


def test_censor_marks_sit_at_each_censoring_time_on_the_curve():
    marks = _curves()["A"]["censor_marks"]
    assert [m["days"] for m in marks] == [3.0, 8.0]
    # At day 3 a death and a censoring coincide: the mark sits on S(3)
    # after the drop, as R's plot.survfit draws it.
    assert marks[0]["survival_probability"] == pytest.approx(2 / 3)
    assert marks[1]["survival_probability"] == pytest.approx(4 / 9)
    assert [m["n_censored"] for m in marks] == [1, 1]


def test_event_and_censored_counts_add_up_to_the_arm():
    a = _curves()["A"]
    assert (a["n_events"], a["n_censored"], a["n"]) == (4, 2, 6)


def test_at_risk_counts_follow_up_at_least_t():
    at_risk = {r["days"]: r["n_at_risk"] for r in _curves()["A"]["at_risk"]}
    assert at_risk == {0.0: 6, 1.0: 6, 3.0: 5, 5.0: 3, 8.0: 2, 10.0: 1}


def test_figure_fields_leave_the_estimate_unchanged():
    points = _curves()["A"]["points"]
    assert [p["days"] for p in points] == [0.0, 1.0, 3.0, 5.0, 8.0, 10.0]
    assert [p["survival_probability"] for p in points] == pytest.approx([1, 5 / 6, 2 / 3, 4 / 9, 4 / 9, 0])


def test_endpoint_is_overall_survival():
    assert ENDPOINT["name"] == "Overall survival"
    assert "death from any cause" in ENDPOINT["event"]
    assert ENDPOINT["time_unit"] == "days"
