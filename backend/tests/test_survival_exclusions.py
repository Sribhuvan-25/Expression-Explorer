"""
Survival must count each exclusion reason separately (METHODS.md 6.1, 8.6).

A TARGET quartile split read "235 excluded (no follow-up time recorded,
or score fell in the excluded middle band ...)": 3 patients without
follow-up and 232 cut by the split, behind one number.
"""
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from app import registry
from app.api import main as api_main
from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

DATASET_ID = "survival_accounting_test"
N = 40
NO_FOLLOW_UP = {"P00", "P01"}


def _build():
    rng = np.random.default_rng(5)
    samples_ids = [f"P{i:02d}" for i in range(N)]
    genes = [f"G{i}" for i in range(398)] + ["SPREAD", "TIED"]
    matrix = pd.DataFrame(rng.uniform(1, 100, size=(len(genes), N)), index=genes, columns=samples_ids)
    # SPREAD: rank i + 1 in sample i (aucThreshold is 100), so every
    # sample gets a different score and no scores tie.
    for i, s in enumerate(samples_ids):
        others = np.sort(matrix.loc[genes[:398], s].to_numpy())[::-1]
        matrix.loc["SPREAD", s] = others[i] + 1e-6
    # TIED: bottom-ranked (score exactly 0) in half the cohort.
    matrix.loc["TIED"] = [0.0 if i < N // 2 else 200.0 + i for i in range(N)]
    samples = pd.DataFrame(
        {
            "sample_id": samples_ids,
            "dataset_id": DATASET_ID,
            "group_columns": [{} for _ in range(N)],
            "vital_status": ["Dead" if i % 2 else "Alive" for i in range(N)],
            "days_to_death": [None if s in NO_FOLLOW_UP or i % 2 == 0 else 300.0 + i for i, s in enumerate(samples_ids)],
            "days_to_last_follow_up": [
                None if s in NO_FOLLOW_UP or i % 2 else 900.0 + i for i, s in enumerate(samples_ids)
            ],
        }
    )
    features = pd.DataFrame({"feature_id": genes, "symbol": genes, "aliases": [[] for _ in genes]})
    source = DatasetSource(
        dataset_id=DATASET_ID,
        display_name="Survival accounting test",
        accession="none",
        repository="test",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.TPM,
        n_samples=N,
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


@pytest.fixture(scope="module")
def client():
    if DATASET_ID not in {d.dataset_id for d in registry.list_descriptors()}:
        registry.register(
            registry.DatasetDescriptor(
                dataset_id=DATASET_ID,
                display_name="Survival accounting test",
                loader=_build,
                group_columns=("vital_status",),
                supports_survival=True,
            )
        )
    api_main._get_dataset.cache_clear()
    return TestClient(api_main.app)


def _run(client, **body):
    r = client.post(f"/datasets/{DATASET_ID}/survival", json={"genes": ["SPREAD"], **body})
    assert r.status_code == 200, r.text
    return r.json()



def test_median_split_excludes_only_patients_without_follow_up(client):
    body = _run(client)
    assert body["exclusion_breakdown"] == [{"reason": "with no follow-up time recorded", "n": 2}]
    assert body["n_excluded"] == 2


def test_quartile_reports_follow_up_and_middle_band_separately(client):
    body = _run(client, cutoff_method="quartile")
    arms = sum(c["n"] for c in body["curves"].values())
    reasons = {e["reason"]: e["n"] for e in body["exclusion_breakdown"]}
    assert reasons == {
        "with no follow-up time recorded": 2,
        "with a score between the low (25%) and high (25%) cutoffs": 38 - arms,
    }
    assert body["n_excluded"] == sum(reasons.values()) == N - arms
    assert body["exclusion_reason"].startswith("2 with no follow-up time recorded; ")


def test_samples_dropped_at_a_tied_cutoff_are_their_own_reason(client):
    trimmed = _run(client, genes=["TIED"], cutoff_method="quartile", tie_policy="trim")
    tie = [e for e in trimmed["exclusion_breakdown"] if "tie policy" in e["reason"]]
    assert len(tie) == 1 and tie[0]["n"] > 0 and "'trim'" in tie[0]["reason"]
    assert trimmed["n_excluded"] == sum(e["n"] for e in trimmed["exclusion_breakdown"])
    assert trimmed["n_excluded"] == N - sum(c["n"] for c in trimmed["curves"].values())

    inclusive = _run(client, genes=["TIED"], cutoff_method="quartile", tie_policy="inclusive")
    assert not [e for e in inclusive["exclusion_breakdown"] if "tie policy" in e["reason"]]


def test_response_names_the_endpoint_and_carries_censoring_and_at_risk(client):
    """The figure must say it is overall survival and draw censoring ticks
    and an at-risk table from the analysis (none of which the API returned)."""
    body = _run(client)
    assert body["endpoint"]["name"] == "Overall survival"
    assert "at least t" in body["at_risk_convention"]
    for c in body["curves"].values():
        # Every patient in an arm is an event or a censoring; in this cohort
        # the alive patients are censored at distinct times 900 + i.
        assert c["n_events"] + c["n_censored"] == c["n"]
        assert sum(m["n_censored"] for m in c["censor_marks"]) == c["n_censored"]
        assert all(m["days"] >= 900 for m in c["censor_marks"])
        assert c["at_risk"][0] == {"days": 0.0, "n_at_risk": c["n"]}
