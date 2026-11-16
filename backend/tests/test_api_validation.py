"""
Input handling and accounting gaps found by an independent audit
(2026-11-16): each returned a 500 or a misleading 200.
"""
import pandas as pd
import pytest

from app.analysis.ranking import rank_by_gene
from app.analysis.survival import kaplan_meier_curves
from app.api import main as api_main
from app.models.contract import AuxLayer
from tests.test_api import DATASET_ID, client  # noqa: F401  (fixture)


def test_unknown_dataset_is_404_on_survival_and_group_values(client):  # noqa: F811
    assert client.post("/datasets/nope/survival", json={"genes": ["MYCN"]}).status_code == 404
    assert client.get("/datasets/nope/group-values", params={"group_column": "x"}).status_code == 404


def test_compare_rejects_a_column_that_is_not_a_grouping(client):  # noqa: F811
    r = client.get(f"/datasets/{DATASET_ID}/compare", params={"gene": "MYCN", "group_column": "days_to_death"})
    assert r.status_code == 422 and "not a grouping column" in r.json()["detail"]


def test_compare_multi_names_an_unknown_dataset(client):  # noqa: F811
    r = client.get(
        "/datasets/compare-multi",
        params={"gene": "MYCN", "group_column": "etp_status", "dataset_ids": f"{DATASET_ID},typo_ds"},
    )
    assert r.status_code == 404 and "typo_ds" in r.json()["detail"]


def test_differential_rejects_a_group_against_itself_and_accounts_for_samples(client):  # noqa: F811
    same = client.get(
        f"/datasets/{DATASET_ID}/differential",
        params={"group_column": "etp_status", "group_a": "ETP", "group_b": "ETP"},
    )
    assert same.status_code == 422
    body = client.get(
        f"/datasets/{DATASET_ID}/differential",
        params={"group_column": "vital_status", "group_a": "Dead", "group_b": "Alive"},
    ).json()
    assert body["n_dataset_total"] == body["n_a"] + body["n_b"] + body["n_excluded"]
    assert body["n_excluded"] == sum(e["n"] for e in body["exclusion_breakdown"])


def test_mutated_gene_with_a_lower_case_symbol_resolves(client):  # noqa: F811
    ds = api_main._get_dataset(DATASET_ID)
    aux = ds.aux[AuxLayer.MUTATION_STATUS]
    original = aux.matrix
    try:
        extra = pd.DataFrame([[i % 3 == 0 for i in range(original.shape[1])]], index=["C11orf53"], columns=original.columns)
        aux.matrix = pd.concat([original, extra])
        r = client.get(
            f"/datasets/{DATASET_ID}/compare-by-mutation", params={"gene": "MYCN", "mutated_gene": "C11orf53"}
        )
        assert r.status_code == 200, r.text
        r = client.get(
            f"/datasets/{DATASET_ID}/compare-by-mutation", params={"gene": "MYCN", "mutated_gene": "c11orf53"}
        )
        assert r.status_code == 200, r.text
    finally:
        aux.matrix = original


def test_rank_rows_omit_survival_times_kept_inside_group_columns():
    samples = pd.DataFrame(
        {
            "sample_id": ["A", "B"],
            "group_columns": [
                {"subtype": "x", "days_to_death": 10.0, "days_to_last_follow_up": None},
                {"subtype": "y", "days_to_death": None, "days_to_last_follow_up": 99.0},
            ],
        }
    )
    matrix = pd.DataFrame([[1.0, 2.0]], index=["G"], columns=["A", "B"])
    rows = rank_by_gene(matrix, samples, "G")
    assert "subtype" in rows.columns
    assert "days_to_death" not in rows.columns and "days_to_last_follow_up" not in rows.columns


def test_no_log_rank_p_for_an_arm_of_one():
    df = pd.DataFrame({"duration": [10.0, 20.0, 30.0], "event": [1, 0, 1]}, index=["a", "b", "c"])
    group = pd.Series(["HIGH", "LOW", None], index=df.index, dtype=object)
    assert kaplan_meier_curves(df, group)["logrank_p_value"] is None
