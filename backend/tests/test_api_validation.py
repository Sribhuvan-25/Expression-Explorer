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


def test_mutation_layer_usable_count_counts_mutants_not_non_nulls(client):  # noqa: F811
    """A true/false layer is never null, so every gene counted as usable
    (TARGET: 5,762 of 5,762 while 2,916 have no mutant)."""
    layers = {l["layer"]: l for l in client.get(f"/datasets/{DATASET_ID}/aux-layers").json()["layers"]}
    mut = layers["mutation_status"]
    # Fixture: NOTCH1 mutated in 8 of 16 covered samples, NEVER in none.
    assert mut["n_features"] == 2
    assert mut["n_usable_features"] == 1


def test_gds4299_records_the_legacy_symbols_it_remaps():
    """KIAA1524 / PME1 are renamed at ingest; without them as aliases the
    resolver's 'older name for CIP2A' hint could never fire on GDS4299."""
    from app.ingest.gds4299 import legacy_aliases

    assert legacy_aliases(["CIP2A", "PPME1", "MYCN"]) == [["KIAA1524"], ["PME1"], []]


def test_legacy_name_hint_from_gene_annotation(client, monkeypatch):  # noqa: F811
    """DepMap and TARGET record no aliases, so a legacy name got a bare
    'not found' while the gene was present under its current symbol."""
    record = {"symbol": "MYCN", "aliases": ["NMYC", "MODED"]}
    monkeypatch.setattr(api_main, "lookup_gene", lambda s: record)
    r = client.get(f"/datasets/{DATASET_ID}/compare", params={"gene": "nmyc", "group_column": "etp_status"})
    assert r.status_code == 404
    assert "alias of MYCN" in r.json()["detail"] and "search for MYCN" in r.json()["detail"]

    # A loose free-text hit that does not list the name as an alias: no hint.
    monkeypatch.setattr(api_main, "lookup_gene", lambda s: {"symbol": "MYCN", "aliases": []})
    r = client.get(f"/datasets/{DATASET_ID}/compare", params={"gene": "ZZZ9", "group_column": "etp_status"})
    assert r.json()["detail"] == "Gene 'ZZZ9' not found in dataset."

    def boom(s):
        raise OSError("offline")

    monkeypatch.setattr(api_main, "lookup_gene", boom)
    r = client.get(f"/datasets/{DATASET_ID}/compare", params={"gene": "NMYC", "group_column": "etp_status"})
    assert r.status_code == 404 and r.json()["detail"] == "Gene 'NMYC' not found in dataset."


def test_gene_annotation_rejects_blank_or_oversized_symbols(client):  # noqa: F811
    """Both answered 500 (empty cache filename; name too long)."""
    assert client.get("/genes/%20%20%20").status_code == 422
    assert client.get("/genes/" + "A" * 5000).status_code == 422
    # Network is blocked in tests: the lookup fails, which is a 503, not a 500.
    assert client.get("/genes/MYCN").status_code in (200, 503)
