"""HTTP-level tests for the design-contrast endpoints (METHODS.md 11)
and the small-n Mann-Whitney floor (1.4a), on a synthetic GSE324984."""
import json
import math

import pytest
from fastapi.testclient import TestClient

from app.api import main as api_main
from app.ingest import gse324984
from tests.test_gse324984 import MYCL_ID, _synthetic, write_dataset  # noqa: F401  (fixture)


@pytest.fixture()
def client(write_dataset, monkeypatch):  # noqa: F811
    ids, names, *_ = _synthetic()
    named = list(names)
    named[50], named[60], named[70] = "CIP2A", "SILENT", "C1orf112"

    def silence(counts):
        counts.loc[60, :] = 0.0  # never expressed -> fails filterByExpr

    write_dataset(tpm_names=named, count_names=named, mutate_counts=silence)
    ds = gse324984.load()
    monkeypatch.setattr(api_main, "_get_dataset", lambda dataset_id: ds)
    return TestClient(api_main.app), ds, ids


def test_lists_both_within_cell_line_contrasts(client):
    c, _, _ = client
    body = c.get("/datasets/gse324984/contrasts").json()
    by_id = {x["contrast_id"]: x for x in body["contrasts"]}
    assert set(by_id) == {"c42b_mycl_oe", "h660_mycl_kd"}
    assert by_id["c42b_mycl_oe"]["stratum_value"] == "C4-2B"
    assert by_id["h660_mycl_kd"]["n_reference"] == by_id["h660_mycl_kd"]["n_treatment"] == 2


def test_panel_rows_missing_and_alias_disclosure(client):
    c, ds, ids = client
    r = c.post(
        "/datasets/gse324984/contrast-panel",
        json={"genes": ["MYCL", "kiaa1524", "NOPE1", "mycl", "SILENT"]},
    )
    assert r.status_code == 200
    # Strict JSON: an untested gene's NaN statistics must arrive as null.
    json.loads(r.text, parse_constant=lambda c: pytest.fail(f"non-JSON constant {c}"))
    body = r.json()
    rows = {row["symbol"]: row for row in body["rows"]}
    assert set(rows) == {"MYCL", "CIP2A", "SILENT"}
    assert rows["CIP2A"]["alias_note"] == "'kiaa1524' is an older name for CIP2A"
    assert rows["MYCL"]["alias_note"] is None
    reasons = {m["gene"]: m["reason"] for m in body["missing"]}
    assert set(reasons) == {"NOPE1", "mycl"}
    assert "Duplicate" in reasons["mycl"]

    oe = next(x for x in rows["MYCL"]["contrasts"] if x["contrast_id"] == "c42b_mycl_oe")
    kd = next(x for x in rows["MYCL"]["contrasts"] if x["contrast_id"] == "h660_mycl_kd")
    assert oe["tested"] and oe["log2_fc"] > 0 > kd["log2_fc"]
    # Replicate values are the dataset's own TPM, for that experiment's samples only.
    assert [v["sample_id"] for v in oe["treatment"]["values"]] == ["GSM9589990", "GSM9589991"]
    assert oe["treatment"]["values"][0]["value"] == pytest.approx(ds.matrix.loc[MYCL_ID, "GSM9589990"])
    assert oe["q_value"] == pytest.approx(ds.contrasts["c42b_mycl_oe"].table.loc[MYCL_ID, "q_value"])

    silent = rows["SILENT"]["contrasts"][0]
    assert silent["tested"] is False and silent["p_value"] is None
    assert "filterByExpr" in silent["not_tested_reason"]


def test_top_genes_agree_with_panel_and_respect_direction(client):
    c, ds, _ = client
    top = c.get("/datasets/gse324984/contrasts/c42b_mycl_oe/top", params={"n": 20}).json()
    ps = [g["p_value"] for g in top["genes"]]
    assert ps == sorted(ps)
    assert top["genes"][0]["symbol"] == "MYCL"
    table = ds.contrasts["c42b_mycl_oe"].table
    for g in top["genes"]:
        assert g["q_value"] == pytest.approx(table.loc[g["feature_id"], "q_value"])
    up = c.get("/datasets/gse324984/contrasts/c42b_mycl_oe/top", params={"direction": "up"}).json()
    assert all(g["log2_fc"] > 0 for g in up["genes"])
    down = c.get("/datasets/gse324984/contrasts/c42b_mycl_oe/top", params={"direction": "down"}).json()
    assert all(g["log2_fc"] < 0 for g in down["genes"])


def test_unknown_contrast_and_contrastless_dataset(client):
    c, ds, _ = client
    r = c.get("/datasets/gse324984/contrasts/nope/top")
    assert r.status_code == 404 and "c42b_mycl_oe" in r.json()["detail"]
    ds.contrasts, saved = {}, ds.contrasts
    try:
        r = c.post("/datasets/gse324984/contrast-panel", json={"genes": ["MYCL"]})
        assert r.status_code == 422
    finally:
        ds.contrasts = saved


def test_legacy_symbol_404_names_the_current_symbol(client):
    c, _, _ = client
    r = c.get("/datasets/gse324984/compare", params={"gene": "KIAA1524", "group_column": "arm"})
    assert r.status_code == 404
    assert "CIP2A" in r.json()["detail"]


def test_compare_reports_the_mann_whitney_floor(client):
    c, _, _ = client
    body = c.get("/datasets/gse324984/compare", params={"gene": "MYCL", "group_column": "arm"}).json()
    floors = [t["min_p_achievable"] for t in body["pairwise_tests"]]
    assert len(floors) == 6 and floors == pytest.approx([1 / 3] * 6)
    # And the actual p never undercuts it on these untied values.
    assert all(t["p_value"] >= t["min_p_achievable"] - 1e-12 for t in body["pairwise_tests"])


def test_alias_hint_works_on_symbol_indexed_datasets():
    """GSE66586/GDS4299 index the matrix by symbol, so feature_id == symbol
    on every row. An earlier version filtered those rows out of the alias
    lookup and the hint silently never appeared there (caught by hand
    against the running API)."""
    import pandas as pd
    from fastapi import HTTPException

    from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

    genes = ["CIP2A", "MYCN"]
    ds = Dataset(
        source=DatasetSource(
            dataset_id="sym", display_name="sym", accession="x", repository="x",
            assay_type=AssayType.MICROARRAY, expression_unit=ExpressionUnit.LOG2_INTENSITY, n_samples=1,
        ),
        matrix=pd.DataFrame({"S1": [1.0, 2.0]}, index=genes),
        samples=pd.DataFrame({"sample_id": ["S1"], "dataset_id": ["sym"], "group_columns": [{}]}),
        features=pd.DataFrame({"feature_id": genes, "symbol": genes, "aliases": [["KIAA1524"], []]}),
    )
    with pytest.raises(HTTPException) as exc:
        api_main._resolve_gene(ds, "kiaa1524")
    assert "older name for CIP2A" in exc.value.detail
    assert api_main._resolve_panel_gene(ds, "KIAA1524") == ("CIP2A", "'KIAA1524' is an older name for CIP2A")


def test_significant_count_follows_direction(client):
    c, ds, _ = client
    t = ds.contrasts["c42b_mycl_oe"].table
    t = t[t["tested"].astype(bool)]
    for direction, mask in (("up", t.log2_fc > 0), ("down", t.log2_fc < 0), ("both", t.log2_fc == t.log2_fc)):
        body = c.get("/datasets/gse324984/contrasts/c42b_mycl_oe/top", params={"direction": direction}).json()
        assert body["n_significant_q05"] == int(((t.q_value < 0.05) & mask).sum())


def test_mixed_case_symbols_resolve_on_an_ensembl_indexed_dataset(client):
    """Upper-casing the query could never match 'C1orf112'; it 404'd on
    every Ensembl-indexed dataset (TARGET, GSE324984) while present."""
    c, _, ids = client
    for query in ("C1orf112", "c1orf112", " C1ORF112 "):
        r = c.get("/datasets/gse324984/compare", params={"gene": query, "group_column": "arm"})
        assert r.status_code == 200, (query, r.json())


def test_case_only_symbol_collision_is_refused_not_guessed():
    import pandas as pd
    from fastapi import HTTPException

    from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

    ids = ["ENSG1", "ENSG2"]
    ds = Dataset(
        source=DatasetSource(dataset_id="x", display_name="x", accession="x", repository="x",
                             assay_type=AssayType.RNA_SEQ, expression_unit=ExpressionUnit.TPM, n_samples=1),
        matrix=pd.DataFrame({"S1": [1.0, 2.0]}, index=ids),
        samples=pd.DataFrame({"sample_id": ["S1"], "dataset_id": ["x"], "group_columns": [{}]}),
        features=pd.DataFrame({"feature_id": ids, "symbol": ["Abc1", "ABc1"], "aliases": [[], []]}),
    )
    assert api_main._resolve_gene(ds, "Abc1") == "ENSG1"  # exact still wins
    with pytest.raises(HTTPException) as exc:
        api_main._resolve_gene(ds, "abc1")
    assert exc.value.status_code == 422


def test_compare_multi_skip_reason_names_the_current_symbol(client):
    c, _, _ = client
    body = c.get("/datasets/compare-multi", params={"gene": "KIAA1524", "group_column": "arm",
                                                     "dataset_ids": "gse324984"}).json()
    entry = next(d for d in body["datasets"] if d["dataset_id"] == "gse324984")
    assert entry["skipped"] and "CIP2A" in entry["skip_reason"]


def test_cross_cell_line_pairs_are_marked(client):
    c, ds, _ = client
    body = c.get("/datasets/gse324984/compare", params={"gene": "MYCL", "group_column": "arm"}).json()
    crossing = [t for t in body["pairwise_tests"] if t["crosses_stratum"]]
    within = [t for t in body["pairwise_tests"] if not t["crosses_stratum"]]
    assert len(crossing) == 4 and len(within) == 2
    for t in within:
        assert t["group_a"].split()[0] == t["group_b"].split()[0]  # same cell line
    assert {crossing[0]["crosses_stratum"]["group_a"], crossing[0]["crosses_stratum"]["group_b"]} == {"C4-2B", "NCI-H660"}
    assert "cell_line" in body["stratum_caveat"]
    # Grouping within one stratum's own column structure: no contrasts -> no caveat.
    ds.contrasts, saved = {}, ds.contrasts
    try:
        body = c.get("/datasets/gse324984/compare", params={"gene": "MYCL", "group_column": "arm"}).json()
        assert body["stratum_caveat"] is None and all(t["crosses_stratum"] is None for t in body["pairwise_tests"])
    finally:
        ds.contrasts = saved


@pytest.fixture()
def review_client(write_dataset, monkeypatch):  # noqa: F811
    """The independent review's cases: a CAST-like name whose expressed
    copy exists only past the TPM truncation, and a mixed-case name shared
    by two TPM rows (the C4orf36 case)."""
    from tests.test_gse324984 import N_TPM

    tail = N_TPM + 20
    _, names, *_ = _synthetic()
    named = list(names)
    named[30] = named[tail] = "CASTLIKE"
    named[40] = named[41] = "C4orf36"

    def shape(counts):
        counts.loc[30, :] = 1.0
        counts.loc[tail, :] = 5000.0
        counts.loc[40, :] = 4000.0  # primary C4orf36
        counts.loc[41, :] = 3.0

    ids, *_ = write_dataset(tpm_names=named, count_names=named, mutate_counts=shape)
    ds = gse324984.load()
    monkeypatch.setattr(api_main, "_get_dataset", lambda dataset_id: ds)
    return TestClient(api_main.app), ds, ids, tail


def test_count_only_gene_is_a_panel_row_with_its_statistics(review_client):
    c, ds, ids, tail = review_client
    body = c.post("/datasets/gse324984/contrast-panel", json={"genes": ["CASTLIKE", ids[tail]]}).json()
    assert len(body["rows"]) == 1 and body["missing"][0]["reason"].startswith("Duplicate")
    row = body["rows"][0]
    assert row["feature_id"] == ids[tail] and row["symbol"] == "CASTLIKE"
    assert row["in_expression_matrix"] is False
    assert any("absent from the dataset's expression matrix" in n for n in row["identity_notes"])
    oe = row["contrasts"][0]
    assert oe["tested"] and oe["p_value"] is not None and oe["reference"]["values"] == []
    assert oe["reference"]["mean_norm_cpm"] is not None  # the values the test used are still shown


def test_compare_on_a_count_only_gene_explains_rather_than_misleads(review_client):
    c, _, ids, tail = review_client
    r = c.get("/datasets/gse324984/compare", params={"gene": "CASTLIKE", "group_column": "arm"})
    assert r.status_code == 404
    assert ids[tail] in r.json()["detail"] and "no value in its expression matrix" in r.json()["detail"]


def test_mixed_case_shared_name_resolves_to_its_primary_locus(review_client):
    """Was: 'C4orf36' upper-cased, missed, then the alias fallback matched
    the secondary row's shared-name alias and served the silent locus."""
    c, _, ids, _ = review_client
    body = c.post("/datasets/gse324984/contrast-panel", json={"genes": ["C4orf36"]}).json()
    row = body["rows"][0]
    assert row["feature_id"] == ids[40] and row["alias_note"] is None
    assert any("shared by 2 loci" in n and "most-expressed" in n for n in row["identity_notes"])
    secondary = c.post("/datasets/gse324984/contrast-panel", json={"genes": [ids[41]]}).json()["rows"][0]
    assert secondary["symbol"] == ids[41] and any("secondary locus" in n for n in secondary["identity_notes"])


def test_top_list_and_panel_show_the_same_identity(review_client):
    c, ds, ids, _ = review_client
    top = c.get("/datasets/gse324984/contrasts/c42b_mycl_oe/top", params={"n": 500}).json()
    by_id = {g["feature_id"]: g for g in top["genes"]}
    for fid in (ids[40], ids[41]):
        if fid in by_id:
            panel = c.post("/datasets/gse324984/contrast-panel", json={"genes": [fid]}).json()["rows"][0]
            assert by_id[fid]["symbol"] == panel["symbol"]
    assert all(g["symbol"] != "C4orf36" or g["feature_id"] == ids[40] for g in top["genes"])


def test_contrast_responses_carry_sample_accounting(review_client):
    c, *_ = review_client
    for body in (
        c.get("/datasets/gse324984/contrasts").json(),
        c.post("/datasets/gse324984/contrast-panel", json={"genes": ["MYCL"]}).json(),
    ):
        for contrast in body["contrasts"]:
            assert contrast["n_dataset_total"] == 8 and contrast["n_excluded"] == 4
            assert "By design" in contrast["exclusion_reason"]
    top = c.get("/datasets/gse324984/contrasts/h660_mycl_kd/top").json()
    assert top["contrast"]["n_excluded"] == 4 and "NCI-H660" in top["contrast"]["exclusion_reason"]


def test_compare_multi_and_pca_carry_the_resolver_explanation(review_client):
    c, _, ids, tail = review_client
    body = c.get("/datasets/compare-multi", params={"gene": "CASTLIKE", "group_column": "arm",
                                                     "dataset_ids": "gse324984"}).json()
    entry = next(d for d in body["datasets"] if d["dataset_id"] == "gse324984")
    assert entry["skipped"] and ids[tail] in entry["skip_reason"]
    pca = c.post("/datasets/gse324984/pca", json={"genes": ["MYCL", "C4orf36", "CASTLIKE"]}).json()
    assert pca["genes_unrecognized"] == ["CASTLIKE"]
    assert ids[tail] in pca["genes_unrecognized_reasons"]["CASTLIKE"]
