"""HTTP-level tests for the API surface.

The rest of the suite tests analysis functions directly, which left
`app/api/main.py` at 26% coverage -- and every bug found in that file
this session (a case-sensitive `aux_feature` while `gene` was not, a 500
on mygene.info's list-shaped `ensembl`, a route shadowed by a dynamic
path) was caught by a human or an agent clicking the UI, never by a
test. Those are all HTTP-shape bugs: the analysis function was correct
and the endpoint around it was not, so testing one layer down could not
have found them.

These tests run against a synthetic dataset registered into the registry,
so they never touch the network or load a real cohort.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from app import registry
from app.api import main as api_main
from app.models.contract import (
    AssayType,
    AuxLayer,
    AuxMatrix,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)

DATASET_ID = "apitest"
N = 24


def _build_dataset() -> Dataset:
    rng = np.random.default_rng(0)
    sample_ids = [f"S{i:02d}" for i in range(N)]
    genes = ["MYCN", "DTX1", "NOTCH1", "MEF2C", "CD34"]
    matrix = pd.DataFrame(
        rng.normal(6.0, 2.0, size=(len(genes), N)).clip(0), index=genes, columns=sample_ids
    )
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "dataset_id": DATASET_ID,
            "etp_status": ["ETP" if i % 2 else "non-ETP" for i in range(N)],
            "vital_status": ["Dead" if i % 3 == 0 else "Alive" for i in range(N)],
            "days_to_death": [100.0 + i if i % 3 == 0 else None for i in range(N)],
            "days_to_last_follow_up": [200.0 + i if i % 3 else None for i in range(N)],
            "group_columns": [{} for _ in sample_ids],
        }
    )
    features = pd.DataFrame(
        {"feature_id": genes, "symbol": genes, "aliases": [[] for _ in genes]}
    )
    source = DatasetSource(
        dataset_id=DATASET_ID,
        display_name="API Test Cohort",
        accession="none",
        repository="test",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.TPM,
        n_samples=N,
    )
    # Aux layer covering only the first 16 samples, and a feature measured
    # on only some of those -- the partial-coverage shape that the
    # exclusion accounting exists for.
    covered = sample_ids[:16]
    crispr = pd.DataFrame(
        {s: [float(i), float(i) if i < 8 else np.nan] for i, s in enumerate(covered)},
        index=["DENSE", "SPARSE"],
    )
    aux = {
        AuxLayer.CRISPR_GENE_EFFECT: AuxMatrix(
            layer=AuxLayer.CRISPR_GENE_EFFECT,
            matrix=crispr,
            value_label="Gene effect",
            value_description="test",
            source_note="test",
        ),
        AuxLayer.MUTATION_STATUS: AuxMatrix(
            layer=AuxLayer.MUTATION_STATUS,
            matrix=pd.DataFrame(
                {s: [i % 2 == 0, False] for i, s in enumerate(covered)},
                index=["NOTCH1", "NEVER"],
            ),
            value_label="Non-silent mutation present",
            value_description="test",
            source_note="test",
        ),
    }
    return Dataset(source=source, matrix=matrix, samples=samples, features=features, aux=aux)


@pytest.fixture(scope="module")
def client():
    if DATASET_ID not in {d.dataset_id for d in registry.list_descriptors()}:
        registry.register(
            registry.DatasetDescriptor(
                dataset_id=DATASET_ID,
                display_name="API Test Cohort",
                loader=_build_dataset,
                # mrd_status is declared but deliberately never populated in
            # _build_dataset() -- this mirrors the real §7.10 incident
            # shape (a supplementary clinical grouping column whose
            # upstream fetch can fail while the column itself stays
            # declared on the dataset) and is what
            # test_group_values_reports_unavailable_for_a_declared_but_empty_column
            # exercises.
            group_columns=("etp_status", "vital_status", "mrd_status"),
                supports_survival=True,
            )
        )
    api_main._get_dataset.cache_clear()
    return TestClient(api_main.app)


def test_health(client):
    assert client.get("/health").status_code == 200


def test_dataset_listing_includes_provenance(client):
    body = client.get("/datasets").json()
    entry = next(d for d in body["datasets"] if d["dataset_id"] == DATASET_ID)
    # Provenance fields exist because a QA pass found a correctly-computed
    # staleness note that no UI ever read (METHODS.md 7.5).
    assert entry["n_samples"] == N
    assert entry["assay_type"] == "rna_seq"
    assert "notes" in entry


def test_unknown_dataset_is_404_not_500(client):
    assert client.get("/datasets/nope/compare?gene=MYCN&group_column=etp_status").status_code == 404


def test_compare_multi_is_not_shadowed_by_the_dynamic_route(client):
    """`/datasets/compare-multi` must not be captured as
    dataset_id="compare-multi". FastAPI matches in registration order, so
    this breaks silently if the routes are ever reordered."""
    r = client.get(f"/datasets/compare-multi?gene=MYCN&group_column=etp_status&dataset_ids={DATASET_ID}")
    assert r.status_code == 200
    assert r.json()["datasets"][0]["dataset_id"] == DATASET_ID


def test_compare_multi_never_pools_across_datasets(client):
    """METHODS.md 7.2 -- results come back per dataset, each with its own
    unit, and no merged statistic anywhere in the payload."""
    body = client.get(
        f"/datasets/compare-multi?gene=MYCN&group_column=etp_status&dataset_ids={DATASET_ID}"
    ).json()
    assert "expression_unit" in body["datasets"][0]
    assert not any("pool" in k.lower() for k in body)


@pytest.mark.parametrize("gene", ["MYCN", "mycn", "  MyCn  "])
def test_gene_lookup_is_case_and_whitespace_insensitive(client, gene):
    """One endpoint uppercased its own parameter while the shared resolver
    did not, so `?gene=dtx1` 404'd on an API the UI happens to uppercase
    for you (METHODS.md 8.6)."""
    r = client.get(f"/datasets/{DATASET_ID}/compare", params={"gene": gene, "group_column": "etp_status"})
    assert r.status_code == 200


def test_unknown_gene_is_404_naming_the_gene(client):
    r = client.get(
        f"/datasets/{DATASET_ID}/compare", params={"gene": "ZZZNOPE", "group_column": "etp_status"}
    )
    assert r.status_code == 404
    assert "ZZZNOPE" in r.json()["detail"]


def test_group_values_404s_on_unknown_column(client):
    """Previously returned an empty list, which reads as "this column
    exists and has no values" rather than "no such column"."""
    ok = client.get(f"/datasets/{DATASET_ID}/group-values?group_column=etp_status")
    assert ok.status_code == 200 and ok.json()["values"]
    assert client.get(f"/datasets/{DATASET_ID}/group-values?group_column=nope").status_code == 404


def test_group_values_reports_unavailable_for_a_declared_but_empty_column(client):
    """§7.10: a column can be genuinely declared on the dataset (it's in
    group_columns, so NOT a 404) while having zero non-null values on
    this particular load -- the exact shape produced by a degraded
    third-party clinical-supplement fetch (mrd_status/etp_status on
    target_all_p2). Distinguishing this from the 404 case, and from
    "still loading", is what lets the UI explain an empty picker rather
    than showing one with no options and no reason why."""
    r = client.get(f"/datasets/{DATASET_ID}/group-values?group_column=mrd_status")
    assert r.status_code == 200
    body = r.json()
    assert body["values"] == []
    assert body["unavailable"] is True


def test_group_values_unavailable_is_false_for_a_populated_column(client):
    r = client.get(f"/datasets/{DATASET_ID}/group-values?group_column=etp_status")
    assert r.status_code == 200
    assert r.json()["unavailable"] is False


def test_correlation_rejects_identical_genes(client):
    """Duplicate column names after the transpose made scipy receive a
    DataFrame rather than a Series -- an unhandled 500."""
    r = client.get(f"/datasets/{DATASET_ID}/correlation?gene_a=MYCN&gene_b=MYCN")
    assert r.status_code == 422


def test_differential_q_is_never_below_p(client):
    body = client.get(
        f"/datasets/{DATASET_ID}/differential",
        params={"group_column": "etp_status", "group_a": "ETP", "group_b": "non-ETP", "limit": 5},
    ).json()
    rows = body["genes"]
    assert rows, "expected differential rows"
    assert all(r["q_value"] + 1e-12 >= r["p_value"] for r in rows)


def test_aux_layers_reports_usable_features_separately(client):
    body = client.get(f"/datasets/{DATASET_ID}/aux-layers").json()
    crispr = next(l for l in body["layers"] if l["layer"] == "crispr_gene_effect")
    # SPARSE is measured on 8 of 16 covered samples, DENSE on all 16 --
    # both clear the minimum here, so usable == total for this fixture.
    assert crispr["n_features"] == 2
    assert crispr["n_samples_covered"] == 16
    assert crispr["n_dataset_total"] == N


@pytest.mark.parametrize("feature", ["DENSE", "dense", "  DeNsE  "])
def test_aux_feature_lookup_is_case_insensitive(client, feature):
    """`?gene=myb` resolved while `?aux_feature=myb` 404'd with an error
    that read as a data-coverage claim (METHODS.md 8.7)."""
    r = client.get(
        f"/datasets/{DATASET_ID}/expression-vs-aux",
        params={"gene": "MYCN", "layer": "crispr_gene_effect", "aux_feature": feature},
    )
    assert r.status_code == 200, r.text
    assert r.json()["aux_feature_label"] == "DENSE"


def test_expression_vs_aux_splits_exclusion_causes_over_http(client):
    """The two causes must stay distinguishable through the API, not just
    in the analysis function (METHODS.md 8.6)."""
    body = client.get(
        f"/datasets/{DATASET_ID}/expression-vs-aux",
        params={"gene": "MYCN", "layer": "crispr_gene_effect", "aux_feature": "SPARSE"},
    ).json()
    assert body["n_missing_layer"] + body["n_missing_feature"] == body["n_excluded"]
    assert body["n"] + body["n_excluded"] == body["n_dataset_total"]


def test_aux_features_picker_only_offers_usable_features(client):
    """The picker and the correlation share MIN_AUX_SAMPLES; if they ever
    disagreed the UI would list features that immediately fail."""
    body = client.get(
        f"/datasets/{DATASET_ID}/aux-features", params={"layer": "crispr_gene_effect"}
    ).json()
    for feature in body["features"]:
        r = client.get(
            f"/datasets/{DATASET_ID}/expression-vs-aux",
            params={
                "gene": "MYCN",
                "layer": "crispr_gene_effect",
                "aux_feature": feature["feature_id"],
            },
        )
        assert r.status_code == 200, f"picker offered {feature['feature_id']} but it failed: {r.text}"


def test_compare_by_mutation_rejects_gene_with_no_mutations_in_cohort(client):
    """2,701 of TARGET's 6,250 mutated genes are mutated only in samples
    with no expression data; these used to return a degenerate 200."""
    r = client.get(
        f"/datasets/{DATASET_ID}/compare-by-mutation", params={"gene": "MYCN", "mutated_gene": "NEVER"}
    )
    assert r.status_code == 422


def test_survival_tie_policies_are_reported_and_change_the_arms(client):
    """METHODS.md 3.1a -- the policy applied must be visible, since it
    changes what the arms mean."""
    seen = {}
    for policy in ("trim", "exclude", "inclusive"):
        r = client.post(
            f"/datasets/{DATASET_ID}/survival",
            json={"genes": ["MYCN"], "cutoff_method": "quartile", "tie_policy": policy},
        )
        assert r.status_code == 200, r.text
        body = r.json()
        if body.get("tie_note"):
            assert body["tie_note"]["tie_policy"] == policy
        seen[policy] = {k: v["n"] for k, v in body["curves"].items()}
    assert seen  # all three policies are accepted by the API


def test_survival_rejects_unknown_tie_policy(client):
    r = client.post(
        f"/datasets/{DATASET_ID}/survival",
        json={"genes": ["MYCN"], "cutoff_method": "quartile", "tie_policy": "nonsense"},
    )
    assert r.status_code == 422


def test_survival_quartile_message_reports_the_real_percentiles(client):
    """A quartile run echoed the unused request defaults and called itself
    a 50%/50% cutoff."""
    body = client.post(
        f"/datasets/{DATASET_ID}/survival",
        json={"genes": ["MYCN"], "cutoff_method": "quartile"},
    ).json()
    if body["n_excluded"]:
        assert "50%" not in body["exclusion_reason"]
        assert "25%" in body["exclusion_reason"]


def test_survival_unsupported_dataset_is_a_clear_400(client):
    registry_ids = {d.dataset_id for d in registry.list_descriptors()}
    non_survival = [
        d for d in registry.list_descriptors() if not d.supports_survival and d.dataset_id in registry_ids
    ]
    if not non_survival:
        pytest.skip("no non-survival dataset registered")
    r = client.post(
        f"/datasets/{non_survival[0].dataset_id}/survival", json={"genes": ["MYCN"]}
    )
    assert r.status_code in (400, 404)


def test_pca_requires_an_explicit_gene_set(client):
    """No auto-default: picking the gene set for the user would decide the
    analysis for them."""
    assert client.post(f"/datasets/{DATASET_ID}/pca", json={"genes": []}).status_code == 422
    ok = client.post(f"/datasets/{DATASET_ID}/pca", json={"genes": ["MYCN", "DTX1", "MEF2C"]})
    assert ok.status_code == 200
    body = ok.json()
    # Variance-explained is a list of fractions, one per component. They
    # are fractions of total variance, so they can never sum above 1 --
    # a UI that renders them as percentages would print >100% if they did.
    assert body["genes_used"] == ["MYCN", "DTX1", "MEF2C"]
    assert len(body["variance_explained"]) == body["n_components"]
    assert all(0.0 <= v <= 1.0 for v in body["variance_explained"])
    assert sum(body["variance_explained"]) <= 1.0 + 1e-9


def test_health_is_independent_of_dataset_readiness(client):
    """/health is the platform's deploy healthcheck. Gating it on data
    being downloaded is what turned a slow GDC pull into a failed
    deployment -- it must answer on process liveness alone."""
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_readiness_reports_per_dataset_state(client):
    r = client.get("/readiness")
    assert r.status_code == 200
    body = r.json()
    for key in ("finished", "n_ready", "n_failed", "n_total", "datasets"):
        assert key in body


def test_datasets_listing_never_blocks_on_a_warming_dataset(client, monkeypatch):
    """/datasets is the frontend's first call. On a cold container a
    dataset can be minutes from ready, so this endpoint must report it as
    warming rather than calling the loader and hanging the whole UI."""
    from app.api import main as api_main

    monkeypatch.setattr(
        api_main.warmup_tracker,
        "snapshot",
        lambda: {
            "finished": False,
            "n_ready": 0,
            "n_failed": 0,
            "n_total": 1,
            "datasets": [{"dataset_id": DATASET_ID, "state": "loading"}],
        },
    )

    def _explode(dataset_id):  # pragma: no cover - must never be reached
        raise AssertionError(f"/datasets called the loader for {dataset_id} while it was warming")

    monkeypatch.setattr(api_main, "_get_dataset", _explode)

    body = client.get("/datasets").json()
    entry = next(d for d in body["datasets"] if d["dataset_id"] == DATASET_ID)
    assert entry["warming"] is True
    # Provenance is absent precisely because it was not waited for.
    assert "n_samples" not in entry


def test_warm_datasets_are_listed_with_provenance_and_not_flagged(client):
    body = client.get("/datasets").json()
    entry = next(d for d in body["datasets"] if d["dataset_id"] == DATASET_ID)
    assert entry["warming"] is False
    assert entry["n_samples"] == N
