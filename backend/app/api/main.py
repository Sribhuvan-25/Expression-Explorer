"""
API surface. Thin by design — every endpoint resolves a dataset through
the registry and calls into app.analysis; no dataset- or disease-specific
logic belongs here. Adding a dataset never touches this file — see
app/registry.py.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Literal

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from app.analysis.compare import expression_by_group, kruskal_wallis, pairwise_tests
from app.analysis.signature import auc_signature_score, log2_mean_signature_score
from app.analysis.survival import binarize_by_median, build_survival_frame, cox_model, kaplan_meier_curves
from app.config import settings
from app.models.contract import Dataset
from app.registry import ensure_loaded, get_descriptor, list_descriptors

ensure_loaded()

app = FastAPI(title="ETP-ALL Expression Explorer API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins,
    allow_methods=["*"],
    allow_headers=["*"],
)


@lru_cache(maxsize=8)
def _get_dataset(dataset_id: str) -> Dataset:
    try:
        descriptor = get_descriptor(dataset_id)
    except KeyError:
        raise HTTPException(404, f"Unknown dataset: {dataset_id}")
    return descriptor.loader()


def _resolve_gene(ds: Dataset, gene: str) -> str:
    """Accept either a symbol or a feature_id; return the feature_id."""
    if gene in ds.matrix.index:
        return gene
    matches = ds.features[ds.features["symbol"] == gene]
    if matches.empty:
        raise HTTPException(404, f"Gene '{gene}' not found in dataset.")
    return matches["feature_id"].iloc[0]


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/datasets")
def list_datasets():
    return {
        "datasets": [
            {
                "dataset_id": d.dataset_id,
                "display_name": d.display_name,
                "group_columns": list(d.group_columns),
                "supports_survival": d.supports_survival,
            }
            for d in list_descriptors()
        ]
    }


@app.get("/datasets/{dataset_id}")
def dataset_info(dataset_id: str):
    descriptor = get_descriptor(dataset_id) if dataset_id in {d.dataset_id for d in list_descriptors()} else None
    if descriptor is None:
        raise HTTPException(404, f"Unknown dataset: {dataset_id}")
    ds = _get_dataset(dataset_id)
    return {
        **ds.source.model_dump(),
        "group_columns": list(descriptor.group_columns),
        "supports_survival": descriptor.supports_survival,
    }


@app.get("/datasets/{dataset_id}/compare")
def compare(dataset_id: str, gene: str, group_column: str):
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    df = expression_by_group(ds.matrix, ds.samples, feature_id, group_column)
    if df.empty:
        raise HTTPException(400, f"No samples with both expression and a value for '{group_column}'.")
    return {
        "gene": gene,
        "group_column": group_column,
        "points": df.to_dict("records"),
        "pairwise_tests": pairwise_tests(df).to_dict("records"),
        "kruskal_wallis": kruskal_wallis(df),
    }


class SignatureRequest(BaseModel):
    genes: list[str]
    method: Literal["auc", "log2_mean"] = "auc"


@app.post("/datasets/{dataset_id}/signature-score")
def signature_score(dataset_id: str, body: SignatureRequest):
    ds = _get_dataset(dataset_id)
    feature_ids = [_resolve_gene(ds, g) for g in body.genes]
    scorer = auc_signature_score if body.method == "auc" else log2_mean_signature_score
    scores = scorer(ds.matrix, feature_ids)
    return {"genes": body.genes, "method": body.method, "scores": scores.to_dict()}


class SurvivalRequest(BaseModel):
    genes: list[str]
    covariates: list[str] = []


@app.post("/datasets/{dataset_id}/survival")
def survival(dataset_id: str, body: SurvivalRequest):
    descriptor = get_descriptor(dataset_id)
    if not descriptor.supports_survival:
        raise HTTPException(
            400,
            f"'{descriptor.display_name}' has no clinical follow-up data — survival analysis isn't available for it.",
        )
    ds = _get_dataset(dataset_id)
    feature_ids = [_resolve_gene(ds, g) for g in body.genes]
    scores = auc_signature_score(ds.matrix, feature_ids)
    surv_df = build_survival_frame(ds.samples, scores)
    if len(surv_df) < 10:
        raise HTTPException(400, "Fewer than 10 complete cases (score + clinical data) for survival analysis.")
    group = binarize_by_median(surv_df)
    result = kaplan_meier_curves(surv_df, group)
    try:
        result["cox"] = cox_model(surv_df, covariates=body.covariates)
    except ValueError:
        result["cox"] = None
    return {"genes": body.genes, "n": len(surv_df), **result}
