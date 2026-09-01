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
from app.analysis.correlation import co_expression_network, gene_correlation
from app.analysis.differential import top_differential_genes
from app.analysis.dimensionality import pca
from app.analysis.ranking import rank_by_gene, rank_by_signature
from app.analysis.signature import auc_signature_score, log2_mean_signature_score
from app.analysis.survival import binarize_by_cutoff, build_survival_frame, cox_model, kaplan_meier_curves
from lifelines.exceptions import ConvergenceError, StatError
from app.config import settings
from app.models.contract import Dataset
from app.registry import ensure_loaded, get_descriptor, list_descriptors
from app.services.gene_info import lookup_gene

ensure_loaded()

app = FastAPI(title="Expression Explorer API")
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


@app.get("/genes/{symbol}")
def gene_info(symbol: str):
    """Read-only gene annotation (name, summary, aliases, Ensembl/Entrez
    IDs) from MyGene.info -- independent of any registered dataset, since
    the annotation is a property of the gene, not of a dataset's own
    feature table. See app/services/gene_info.py for the alias-resolution
    rationale and why this isn't wired into gene lookup elsewhere yet."""
    info = lookup_gene(symbol)
    if info is None:
        raise HTTPException(404, f"No annotation found for '{symbol}'.")
    return info


@app.get("/datasets")
def list_datasets():
    # n_samples/assay/accession come from the loaded dataset's own source
    # record rather than the descriptor, so the listing can show provenance
    # (how many samples, which assay, which accession) without the frontend
    # making a follow-up request per dataset. _get_dataset is lru_cached and
    # every dataset is pre-warmed at startup, so this stays cheap.
    out = []
    for d in list_descriptors():
        entry = {
            "dataset_id": d.dataset_id,
            "display_name": d.display_name,
            "group_columns": list(d.group_columns),
            "supports_survival": d.supports_survival,
        }
        try:
            source = _get_dataset(d.dataset_id).source
            entry["n_samples"] = source.n_samples
            entry["assay_type"] = source.assay_type
            entry["accession"] = source.accession
            # Provenance caveats (e.g. DepMap's release-staleness disclosure,
            # see METHODS.md 7.5) live in `notes` and are otherwise invisible
            # -- nothing in the frontend called the per-dataset /datasets/{id}
            # endpoint this would have come from, so a real, correctly-computed
            # warning sat in the API response with no UI ever reading it.
            entry["notes"] = source.notes
        except Exception:
            # A dataset that fails to load should still be listed (so the
            # UI can show it exists) rather than taking down the whole
            # registry listing.
            pass
        out.append(entry)
    return {"datasets": out}


@app.get("/datasets/compare-multi")
def compare_multi(gene: str, group_column: str, dataset_ids: str | None = None):
    """Same comparison as /datasets/{id}/compare, run independently across
    several datasets and returned side by side -- one block per dataset,
    each with its own points/stats/n-accounting. Never pools raw values
    across datasets (see METHODS.md 7.2): different assay types and units
    aren't on a comparable scale, so each dataset keeps its own axis.

    `dataset_ids` is a comma-separated allowlist; omit it to run on every
    registered dataset that exposes `group_column` (the "use everything
    compatible" default). A dataset missing the column, or missing the
    gene entirely, is reported as skipped rather than raising -- the
    point of asking for "all compatible datasets" is that incompatible
    ones are quietly left out, not an error per dataset.

    Registered ABOVE /datasets/{dataset_id} deliberately: FastAPI matches
    routes in registration order, so this static path must be declared
    before the dynamic one or "/datasets/compare-multi" would be captured
    as dataset_id="compare-multi" and 404 as an unknown dataset.
    """
    requested = set(dataset_ids.split(",")) if dataset_ids else None
    results = []
    for d in list_descriptors():
        if requested is not None and d.dataset_id not in requested:
            continue
        if group_column not in d.group_columns:
            results.append(
                {"dataset_id": d.dataset_id, "display_name": d.display_name, "skipped": True,
                 "skip_reason": f"does not expose grouping column '{group_column}'"}
            )
            continue
        ds = _get_dataset(d.dataset_id)
        try:
            feature_id = _resolve_gene(ds, gene)
        except HTTPException:
            results.append(
                {"dataset_id": d.dataset_id, "display_name": d.display_name, "skipped": True,
                 "skip_reason": f"gene '{gene}' not found in this dataset"}
            )
            continue
        df = expression_by_group(ds.matrix, ds.samples, feature_id, group_column)
        if df.empty:
            results.append(
                {"dataset_id": d.dataset_id, "display_name": d.display_name, "skipped": True,
                 "skip_reason": f"no samples with both expression and a value for '{group_column}'"}
            )
            continue
        n_total = len(ds.samples)
        n_excluded = n_total - len(df)
        results.append({
            "dataset_id": d.dataset_id,
            "display_name": d.display_name,
            "skipped": False,
            "assay_type": ds.source.assay_type,
            "expression_unit": ds.source.expression_unit,
            "n_dataset_total": n_total,
            "n_excluded": n_excluded,
            # None (not a reason string) when nothing was actually excluded --
            # this used to be unconditional, so e.g. DepMap grouped by
            # lineage (0 excluded, 100% coverage) still claimed "no value
            # recorded for 'lineage'", implying an exclusion that never
            # happened. Caught in QA, not by a user, but a real defect.
            "exclusion_reason": f"no value recorded for '{group_column}'" if n_excluded > 0 else None,
            "points": df.to_dict("records"),
            "pairwise_tests": pairwise_tests(df).to_dict("records"),
            "kruskal_wallis": kruskal_wallis(df),
        })
    return {"gene": gene, "group_column": group_column, "datasets": results}


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


@app.get("/datasets/{dataset_id}/group-values")
def group_values(dataset_id: str, group_column: str):
    """Distinct values a grouping column actually has for this dataset,
    with counts -- e.g. etp_status -> {"ETP": 19, "non-ETP": 146,
    "near-ETP": 25}. Exists specifically so a UI can populate a "pick two
    groups to compare" dropdown without guessing at values or requiring a
    gene to already be entered (an earlier draft of the Differential
    Genes UI probed this via /compare with a hardcoded gene name, which
    is fragile -- not every dataset is guaranteed to have that gene)."""
    descriptor = get_descriptor(dataset_id)
    if group_column not in descriptor.group_columns:
        # Sibling endpoints (/compare -> 400, /differential -> 422) both
        # reject an unrecognised group_column outright; this endpoint used
        # to silently return {"values": []} instead, which looks
        # identical to "this column legitimately has zero values" -- a UI
        # can't tell "you typo'd the column name" from "this column is
        # empty" without an explicit error (caught in QA).
        raise HTTPException(
            404, f"'{group_column}' is not a grouping column on '{descriptor.display_name}'."
        )
    ds = _get_dataset(dataset_id)
    if group_column in ds.samples.columns:
        values = ds.samples[group_column]
    else:
        values = ds.samples["group_columns"].apply(lambda d: d.get(group_column))
    counts = values.dropna().value_counts()
    return {"group_column": group_column, "values": [{"value": v, "n": int(n)} for v, n in counts.items()]}


@app.get("/datasets/{dataset_id}/compare")
def compare(dataset_id: str, gene: str, group_column: str):
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    df = expression_by_group(ds.matrix, ds.samples, feature_id, group_column)
    if df.empty:
        raise HTTPException(400, f"No samples with both expression and a value for '{group_column}'.")
    # Grouping columns are often populated for only part of a cohort --
    # TARGET-ALL-P2 classifies etp_status for 190 of its 469 samples, so a
    # comparison "on TARGET" is silently running on 40% of it. That's a
    # material caveat for anyone reading the p-values, so report what was
    # left out rather than only what was included.
    n_total = len(ds.samples)
    n_excluded = n_total - len(df)
    return {
        "gene": gene,
        "group_column": group_column,
        "n_dataset_total": n_total,
        "n_excluded": n_excluded,
        # None when nothing was excluded -- see compare_multi's identical fix
        # above for why this can't be an unconditional string.
        "exclusion_reason": f"no value recorded for '{group_column}'" if n_excluded > 0 else None,
        "points": df.to_dict("records"),
        "pairwise_tests": pairwise_tests(df).to_dict("records"),
        "kruskal_wallis": kruskal_wallis(df),
    }


@app.get("/datasets/{dataset_id}/differential")
def differential(dataset_id: str, group_column: str, group_a: str, group_b: str, top_n: int = 50):
    """Genome-wide: every gene in this dataset ranked by how strongly it
    separates group_a from group_b -- the inverse of /compare (which
    starts from one gene and asks about groups; this starts from a
    grouping and asks about every gene). See
    app/analysis/differential.py for why this is two-group-only and why
    it's computed as one vectorised scipy call rather than a per-gene
    loop (a real, measured ~15x difference at TARGET's actual gene
    count, not a micro-optimisation)."""
    if not (1 <= top_n <= 500):
        raise HTTPException(422, "top_n must be between 1 and 500.")
    ds = _get_dataset(dataset_id)
    try:
        result = top_differential_genes(ds.matrix, ds.samples, group_column, group_a, group_b, top_n=top_n)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # Attach display symbols the same way /co-expression does -- the
    # matrix's own index may be an Ensembl id (TARGET) or a symbol
    # (DepMap/GDS4299) depending on the dataset, but a reader always wants
    # a human-readable symbol regardless of which id space is underneath.
    symbol_by_feature = ds.features["symbol"]
    for gene in result["genes"]:
        gene["symbol"] = symbol_by_feature.get(gene["feature_id"], gene["feature_id"])
    return result


class SignatureRequest(BaseModel):
    genes: list[str]
    method: Literal["auc", "log2_mean"] = "auc"


@app.post("/datasets/{dataset_id}/signature-score")
def signature_score(dataset_id: str, body: SignatureRequest):
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
    ds = _get_dataset(dataset_id)
    feature_ids = [_resolve_gene(ds, g) for g in body.genes]
    scorer = auc_signature_score if body.method == "auc" else log2_mean_signature_score
    scores = scorer(ds.matrix, feature_ids)
    return {"genes": body.genes, "method": body.method, "scores": scores.to_dict()}


@app.get("/datasets/{dataset_id}/rank")
def rank(dataset_id: str, gene: str):
    """Every sample sorted high to low by one gene's expression, with
    sample metadata attached — e.g. 'which DepMap cell lines have the
    highest MYCN expression', for picking cell lines to purchase."""
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    df = rank_by_gene(ds.matrix, ds.samples, feature_id)
    return {"gene": gene, "n": len(df), "rows": df.to_dict("records")}


class RankSignatureRequest(BaseModel):
    genes: list[str]
    method: Literal["auc", "log2_mean"] = "auc"


@app.post("/datasets/{dataset_id}/rank-signature")
def rank_signature(dataset_id: str, body: RankSignatureRequest):
    """Same as /rank, but ordered by a multi-gene signature score."""
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
    ds = _get_dataset(dataset_id)
    feature_ids = [_resolve_gene(ds, g) for g in body.genes]
    df = rank_by_signature(ds.matrix, ds.samples, feature_ids, method=body.method)
    return {"genes": body.genes, "method": body.method, "n": len(df), "rows": df.to_dict("records")}


@app.get("/datasets/{dataset_id}/correlation")
def correlation(dataset_id: str, gene_a: str, gene_b: str, method: Literal["pearson", "spearman", "kendall"] = "pearson"):
    """Correlation between two genes' expression across this dataset's
    samples -- e.g. does MYCN track with a specific transcription factor.
    Never pools across datasets (matches the Compare page's own rule,
    METHODS.md 7.2): a caller comparing two genes across cohorts runs this
    once per dataset and reads the two results side by side."""
    ds = _get_dataset(dataset_id)
    feature_a = _resolve_gene(ds, gene_a)
    feature_b = _resolve_gene(ds, gene_b)
    try:
        result = gene_correlation(ds.matrix, feature_a, feature_b, method=method)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    result["gene_a"] = gene_a
    result["gene_b"] = gene_b
    return result


@app.get("/datasets/{dataset_id}/co-expression")
def co_expression(
    dataset_id: str,
    gene: str,
    top_n: int = 20,
    method: Literal["pearson", "spearman"] = "pearson",
    direction: Literal["positive", "negative"] = "positive",
):
    """Top-N genes most (anti-)correlated with `gene` across this
    dataset's samples -- a co-expression network rooted at one gene."""
    if not (1 <= top_n <= 200):
        raise HTTPException(422, "top_n must be between 1 and 200.")
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    network = co_expression_network(ds.matrix, feature_id, top_n=top_n, method=method, direction=direction)
    # Map each result's feature_id back to its display symbol -- the
    # network is computed on the matrix's own index (Ensembl IDs for
    # TARGET, symbols for DepMap/GDS4299), but a reader wants gene symbols
    # regardless of which id space a given dataset happens to key on.
    symbol_by_feature = ds.features["symbol"]
    for entry in network:
        entry["symbol"] = symbol_by_feature.get(entry["gene"], entry["gene"])
    return {"gene": gene, "method": method, "direction": direction, "network": network}


class PCARequest(BaseModel):
    genes: list[str]
    n_components: int = 2


@app.post("/datasets/{dataset_id}/pca")
def dataset_pca(dataset_id: str, body: PCARequest):
    """PCA of this dataset's samples over a user-supplied gene set --
    e.g. does ETP vs non-ETP actually separate on a signature's genes
    across more than one axis, not just the single-gene view Compare
    shows. No auto-selected default gene set (deliberate -- see
    METHODS.md): the caller always supplies the genes."""
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
    if not (1 <= body.n_components <= 10):
        raise HTTPException(422, "n_components must be between 1 and 10.")
    ds = _get_dataset(dataset_id)
    feature_ids = []
    missing_input_genes = []
    for g in body.genes:
        try:
            feature_ids.append(_resolve_gene(ds, g))
        except HTTPException:
            missing_input_genes.append(g)
    # If too few genes resolved at all, say so in terms of what the user
    # actually typed (e.g. "1 of 2 genes you entered"), not pca()'s own
    # "1 of 1" -- pca() only ever sees the genes that already resolved, so
    # its own error message has no way to know 2 were originally
    # requested and would otherwise misreport the count back to the user.
    if len(feature_ids) < 2:
        raise HTTPException(
            422,
            f"Only {len(feature_ids)} of {len(body.genes)} requested genes are present in this "
            f"dataset -- need at least 2 to run PCA."
            + (f" Not found: {', '.join(missing_input_genes)}." if missing_input_genes else ""),
        )
    try:
        result = pca(ds.matrix, feature_ids, ds.samples, n_components=body.n_components)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # genes_missing from pca() reflects feature_ids that resolved but
    # weren't in the matrix -- shouldn't normally happen since
    # _resolve_gene already checked the matrix, but a gene symbol that
    # never resolved at all (unrecognised input) is a distinct kind of
    # "missing" a caller needs to see, so both lists are surfaced rather
    # than the input-level failures being silently swallowed.
    result["genes_unrecognized"] = missing_input_genes
    return result


class SurvivalRequest(BaseModel):
    genes: list[str]
    covariates: list[str] = []
    # "quartile"/"custom" deliberately drop a middle band of samples (see
    # binarize_by_cutoff) -- a materially different comparison from
    # "median", not just a different threshold, so it's exposed as its
    # own explicit choice rather than a numeric-only parameter.
    cutoff_method: Literal["median", "quartile", "custom"] = "median"
    cutoff_high_pct: float = 50.0
    cutoff_low_pct: float = 50.0


@app.post("/datasets/{dataset_id}/survival")
def survival(dataset_id: str, body: SurvivalRequest):
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
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
    # Survival runs on fewer samples than the dataset holds: a patient with
    # neither days_to_death nor days_to_last_follow_up has no time axis to
    # be plotted on and is dropped by build_survival_frame. Reporting only
    # the surviving count made that look like a different dataset entirely
    # -- a reviewer sees 469 samples listed but "n = 466" on the curve and
    # reasonably asks which cohort this is. Return the accounting so the UI
    # can say so outright.
    n_total = len(ds.samples)
    try:
        group = binarize_by_cutoff(
            surv_df,
            method=body.cutoff_method,
            custom_high_pct=body.cutoff_high_pct,
            custom_low_pct=body.cutoff_low_pct,
        )
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # Quartile/custom cutoffs exclude a middle band on top of whatever
    # build_survival_frame already dropped -- both are real, distinct
    # reasons a sample doesn't appear on the curve, so n_excluded/
    # exclusion_reason must account for the union of both, not just the
    # clinical-data gap alone (that was the whole point of the sample-
    # accounting work in METHODS.md 6.1; a cutoff mode that silently
    # under-reported exclusions would be regressing exactly that fix).
    n_grouped = int(group.notna().sum())
    n_excluded = n_total - n_grouped
    grouped_df = surv_df.loc[group.notna()]
    grouped_group = group.dropna()
    result = kaplan_meier_curves(grouped_df, grouped_group)
    try:
        result["cox"] = cox_model(grouped_df, covariates=body.covariates)
    # Too few complete cases or an unknown covariate name (ValueError, both
    # raised deliberately by cox_model), and a degenerate fit -- a sparsely
    # populated or collinear covariate -- that lifelines itself rejects
    # (ConvergenceError/StatError) all mean the same thing to a caller: the
    # KM curves above are still valid, but no Cox model could be fit for
    # this covariate choice. Previously only ValueError was caught, so a
    # sparse covariate (e.g. etp_status, populated for a subset of TARGET
    # samples) surfaced as an unhandled 500 instead of cox: null.
    except (ValueError, ConvergenceError, StatError):
        result["cox"] = None
    exclusion_reason = None
    if n_excluded > 0:
        if body.cutoff_method == "median":
            exclusion_reason = "no follow-up time recorded"
        else:
            # quartile fixes both percentiles at 25 regardless of what the
            # request body's cutoff_high_pct/cutoff_low_pct happened to be
            # (those fields only take effect for "custom") -- the message
            # must report what binarize_by_cutoff actually used, not the
            # unused request defaults, or a quartile run would misreport
            # itself as a "50%/50%" cutoff.
            high_pct, low_pct = (25.0, 25.0) if body.cutoff_method == "quartile" else (
                body.cutoff_high_pct, body.cutoff_low_pct,
            )
            exclusion_reason = (
                "no follow-up time recorded, or score fell in the excluded middle band "
                f"between the low ({low_pct:.0f}%) and high ({high_pct:.0f}%) cutoffs"
            )
    return {
        "genes": body.genes,
        "n": n_grouped,
        "n_dataset_total": n_total,
        "n_excluded": n_excluded,
        "exclusion_reason": exclusion_reason,
        **result,
    }
