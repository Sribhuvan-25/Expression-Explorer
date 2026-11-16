"""
API surface. Thin by design — every endpoint resolves a dataset through
the registry and calls into app.analysis; no dataset- or disease-specific
logic belongs here. Adding a dataset never touches this file — see
app/registry.py.
"""
from __future__ import annotations

from contextlib import asynccontextmanager
from functools import lru_cache
from typing import Literal

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from app.analysis.aux_analysis import (
    expression_by_mutation_status,
    expression_vs_aux_correlation,
    most_mutated_genes,
)
from app.analysis.compare import (
    clean_group_label,
    expression_by_group,
    group_labels,
    kruskal_wallis,
    missing_group_reason,
    pairwise_tests,
)
from app.analysis.contrast import contrast_gene_stats, top_contrast_genes
from app.analysis.correlation import co_expression_network, gene_correlation
from app.analysis.differential import top_differential_genes
from app.analysis.dimensionality import pca
from app.analysis.ranking import rank_by_gene, rank_by_signature
from app.analysis.signature import auc_signature_score, is_log2_unit, log2_mean_signature_score
from app.analysis.scale import on_analysis_scale
from app.analysis.survival import tie_sensitivity, binarize_by_cutoff, build_survival_frame, cox_model, kaplan_meier_curves
from lifelines.exceptions import ConvergenceError, StatError
from app.config import settings
from app.models.contract import MIN_AUX_SAMPLES, AuxLayer, AuxMatrix, Dataset
from app.registry import ensure_loaded, get_descriptor, list_descriptors
from app.services.gene_info import lookup_gene
from app.services.warmup import start_background_warmup, tracker as warmup_tracker

ensure_loaded()


@asynccontextmanager
async def lifespan(_app: FastAPI):
    """Warm every dataset's cache on a background thread at startup.

    Deliberately NOT awaited: the container used to run the same work to
    completion before uvicorn bound, which made deploy success a race
    against a multi-minute download and cost a failed deploy once the aux
    layers pushed cold start to ~848s of a 900s healthcheck budget. See
    app/services/warmup.py for the full history.

    Disabled by WARMUP_ON_STARTUP=0 -- tests and local runs generally
    want the lazy path, not ~500MB of downloads.
    """
    if settings.warmup_on_startup:
        start_background_warmup(
            loader_for=_get_dataset,
            dataset_ids=[d.dataset_id for d in list_descriptors()],
        )
    yield


app = FastAPI(title="Expression Explorer API", lifespan=lifespan)
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins,
    allow_methods=["*"],
    allow_headers=["*"],
)


def _descriptor(dataset_id: str):
    """The dataset's descriptor, or a 404 (a bare get_descriptor raised
    KeyError, and /survival and /group-values answered 500)."""
    try:
        return get_descriptor(dataset_id)
    except KeyError:
        raise HTTPException(404, f"Unknown dataset: {dataset_id}")


@lru_cache(maxsize=8)
def _get_dataset(dataset_id: str) -> Dataset:
    return _descriptor(dataset_id).loader()


def _require_group_column(dataset_id: str, group_column: str) -> None:
    """Only declared grouping columns can group samples: /compare accepted
    days_to_death (about 100 'groups') or sample_id."""
    descriptor = _descriptor(dataset_id)
    if group_column not in descriptor.group_columns:
        raise HTTPException(
            422,
            f"'{group_column}' is not a grouping column on '{descriptor.display_name}' "
            f"(available: {', '.join(descriptor.group_columns) or 'none'}).",
        )


def _sample_groups(dataset_id: str, ds: Dataset) -> dict[str, dict]:
    """sample_id -> {declared group column: value} for colouring points.

    A declared column can live either inside each sample's
    `group_columns` dict or as a top-level samples column -- TARGET keeps
    `etp_status`/`mrd_status`/`vital_status` top-level, which is why
    `expression_by_group` checks both. PCA used to return only the dict,
    so its "Color by etp_status" default coloured no TARGET point at all.
    Missing values come back as None, never NaN (strict JSON).
    """
    declared = get_descriptor(dataset_id).group_columns
    samples = ds.samples
    out: dict[str, dict] = {}
    for sid, row in samples.iterrows():
        nested = row.get("group_columns") if "group_columns" in samples.columns else None
        nested = nested if isinstance(nested, dict) else {}
        groups = {}
        for col in declared:
            v = row[col] if col in samples.columns else nested.get(col)
            groups[col] = clean_group_label(v)
        out[str(sid)] = groups
    return out


def _resolve_gene(ds: Dataset, gene: str) -> str:
    """Accept either a symbol or a feature_id; return the feature_id.

    Whitespace and case are normalised here rather than at each call site:
    `mutated_gene` was being `.strip().upper()`'d by its own endpoint while
    `gene` wasn't, so `?gene=dtx1` 404'd on an API the UI happens to
    uppercase for you. Normalising centrally makes every endpoint behave
    the same way (caught in QA). Ensembl feature_ids are matched before
    upper-casing since they're already canonical.
    """
    if gene in ds.matrix.index:
        return gene
    stripped = gene.strip()
    if stripped in ds.matrix.index:
        return stripped
    gene = stripped.upper()
    if gene in ds.matrix.index:
        return gene
    symbols = ds.features["symbol"]
    matches = ds.features[symbols == stripped]
    if matches.empty:
        matches = ds.features[symbols == gene]
    if matches.empty:
        # Case-insensitive last resort, for symbols that are not all
        # upper-case (C1orf112, the "orf" family): upper-casing the query
        # alone could never match them on an Ensembl-indexed dataset, so
        # `?gene=C1orf112` 404'd on TARGET and GSE324984 while the gene was
        # present. Only an unambiguous hit resolves.
        folded = symbols.astype(str).str.upper() == gene
        hits = sorted(set(symbols[folded]))
        if len(hits) > 1:
            raise HTTPException(
                422, f"'{stripped}' matches several symbols differing only by case ({', '.join(hits)}) -- use the exact symbol."
            )
        matches = ds.features[folded]
    if matches.empty:
        # A legacy symbol this dataset records as an alias (e.g. KIAA1524
        # for CIP2A) is named in the error rather than silently matched:
        # alias *resolution* is still undecided (METHODS.md 7.4), but a bare
        # "not found" for a gene that is present under its current name
        # reads as a coverage claim and is wrong.
        current = _current_symbols_for_alias(ds, gene)
        if current:
            raise HTTPException(
                404,
                f"Gene '{gene}' not found in dataset under that name; this dataset records it "
                f"as an older name for {', '.join(current)} -- search for that instead.",
            )
        count_only = _contrast_only_feature(ds, stripped)
        if count_only is not None:
            raise HTTPException(
                404,
                f"Gene '{stripped}' ({count_only}) is in this dataset's counts but has no value in its "
                f"expression matrix ({ds.source.expression_unit}) -- see the dataset notes. Its "
                "statistics are in the design contrasts (Perturbation Contrasts).",
            )
        current = _external_current_symbol(ds, stripped)
        if current:
            raise HTTPException(
                404,
                f"Gene '{stripped}' not found in dataset under that name; MyGene.info lists it as an "
                f"older name or alias of {current}, which this dataset has -- search for {current} instead.",
            )
        raise HTTPException(404, f"Gene '{gene}' not found in dataset.")
    ids = list(matches["feature_id"])
    in_matrix = [f for f in ids if f in ds.matrix.index]
    if len(in_matrix) > 1:
        # One symbol, several rows (TARGET: 110 symbols over two or more
        # Ensembl ids). The first-listed row was used, and in 42 of them it
        # is not the most expressed: RGS5 got ENSG00000143248 (mean 0.28
        # TPM) over ENSG00000232995 (1.29), DGCR5 an all-zero row. Use the
        # highest-mean row -- the same label-independent rule as microarray
        # probe collapse (METHODS.md 7.1, 7.4a).
        return ds.matrix.loc[in_matrix].mean(axis=1).idxmax()
    return in_matrix[0] if in_matrix else ids[0]


def _resolve_genes(ds: Dataset, genes: list[str]) -> list[str]:
    """Resolve a gene list, each feature once. "MYB" and "MYB (4602)", or
    a symbol and its Ensembl id, are one gene: PCA took a gene given k
    times k^2 times (TARGET MYCN twice: PC1 0.543 vs 0.375 deduplicated)
    and log2-mean weighted it double."""
    return list(dict.fromkeys(_resolve_gene(ds, g) for g in genes))


def _external_current_symbol(ds: Dataset, gene: str) -> str | None:
    """The current symbol for a legacy name or alias, from MyGene.info, if
    this dataset has it -- a hint only, never a silent match (alias
    resolution is undecided, METHODS.md 7.4). DepMap and TARGET record no
    aliases, so KIAA1524 got a bare "not found" there while CIP2A is
    present. Only when the queried name is listed among the record's
    aliases (the lookup's free-text fallback can match loosely), and any
    lookup failure means no hint."""
    try:
        info = lookup_gene(gene)
    except Exception:
        return None
    if not info or not info.get("symbol"):
        return None
    symbol = str(info["symbol"])
    if symbol.upper() == gene.upper():
        return None
    if gene.upper() not in {str(a).upper() for a in info.get("aliases") or []}:
        return None
    if not ds.features["symbol"].astype(str).str.upper().eq(symbol.upper()).any():
        return None
    return symbol


def _contrast_only_feature(ds: Dataset, gene: str) -> str | None:
    """A feature present in a design-contrast table but not in the
    expression matrix (GSE324984: genes past the TPM file's truncation),
    matched by id or by its displayed symbol, case-insensitively and only
    if unambiguous. None if there is no such feature."""
    for contrast in ds.contrasts.values():
        table = contrast.table
        if gene in table.index and gene not in ds.matrix.index:
            return gene
        if "symbol" in table.columns:
            hits = table.index[table["symbol"].astype(str).str.upper() == gene.upper()]
            hits = [h for h in hits if h not in ds.matrix.index]
            if len(hits) == 1:
                return hits[0]
    return None


def _current_symbols_for_alias(ds: Dataset, alias: str) -> list[str]:
    if "aliases" not in ds.features.columns:
        return []
    folded = alias.strip().upper()
    hits = ds.features[ds.features["aliases"].apply(
        lambda a: isinstance(a, (list, tuple)) and folded in {str(x).upper() for x in a}
    )]
    # No filtering on symbol == feature_id: symbol-indexed datasets
    # (GSE66586, GDS4299, DepMap) have that for every row. GSE324984's
    # duplicate-symbol secondary rows carry the shared name as an alias,
    # but that name is always some primary row's symbol, so a query for it
    # resolves above and never reaches this lookup.
    return sorted(set(hits["symbol"]))


def _annotate_strata(ds: Dataset, df, tests: list[dict]) -> str | None:
    """Mark pairwise tests whose two groups sit in different sub-experiments.

    A dataset with design contrasts names the column separating its
    experiments (GSE324984: `cell_line`). A C4-2B arm vs an NCI-H660 arm is
    a valid test of *something*, but mostly of the cell line, not the
    perturbation -- and Compare showed those p-values beside the
    within-line ones with nothing to tell them apart (METHODS.md 11.1).
    Each test gets `crosses_stratum` ({column, group_a, group_b} or None);
    returns a one-line caveat, or None when nothing crosses.
    """
    columns = sorted({c.spec.stratum_column for c in ds.contrasts.values() if c.spec.stratum_column})
    for t in tests:
        t["crosses_stratum"] = None
    crossed = False
    for column in columns:
        if column in ds.samples.columns:
            strata = ds.samples[column]
        else:
            strata = ds.samples["group_columns"].apply(lambda d: d.get(column))
        by_group = (
            df.assign(_stratum=strata.reindex(df["sample_id"]).to_numpy())
            .groupby("group")["_stratum"]
            .agg(lambda v: set(v.dropna()))
        )
        for t in tests:
            sa, sb = by_group.get(t["group_a"], set()), by_group.get(t["group_b"], set())
            if len(sa) == 1 and len(sb) == 1 and sa != sb:
                t["crosses_stratum"] = {"column": column, "group_a": next(iter(sa)), "group_b": next(iter(sb))}
                crossed = True
    if not crossed:
        return None
    return (
        f"Some pairs compare groups from different {' / '.join(columns)} -- those differences "
        "mostly reflect that, not the perturbation. Use the design contrasts for perturbation effects."
    )


def _resolve_aux_feature(aux: AuxMatrix, feature: str) -> str:
    """Same normalisation contract as `_resolve_gene`, for an aux layer's
    own row labels.

    Kept separate rather than folded into `_resolve_gene` because an aux
    index is not always a gene symbol: the drug layer is keyed by Broad
    compound id ("BRD:BRD-K00104122-001-01-9"), so the blind `.upper()`
    that's right for symbols is not safe here. Exact match wins first, and
    the case-insensitive pass is a fallback that only resolves when it is
    unambiguous.

    Without this, one pane disagreed with itself -- `?gene=myb` resolved
    but `?aux_feature=myb` 404'd with "myb not found in the
    crispr_gene_effect layer", which reads as a coverage claim about the
    data rather than a case mismatch, and could talk a user out of a
    perfectly valid analysis (caught in QA).
    """
    if feature in aux.matrix.index:
        return feature
    stripped = feature.strip()
    if stripped in aux.matrix.index:
        return stripped
    folded = stripped.casefold()
    hits = [str(label) for label in aux.matrix.index if str(label).casefold() == folded]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise HTTPException(
            422,
            f"'{feature}' matches {len(hits)} features in this layer "
            f"({', '.join(sorted(hits)[:5])}) -- use the exact identifier.",
        )
    raise HTTPException(404, f"'{feature}' not found in this layer.")


@app.get("/health")
def health():
    """Liveness only -- deliberately independent of dataset readiness.

    This is the platform's deploy healthcheck. Gating it on data being
    downloaded is what turned a slow GDC pull into a failed deployment;
    the process being up and able to route is the thing this endpoint is
    actually asked about. Use /readiness for data state.
    """
    return {"status": "ok"}


@app.get("/readiness")
def readiness():
    """Which datasets are warm, which are still downloading, which failed.

    Split from /health on purpose: a half-warm server is genuinely
    serving (every ready dataset works, and an unwarmed one still
    lazy-loads), so reporting it as unhealthy would be wrong -- but
    pretending it is fully warm would hide a real, visible delay from
    whoever is looking. Both facts get said, separately.
    """
    return warmup_tracker.snapshot()


@app.get("/genes/{symbol}")
def gene_info(symbol: str):
    """Read-only gene annotation (name, summary, aliases, Ensembl/Entrez
    IDs) from MyGene.info -- independent of any registered dataset, since
    the annotation is a property of the gene, not of a dataset's own
    feature table. See app/services/gene_info.py for the alias-resolution
    rationale and why this isn't wired into gene lookup elsewhere yet."""
    symbol = symbol.strip()
    # A blank symbol became an empty cache filename and a 5,000-character
    # one exceeded the filesystem's name limit: both answered 500.
    if not symbol or len(symbol) > 64:
        raise HTTPException(422, "A gene symbol of 1 to 64 characters is required.")
    try:
        info = lookup_gene(symbol)
    except Exception:
        raise HTTPException(503, "The gene annotation service could not be reached; try again later.")
    if info is None:
        raise HTTPException(404, f"No annotation found for '{symbol}'.")
    return info


@app.get("/datasets")
def list_datasets():
    # n_samples/assay/accession come from the loaded dataset's own source
    # record rather than the descriptor, so the listing can show provenance
    # (how many samples, which assay, which accession) without the frontend
    # making a follow-up request per dataset.
    #
    # This endpoint must NEVER block on a download. It is the first call
    # the frontend makes, and on a cold container a dataset can be minutes
    # from ready -- calling the loader here would hang the whole UI behind
    # a ~14-minute GDC pull, which is the failure this background-warmup
    # work exists to remove. So provenance is read only for datasets
    # already warm; the rest are listed with `warming: true` and the UI
    # says so. `_get_dataset` is lru_cached, so a ready dataset is a dict
    # lookup.
    out = []
    warm_states = {d["dataset_id"]: d["state"] for d in warmup_tracker.snapshot()["datasets"]}
    for d in list_descriptors():
        state = warm_states.get(d.dataset_id)
        entry = {
            "dataset_id": d.dataset_id,
            "display_name": d.display_name,
            "group_columns": list(d.group_columns),
            "supports_survival": d.supports_survival,
            # "warming" means: real dataset, not queryable *yet*. Distinct
            # from a dataset that failed to load, which is queryable via
            # the lazy path and will surface its own error.
            "warming": state in ("pending", "loading"),
        }
        if entry["warming"]:
            out.append(entry)
            continue
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
            entry["expression_unit"] = source.expression_unit
            # Lets the UI offer the design-contrast pane only where it
            # applies, without probing every dataset (METHODS.md 11.1).
            entry["n_contrasts"] = len(_get_dataset(d.dataset_id).contrasts)
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
    requested = {d.strip() for d in dataset_ids.split(",") if d.strip()} if dataset_ids else None
    if requested:
        unknown = sorted(requested - {d.dataset_id for d in list_descriptors()})
        if unknown:
            # A typo was silently dropped, so the user saw fewer datasets
            # than they asked for with no reason.
            raise HTTPException(404, f"Unknown dataset(s): {', '.join(unknown)}")
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
        except HTTPException as exc:
            # Pass the resolver's own message through (legacy-symbol hint,
            # count-only gene, case ambiguity) rather than flattening it to
            # "not found" -- the Compare page reads this endpoint, so a
            # flattened reason hid exactly the explanations /compare gives.
            reason = exc.detail
            results.append(
                {"dataset_id": d.dataset_id, "display_name": d.display_name, "skipped": True,
                 "skip_reason": reason}
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
        tests = pairwise_tests(df).to_dict("records")
        stratum_caveat = _annotate_strata(ds, df, tests)
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
            "exclusion_reason": missing_group_reason(
            ds.samples, group_column, set(ds.samples["sample_id"]) - set(df["sample_id"])
        ),
            "points": df.to_dict("records"),
            "pairwise_tests": tests,
            "stratum_caveat": stratum_caveat,
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
    descriptor = _descriptor(dataset_id)
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
    counts = group_labels(ds.samples, group_column).dropna().value_counts()
    return {
        "group_column": group_column,
        "values": [{"value": v, "n": int(n)} for v, n in counts.items()],
        # Distinct from a 404: the column is genuinely declared on this
        # dataset, but the current load has zero non-null values for it.
        # Without this the UI can't tell "typo'd column name" (404),
        # "genuinely empty right now" (this), and "the dropdown just
        # hasn't loaded yet" apart -- all three look identical as a bare
        # []. The realistic cause is an upstream fetch for a
        # supplementary grouping column (e.g. target_all_p2's
        # etp_status/mrd_status, sourced from a third-party file this app
        # doesn't control) degrading on this particular load -- see
        # app/ingest/_http.py and METHODS.md 7.10.
        "unavailable": len(counts) == 0,
    }


@app.get("/datasets/{dataset_id}/compare")
def compare(dataset_id: str, gene: str, group_column: str):
    _require_group_column(dataset_id, group_column)
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
    tests = pairwise_tests(df).to_dict("records")
    stratum_caveat = _annotate_strata(ds, df, tests)
    return {
        "gene": gene,
        "group_column": group_column,
        "n_dataset_total": n_total,
        "n_excluded": n_excluded,
        # None when nothing was excluded -- see compare_multi's identical fix
        # above for why this can't be an unconditional string.
        "exclusion_reason": missing_group_reason(
            ds.samples, group_column, set(ds.samples["sample_id"]) - set(df["sample_id"])
        ),
        "points": df.to_dict("records"),
        "pairwise_tests": tests,
        "stratum_caveat": stratum_caveat,
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
    _require_group_column(dataset_id, group_column)
    if group_a == group_b:
        raise HTTPException(422, f"Group A and group B must differ (both were {group_a!r}).")
    ds = _get_dataset(dataset_id)
    try:
        result = top_differential_genes(ds.matrix, ds.samples, group_column, group_a, group_b, top_n=top_n)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # Sample accounting (METHODS.md 6.1), each reason on its own: samples
    # in neither compared group vs samples with no value at all.
    labels = group_labels(ds.samples, group_column)
    n_total = len(ds.samples)
    n_missing = int(labels.isna().sum())
    n_other = n_total - n_missing - result["n_a"] - result["n_b"]
    breakdown = []
    if n_other:
        breakdown.append({"reason": f"in another '{group_column}' group", "n": n_other})
    if n_missing:
        breakdown.append({"reason": f"with no value recorded for '{group_column}'", "n": n_missing})
    result["n_dataset_total"] = n_total
    result["n_excluded"] = n_total - result["n_a"] - result["n_b"]
    result["exclusion_breakdown"] = breakdown
    result["exclusion_reason"] = "; ".join(f"{e['n']} {e['reason']}" for e in breakdown) or None
    # Attach display symbols the same way /co-expression does -- the
    # matrix's own index may be an Ensembl id (TARGET) or a symbol
    # (DepMap/GDS4299) depending on the dataset, but a reader always wants
    # a human-readable symbol regardless of which id space is underneath.
    symbol_by_feature = ds.features["symbol"]
    for gene in result["genes"]:
        gene["symbol"] = symbol_by_feature.get(gene["feature_id"], gene["feature_id"])
    return result


# --------------------------------------------------------------------------
# Pre-specified design contrasts (METHODS.md 11)
# --------------------------------------------------------------------------

def _require_contrast(ds: Dataset, contrast_id: str):
    if contrast_id not in ds.contrasts:
        available = ", ".join(sorted(ds.contrasts)) or "none"
        raise HTTPException(
            404, f"No contrast '{contrast_id}' on {ds.source.display_name} (available: {available})."
        )
    return ds.contrasts[contrast_id]


@app.get("/datasets/{dataset_id}/contrasts")
def list_contrasts(dataset_id: str):
    """The design contrasts this dataset declares -- empty for cohort
    datasets, which compare whatever groups the user picks."""
    ds = _get_dataset(dataset_id)
    return {
        "dataset_id": dataset_id,
        "expression_unit": ds.source.expression_unit,
        "contrasts": [_contrast_summary(ds, c) for c in ds.contrasts.values()],
    }


class ContrastPanelRequest(BaseModel):
    genes: list[str]
    contrast_ids: list[str] | None = None


def _resolve_panel_gene(ds: Dataset, gene: str) -> tuple[str, str | None]:
    """(feature_id, alias_note). Exact symbol/id first, via `_resolve_gene`.

    Unlike every other endpoint, the panel then accepts a legacy symbol the
    dataset itself records as an alias (KIAA1524 -> CIP2A) -- and says so on
    that row. A researcher's gene list is written in whichever names they
    learned; rejecting the whole 25-gene panel over a rename, or dropping
    the gene, would both be worse than a disclosed resolution. Scoped to
    this endpoint deliberately: METHODS.md 7.4's question (should alias
    matching change *which gene a query hits* everywhere) is still open.
    """
    try:
        return _resolve_gene(ds, gene), None
    except HTTPException:
        current = _current_symbols_for_alias(ds, gene)
        if len(current) == 1:
            # Same row `_resolve_gene` would pick for that symbol.
            fid = _resolve_gene(ds, current[0])
            return fid, f"'{gene.strip()}' is an older name for {current[0]}"
        # A gene with statistics but no expression-matrix row (beyond the
        # truncated TPM file) is still a valid panel row -- its numbers
        # exist; only the TPM column is empty, and the row says why.
        count_only = _contrast_only_feature(ds, gene.strip())
        if count_only is not None:
            return count_only, None
        raise


def _feature_identity(ds: Dataset, feature_id: str, contrasts) -> dict:
    """symbol / gene_name / aliases / shared-name and matrix-presence notes
    for one panel row -- the same identity the top-genes list shows."""
    in_matrix = feature_id in ds.matrix.index
    row = next((c.table.loc[feature_id] for c in contrasts if feature_id in c.table.index), None)
    symbol = ds.features.loc[feature_id, "symbol"] if in_matrix else (row["symbol"] if row is not None else feature_id)
    gene_name = str(row["gene_name"]) if row is not None and "gene_name" in row.index else symbol
    n_loci = int(row["n_loci_with_name"]) if row is not None and "n_loci_with_name" in row.index else 1
    notes = []
    if not in_matrix:
        notes.append(
            f"No {ds.source.expression_unit} value: this gene is absent from the dataset's expression "
            "matrix (see dataset notes); statistics come from the counts."
        )
    if n_loci > 1:
        notes.append(
            f"The name {gene_name} is shared by {n_loci} loci; "
            + (
                "this is the most-expressed one (mean TPM, or mean count where a copy has no TPM), which a search for the name reaches."
                if symbol == gene_name
                else f"this is a secondary locus, shown by its Ensembl id."
            )
            + " Reads can be split arbitrarily between copies (e.g. PAR X/Y), so read its fold-change with care."
        )
    return {
        "symbol": symbol,
        "gene_name": gene_name,
        "aliases": list(ds.features.loc[feature_id, "aliases"] or []) if in_matrix else [],
        "in_expression_matrix": in_matrix,
        "n_loci_with_name": n_loci,
        "identity_notes": notes,
    }


def _contrast_summary(ds: Dataset, contrast) -> dict:
    """Contrast description plus the sample accounting CLAUDE.md requires:
    a contrast runs on 4 of GSE324984's 8 samples, by design."""
    summary = contrast.summary()
    n_total = len(ds.samples)
    used = summary["n_reference"] + summary["n_treatment"]
    spec = contrast.spec
    summary.update(
        n_dataset_total=n_total,
        n_excluded=n_total - used,
        exclusion_reason=(
            f"By design: only the '{spec.reference}' and '{spec.treatment}' samples"
            + (f" ({spec.stratum_column} = {spec.stratum_value})" if spec.stratum_value else "")
            + " enter this contrast; the other samples belong to other arms."
            if n_total > used else None
        ),
    )
    return summary


@app.post("/datasets/{dataset_id}/contrast-panel")
def contrast_panel(dataset_id: str, body: ContrastPanelRequest):
    """A gene list x the dataset's design contrasts: stored log2FC / p / q
    per gene per contrast, with each replicate's expression value.

    Genes that cannot be found are listed in `missing` with the reason,
    not silently dropped and not fatal -- each row is independent, unlike
    a signature score where one missing gene changes everyone's number
    (METHODS.md 2.5, 11.4)."""
    genes = [g for g in (s.strip() for s in body.genes) if g]
    if not genes:
        raise HTTPException(422, "At least one gene is required.")
    if len(genes) > 200:
        raise HTTPException(422, "At most 200 genes per panel.")
    ds = _get_dataset(dataset_id)
    if not ds.contrasts:
        raise HTTPException(
            422, f"{ds.source.display_name} declares no design contrasts; use Compare instead."
        )
    contrast_ids = body.contrast_ids or list(ds.contrasts)
    contrasts = [_require_contrast(ds, cid) for cid in contrast_ids]

    rows, missing, seen = [], [], set()
    for gene in genes:
        try:
            feature_id, alias_note = _resolve_panel_gene(ds, gene)
        except HTTPException as exc:
            missing.append({"gene": gene, "reason": exc.detail})
            continue
        if feature_id in seen:
            missing.append({"gene": gene, "reason": "Duplicate of a gene already in the panel."})
            continue
        seen.add(feature_id)
        rows.append(
            {
                "gene": gene,
                "feature_id": feature_id,
                "alias_note": alias_note,
                **_feature_identity(ds, feature_id, contrasts),
                "contrasts": [contrast_gene_stats(ds.matrix, c, feature_id) for c in contrasts],
            }
        )
    return {
        "dataset_id": dataset_id,
        "expression_unit": ds.source.expression_unit,
        "contrasts": [_contrast_summary(ds, c) for c in contrasts],
        "rows": rows,
        "missing": missing,
    }


@app.get("/datasets/{dataset_id}/contrasts/{contrast_id}/top")
def contrast_top_genes(
    dataset_id: str, contrast_id: str, n: int = 50, direction: Literal["both", "up", "down"] = "both"
):
    if not (1 <= n <= 500):
        raise HTTPException(422, "n must be between 1 and 500.")
    ds = _get_dataset(dataset_id)
    contrast = _require_contrast(ds, contrast_id)
    result = top_contrast_genes(ds.matrix, contrast, n=n, direction=direction)
    result["contrast"] = _contrast_summary(ds, contrast)
    return {"expression_unit": ds.source.expression_unit, **result}


class SignatureRequest(BaseModel):
    genes: list[str]
    method: Literal["auc", "log2_mean"] = "auc"


@app.post("/datasets/{dataset_id}/signature-score")
def signature_score(dataset_id: str, body: SignatureRequest):
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
    ds = _get_dataset(dataset_id)
    feature_ids = _resolve_genes(ds, body.genes)
    if body.method == "auc":
        scores = auc_signature_score(ds.matrix, feature_ids)
    else:
        scores = log2_mean_signature_score(
            ds.matrix, feature_ids, already_log2=is_log2_unit(ds.source.expression_unit)
        )
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
    feature_ids = _resolve_genes(ds, body.genes)
    df = rank_by_signature(
        ds.matrix, ds.samples, feature_ids, method=body.method,
        already_log2=is_log2_unit(ds.source.expression_unit),
    )
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
    pair, value_unit = on_analysis_scale(ds.matrix.loc[[feature_a, feature_b]], ds.source.expression_unit)
    try:
        result = gene_correlation(pair, feature_a, feature_b, method=method)
    except ValueError as exc:
        # Messages name the matrix's feature ids (Ensembl on TARGET); say
        # them the way the user typed them.
        raise HTTPException(422, str(exc).replace(feature_a, gene_a).replace(feature_b, gene_b))
    result["gene_a"] = gene_a
    result["gene_b"] = gene_b
    result["value_unit"] = value_unit
    # Per-point groups so the scatter can be coloured: r pools every
    # sample, and a between-group difference alone can produce a strong r
    # (GSE66586 MYCN vs PPP2R5B, r = -0.85, is mostly amplified vs not).
    groups = _sample_groups(dataset_id, ds)
    for point in result["points"]:
        point["group_columns"] = groups.get(point["sample_id"], {})
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
    try:
        matrix, value_unit = on_analysis_scale(ds.matrix, ds.source.expression_unit)
        network = co_expression_network(matrix, feature_id, top_n=top_n, method=method, direction=direction)
    except ValueError as exc:
        raise HTTPException(422, str(exc).replace(feature_id, gene))
    # Map each result's feature_id back to its display symbol -- the
    # network is computed on the matrix's own index (Ensembl IDs for
    # TARGET, symbols for DepMap/GDS4299), but a reader wants gene symbols
    # regardless of which id space a given dataset happens to key on.
    symbol_by_feature = ds.features["symbol"]
    for entry in network:
        entry["symbol"] = symbol_by_feature.get(entry["gene"], entry["gene"])
    # n and scale, so an r over 8 samples is not read like one over 469.
    return {
        "gene": gene,
        "method": method,
        "direction": direction,
        "n": int(ds.matrix.shape[1]),
        "value_unit": value_unit,
        "network": network,
    }


@app.get("/datasets/{dataset_id}/aux-layers")
def aux_layers(dataset_id: str):
    """What auxiliary measurement layers this dataset carries, and how much
    of the cohort each one covers. Exists so a UI can offer only the
    analyses a given dataset can actually support, and state the coverage
    up front rather than after a query returns a surprising n."""
    ds = _get_dataset(dataset_id)
    out = []
    for layer, aux in ds.aux.items():
        covered, total = ds.aux_sample_overlap(layer)
        shared = [c for c in ds.matrix.columns if c in aux.matrix.columns]
        out.append({
            "layer": layer.value,
            "value_label": aux.value_label,
            "value_description": aux.value_description,
            "source_note": aux.source_note,
            "n_features": aux.n_features,
            # n_features alone oversells a sparse layer -- PRISM ships 6,790
            # compounds but only ~22% were assayed on enough of this
            # dataset's samples to be analysable. See
            # AuxMatrix.n_usable_features.
            "n_usable_features": aux.n_usable_features(shared),
            "n_samples_covered": covered,
            "n_dataset_total": total,
        })
    return {"dataset_id": dataset_id, "layers": out}


@app.get("/datasets/{dataset_id}/aux-features")
def aux_features(
    dataset_id: str,
    layer: Literal["crispr_gene_effect", "drug_sensitivity", "mutation_status"],
    q: str | None = None,
    limit: int = 50,
    usable_only: bool = True,
):
    """Searchable list of features in an aux layer that actually have data
    on THIS dataset's samples, newest-coverage first.

    Exists for the same reason `/mutated-genes` does: a user cannot guess
    a valid value. It matters far more here -- PRISM ships 6,790 compounds
    but 5,272 of them (77.6%) have zero measurements on the 79 covered
    Lymphoid lines, so free-texting a compound id fails ~78% of the time,
    and the ids are opaque Broad strings ("BRD:BRD-K05804044-001-18-5")
    that nobody types from memory. Results carry the human drug name and
    mechanism where the layer has a feature table, so the picker can be
    searched by name rather than id.

    `usable_only` (default) hides features with too few values to
    correlate at all, which is the difference between a 6,790-entry list
    that mostly 422s and a 1,518-entry list that works.
    """
    if not (1 <= limit <= 200):
        raise HTTPException(422, "limit must be between 1 and 200.")
    ds = _get_dataset(dataset_id)
    try:
        aux = ds.require_aux(AuxLayer(layer))
    except ValueError as exc:
        raise HTTPException(404, str(exc))

    shared = [c for c in ds.matrix.columns if c in aux.matrix.columns]
    if not shared:
        return {"dataset_id": dataset_id, "layer": layer, "n_matching": 0, "features": []}

    counts = aux.usable_counts(shared)
    if usable_only:
        counts = counts[counts >= MIN_AUX_SAMPLES]

    # Search over the readable name when the layer has one (drug names),
    # falling back to the raw id -- searching only ids would defeat the
    # point for exactly the layer that needs this most.
    labels = {}
    if aux.features is not None and "symbol" in aux.features.columns:
        labels = dict(zip(aux.features["feature_id"], aux.features["symbol"]))
    moas = {}
    if aux.features is not None and "moa" in aux.features.columns:
        moas = dict(zip(aux.features["feature_id"], aux.features["moa"]))
    targets = {}
    if aux.features is not None and "target" in aux.features.columns:
        targets = dict(zip(aux.features["feature_id"], aux.features["target"]))

    if q:
        # Search name, id, mechanism AND target: a researcher looks for
        # "MEK inhibitor" or "BRAF", not a drug's trade name, and matching
        # only name+id returned nothing for exactly those queries.
        needle = q.strip().casefold()
        counts = counts[
            [
                needle in str(fid).casefold()
                or needle in str(labels.get(fid, "")).casefold()
                or needle in str(moas.get(fid, "")).casefold()
                or needle in str(targets.get(fid, "")).casefold()
                for fid in counts.index
            ]
        ]

    n_matching = int(len(counts))
    top = counts.sort_values(ascending=False).head(limit)

    return {
        "dataset_id": dataset_id,
        "layer": layer,
        "n_matching": n_matching,
        "n_samples_covered": len(shared),
        "features": [
            {
                "feature_id": str(fid),
                "label": str(labels.get(fid, fid)) or str(fid),
                "moa": (str(moas.get(fid, "")) or None),
                "n_samples_with_value": int(n),
            }
            for fid, n in top.items()
        ],
    }


@app.get("/datasets/{dataset_id}/mutated-genes")
def mutated_genes(dataset_id: str, top_n: int = 25):
    """Most recurrently mutated genes in this cohort -- the entry point for
    picking something to stratify by, since a user can't guess which genes
    are actually mutated often enough to compare."""
    if not (1 <= top_n <= 200):
        raise HTTPException(422, "top_n must be between 1 and 200.")
    ds = _get_dataset(dataset_id)
    try:
        return most_mutated_genes(ds, top_n=top_n)
    except ValueError as exc:
        raise HTTPException(422, str(exc))


def _resolve_mutated_gene(ds: Dataset, gene: str) -> str:
    """Match a mutation-layer gene exactly, else case-insensitively if
    unambiguous. Upper-casing alone turned C11orf53 (mutated in 2 covered
    TARGET samples) into "C11ORF53 has no recorded mutation" -- 40 layer
    genes have non-upper-case symbols."""
    stripped = gene.strip()
    aux = ds.aux.get(AuxLayer.MUTATION_STATUS) if ds.aux else None
    if aux is None or stripped in aux.matrix.index:
        return stripped
    hits = [g for g in aux.matrix.index if str(g).upper() == stripped.upper()]
    return hits[0] if len(hits) == 1 else stripped.upper()


@app.get("/datasets/{dataset_id}/compare-by-mutation")
def compare_by_mutation(dataset_id: str, gene: str, mutated_gene: str):
    """Expression of `gene`, split by whether each sample carries a
    non-silent mutation in `mutated_gene` -- GEPIA3's "Hotspot Mutation"
    equivalent. Deliberately a separate endpoint from /compare rather than
    a mode flag: the grouping comes from a different measurement layer with
    its own partial coverage, which needs its own exclusion accounting."""
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    try:
        result = expression_by_mutation_status(ds, feature_id, _resolve_mutated_gene(ds, mutated_gene))
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    result["gene"] = gene
    result["expression_unit"] = ds.source.expression_unit
    return result


@app.get("/datasets/{dataset_id}/expression-vs-aux")
def expression_vs_aux(
    dataset_id: str,
    gene: str,
    layer: Literal["crispr_gene_effect", "drug_sensitivity"],
    aux_feature: str,
    method: Literal["pearson", "spearman"] = "pearson",
):
    """Correlate a gene's expression against an auxiliary measurement on the
    same samples -- e.g. does MYCN expression track with dependency on a
    gene (CRISPR) or with sensitivity to a compound (PRISM drug screen)."""
    ds = _get_dataset(dataset_id)
    feature_id = _resolve_gene(ds, gene)
    try:
        aux_layer = ds.require_aux(AuxLayer(layer))
    except ValueError as exc:
        raise HTTPException(404, str(exc))
    aux_feature = _resolve_aux_feature(aux_layer, aux_feature)
    try:
        result = expression_vs_aux_correlation(
            ds, AuxLayer(layer), feature_id, aux_feature, method=method
        )
    except KeyError as exc:
        raise HTTPException(404, str(exc).strip("'"))
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    result["gene"] = gene
    result["aux_feature"] = aux_feature
    # Compound ids ("BRD:BRD-K00104122-001-01-9") mean nothing to a reader,
    # so resolve to the drug name where the layer carries a feature table.
    aux = ds.require_aux(AuxLayer(layer))
    if aux.features is not None and aux_feature in aux.features.index:
        row = aux.features.loc[aux_feature]
        result["aux_feature_label"] = str(row.get("symbol", aux_feature))
        result["aux_feature_moa"] = str(row.get("moa", "")) or None
    else:
        result["aux_feature_label"] = aux_feature
        result["aux_feature_moa"] = None
    return result


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
    missing_reasons = {}
    for g in body.genes:
        try:
            feature_ids.append(_resolve_gene(ds, g))
        except HTTPException as exc:
            missing_input_genes.append(g)
            missing_reasons[g] = exc.detail
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
    feature_ids = list(dict.fromkeys(feature_ids))
    subset, value_unit = on_analysis_scale(ds.matrix.loc[feature_ids], ds.source.expression_unit)
    try:
        result = pca(subset, feature_ids, ds.samples, n_components=body.n_components)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # genes_missing from pca() reflects feature_ids that resolved but
    # weren't in the matrix -- shouldn't normally happen since
    # _resolve_gene already checked the matrix, but a gene symbol that
    # never resolved at all (unrecognised input) is a distinct kind of
    # "missing" a caller needs to see, so both lists are surfaced rather
    # than the input-level failures being silently swallowed.
    groups = _sample_groups(dataset_id, ds)
    for point in result["points"]:
        point["group_columns"] = groups.get(point["sample_id"], {})
    # Feature ids are Ensembl on TARGET; the title, export name and
    # zero-variance note read them, so give symbols.
    symbol = ds.features.set_index("feature_id")["symbol"] if "feature_id" in ds.features.columns else ds.features["symbol"]
    for key in ("genes_used", "genes_missing", "genes_zero_variance"):
        result[key] = [str(symbol.get(f, f)) for f in result[key]]
    result["value_unit"] = value_unit
    result["genes_unrecognized"] = missing_input_genes
    # Why each one failed, in the resolver's words (e.g. "older name for
    # CIP2A"), so PCA explains what /compare already explains.
    result["genes_unrecognized_reasons"] = missing_reasons
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
    # What to do when a percentile threshold lands exactly on a block of
    # tied scores -- see binarize_by_cutoff. Exposed rather than decided
    # internally because the right answer is a domain judgement: "trim"
    # keeps the arms the caption promises, "exclude" never splits equal
    # scores across arms but can empty an arm entirely (MEF2C/TAL1 on
    # TARGET both give LOW n=0), and "inclusive" reproduces the previous
    # behaviour. The response's `tie_note` says which was applied and how
    # many samples were tied.
    tie_policy: Literal["trim", "exclude", "inclusive"] = "trim"


@app.post("/datasets/{dataset_id}/survival")
def survival(dataset_id: str, body: SurvivalRequest):
    if not body.genes:
        raise HTTPException(422, "At least one gene is required.")
    descriptor = _descriptor(dataset_id)
    if not descriptor.supports_survival:
        raise HTTPException(
            400,
            f"'{descriptor.display_name}' has no clinical follow-up data — survival analysis isn't available for it.",
        )
    ds = _get_dataset(dataset_id)
    feature_ids = _resolve_genes(ds, body.genes)
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
            tie_policy=body.tie_policy,
        )
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    # Report the tie situation regardless of which policy was applied: the
    # user is a domain expert who can only choose sensibly if the app says
    # what it actually did. `inclusive` silently produced arms that weren't
    # the quartiles the caption promised (MEF2C: LOW n=217 of 466), so the
    # fix isn't only picking a better default -- it's surfacing the tie.
    tie_note = None
    if body.cutoff_method != "median":
        scores = surv_df["signature_score"]
        high_pct, low_pct = (25.0, 25.0) if body.cutoff_method == "quartile" else (
            body.cutoff_high_pct, body.cutoff_low_pct,
        )
        n_tied_high = int((scores == scores.quantile(1 - high_pct / 100)).sum())
        n_tied_low = int((scores == scores.quantile(low_pct / 100)).sum())
        n_tied = max(n_tied_high, n_tied_low)
        if n_tied > 1:
            tie_note = {
                "n_tied_at_cutoff": n_tied,
                "tie_policy": body.tie_policy,
                "message": (
                    f"{n_tied} samples share the exact score at a cutoff boundary "
                    f"(common for a single gene: every sample where it falls outside "
                    f"the top 25% of genes scores exactly 0). Applied tie policy: "
                    f"'{body.tie_policy}'."
                    + (
                        " Trim keeps a fixed pseudo-random subset of the tied samples "
                        "to reach the requested arm size."
                        if body.tie_policy == "trim"
                        else ""
                    )
                ),
            }
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
    # Each reason counted on its own (METHODS.md 6.1, 8.6). One merged
    # "235 excluded (no follow-up time recorded, or score fell in the
    # excluded middle band)" put 3 patients without follow-up and 232 cut
    # by the quartile split behind the same word.
    exclusion_breakdown = []
    n_no_followup = n_total - len(surv_df)
    if n_no_followup > 0:
        exclusion_breakdown.append({"reason": "with no follow-up time recorded", "n": n_no_followup})
    if body.cutoff_method != "median":
        high_pct, low_pct = (25.0, 25.0) if body.cutoff_method == "quartile" else (
            body.cutoff_high_pct, body.cutoff_low_pct,
        )
        scores = surv_df["signature_score"]
        in_band = (scores < scores.quantile(1 - high_pct / 100)) & (scores > scores.quantile(low_pct / 100))
        n_band = int((in_band & group.isna()).sum())
        n_tie_dropped = int((~in_band & group.isna()).sum())
        if n_band > 0:
            exclusion_breakdown.append({
                "reason": f"with a score between the low ({low_pct:g}%) and high ({high_pct:g}%) cutoffs",
                "n": n_band,
            })
        if n_tie_dropped > 0:
            exclusion_breakdown.append({
                "reason": (
                    f"with a score tied at the low ({low_pct:g}%) or high ({high_pct:g}%) cutoff, "
                    f"dropped by the '{body.tie_policy}' tie policy"
                ),
                "n": n_tie_dropped,
            })
    exclusion_reason = "; ".join(f"{e['n']} {e['reason']}" for e in exclusion_breakdown) or None
    # When trim had to choose which tied patients to keep, say how much
    # that arbitrary choice moves the answer (METHODS.md 3.1a).
    if tie_note is not None and body.tie_policy == "trim" and any(
        "tie policy" in e["reason"] for e in exclusion_breakdown
    ):
        tie_note["sensitivity"] = tie_sensitivity(
            surv_df,
            method=body.cutoff_method,
            custom_high_pct=body.cutoff_high_pct,
            custom_low_pct=body.cutoff_low_pct,
        )
    return {
        "genes": body.genes,
        "n": n_grouped,
        "n_dataset_total": n_total,
        "n_excluded": n_excluded,
        "exclusion_reason": exclusion_reason,
        "exclusion_breakdown": exclusion_breakdown,
        "tie_note": tie_note,
        **result,
    }
