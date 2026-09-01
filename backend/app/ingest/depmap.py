"""
DepMap loader — cell-line RNA-seq, for the "which lines to buy" question
(PDF p.3: MYCN/PP2A expression across cell lines).

DepMap's own /portal/api/download/files endpoint is behind Cloudflare bot
verification (confirmed directly: a plain request gets a Turnstile
challenge page, not data) and cannot be scraped programmatically. The
legitimate, scriptable distribution channel is the versioned Figshare
release; we resolve the current release by searching Figshare's public API
rather than hardcoding a release id, so a new quarterly release doesn't
silently break this loader.

Known, bounded limitation (see METHODS.md 7.5): DepMap stopped publishing
full quarterly bundles to Figshare after 24Q4 (Dec 2024) -- newer releases
(25Q2 onward) are portal-only, behind the same bot gate. This loader can
therefore only ever find 24Q4 on Figshare going forward, and has no way to
detect a newer release existing elsewhere. Rather than silently serve an
ever-staler snapshot under a "current data" claim, `_release_age_warning`
below computes the release's age from its own title and surfaces it -- in
the dataset's `notes` (so the UI shows it) and as a log warning past a
threshold -- so staleness is visible instead of assumed away.
"""
from __future__ import annotations

import logging
import re
from datetime import date

import httpx
import pandas as pd

from app.config import settings
from app.models.contract import (
    AssayType,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)

FIGSHARE_SEARCH = "https://api.figshare.com/v2/articles/search"
CACHE_DIR = settings.cache_dir / "depmap"

EXPRESSION_FILENAME = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL_FILENAME = "Model.csv"

# DepMap release articles are titled "DepMap <YY>Q<N> Public" (quarter comes
# before "Public"). Figshare's free-text search is relevance-ranked and not
# exhaustive for a single query, so we page through results and keep every
# title matching this pattern, then take the most recently published one.
_RELEASE_TITLE = re.compile(r"^DepMap (\d\d)Q(\d) Public$")

log = logging.getLogger(__name__)

# A quarterly cadence release older than this is flagged, in both the
# dataset's own `notes` and the app log -- verified empirically (2026-08-29)
# that DepMap 24Q4 is the newest release matching this title on Figshare;
# every release since (25Q1 onward) is portal-only, behind the same
# Cloudflare bot gate the module docstring describes. Two quarters is a
# deliberately loose threshold: Figshare lagging DepMap's own release
# cadence by a quarter isn't unusual on its own, but anything past that is
# worth a human looking at, not a code assumption.
_STALE_AFTER_QUARTERS = 2


def _release_age_warning(release_title: str, today: date | None = None) -> str | None:
    """None if the release looks current; otherwise a human-readable note
    on how far behind it is, for both the dataset's notes and a log line.

    Parses "DepMap YYQN Public" rather than trusting Figshare's own
    published_date, since a cached release note (read from disk, not
    re-fetched from Figshare) only ever has the title string to go on.
    """
    match = _RELEASE_TITLE.match(release_title)
    if not match:
        return None
    yy, q = int(match.group(1)), int(match.group(2))
    release_quarters = yy * 4 + (q - 1)  # e.g. 24Q4 -> 24*4+3 = 99

    today = today or date.today()
    current_quarters = today.year % 100 * 4 + (today.month - 1) // 3

    behind = current_quarters - release_quarters
    if behind <= _STALE_AFTER_QUARTERS:
        return None
    return (
        f"This is DepMap {yy:02d}Q{q} Public, {behind} quarters behind the current "
        "calendar quarter. DepMap stopped publishing full releases to Figshare "
        "after 24Q4; newer releases are portal-only and blocked by DepMap's "
        "own bot verification (Cloudflare Turnstile) from this app's "
        "programmatic loader -- see METHODS.md 7.5. Not a bug in this app; "
        "a known, monitored gap in the upstream distribution channel."
    )


def _resolve_release_urls(query: str = "DepMap Public", max_pages: int = 5) -> tuple[dict[str, str], str]:
    """Find the newest DepMap quarterly release on Figshare and return
    {filename: download_url} for the files this loader needs."""
    found: dict[int, dict] = {}
    for page in range(1, max_pages + 1):
        resp = httpx.post(
            FIGSHARE_SEARCH,
            json={"search_for": query, "page_size": 100, "page": page},
            timeout=30.0,
        )
        resp.raise_for_status()
        batch = resp.json()
        if not batch:
            break
        for a in batch:
            if _RELEASE_TITLE.match(a["title"]):
                found[a["id"]] = a
    if not found:
        raise RuntimeError("No 'DepMap <Q> Public' release found on Figshare.")
    latest = max(found.values(), key=lambda a: a["published_date"])

    detail = httpx.get(latest["url_public_api"], timeout=30.0)
    detail.raise_for_status()
    files = {f["name"]: f["download_url"] for f in detail.json()["files"]}
    missing = {EXPRESSION_FILENAME, MODEL_FILENAME} - files.keys()
    if missing:
        raise RuntimeError(f"Release {latest['title']} is missing expected files: {missing}")
    return files, latest["title"]


def load(lineage_filter: str | None = "Lymphoid", use_cache: bool = True) -> Dataset:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    matrix_path = CACHE_DIR / "matrix.parquet"
    model_path = CACHE_DIR / "model.parquet"
    release_note_path = CACHE_DIR / "release.txt"

    if use_cache and matrix_path.exists() and model_path.exists():
        matrix = pd.read_parquet(matrix_path)
        model = pd.read_parquet(model_path)
        release_title = release_note_path.read_text().strip() if release_note_path.exists() else "cached"
    else:
        urls, release_title = _resolve_release_urls()
        # rows=ModelID, cols="SYMBOL (ENTREZID)"; values are log2(TPM+1)
        expr = pd.read_csv(urls[EXPRESSION_FILENAME], index_col=0)
        model = pd.read_csv(urls[MODEL_FILENAME])
        matrix = expr.T  # -> rows=genes, cols=ModelID
        matrix.to_parquet(matrix_path)
        model.to_parquet(model_path)
        release_note_path.write_text(release_title)

    if lineage_filter and "OncotreeLineage" in model.columns:
        keep_ids = set(model.loc[model["OncotreeLineage"] == lineage_filter, "ModelID"])
        matrix = matrix[[c for c in matrix.columns if c in keep_ids]]

    samples = model[model["ModelID"].isin(matrix.columns)].copy()
    samples = samples.rename(columns={"ModelID": "sample_id"})
    samples["dataset_id"] = "depmap"
    samples["group_columns"] = samples.apply(
        lambda r: {
            "lineage": r.get("OncotreeLineage"),
            "subtype": r.get("OncotreeSubtype"),
            "cell_line_name": r.get("StrippedCellLineName", r.get("CellLineName")),
        },
        axis=1,
    )
    samples = samples[["sample_id", "dataset_id", "group_columns"]].reset_index(drop=True)
    matrix = matrix[samples["sample_id"].tolist()]

    features = pd.DataFrame({"feature_id": matrix.index})
    features["symbol"] = features["feature_id"].str.replace(r"\s*\(\d+\)$", "", regex=True)
    features["aliases"] = [[] for _ in range(len(features))]
    matrix.index = features["feature_id"]

    staleness = _release_age_warning(release_title)
    if staleness:
        log.warning("depmap dataset: %s", staleness)

    source = DatasetSource(
        dataset_id="depmap",
        display_name="DepMap cell line RNA-seq",
        accession=release_title,
        repository="DepMap (via Figshare)",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.LOG2_INTENSITY,  # log2(TPM+1), not raw TPM
        gene_identifier="symbol",
        n_samples=matrix.shape[1],
        notes=(
            (f"Filtered to lineage={lineage_filter!r}. " if lineage_filter else "All lineages. ")
            + "Values are log2(TPM+1), matching DepMap's native scale."
            + (f" {staleness}" if staleness else "")
        ),
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


from app.registry import DatasetDescriptor, register  # noqa: E402

register(
    DatasetDescriptor(
        dataset_id="depmap",
        display_name="DepMap cell line RNA-seq",
        loader=load,
        group_columns=("lineage", "subtype", "cell_line_name"),
        supports_survival=False,  # no clinical follow-up data — cell lines, not patients
    )
)
