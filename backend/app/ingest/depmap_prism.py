"""
DepMap PRISM Repurposing drug-sensitivity screen -- an auxiliary layer on
the existing `depmap` dataset.

Same cell lines (ACH- ModelIDs) as the expression matrix, measured a third
way: expose the line to a compound, how much does it die? Attached as an
`AuxLayer` rather than registered as its own dataset for the same reason
as the CRISPR layer -- see METHODS.md 8.1.

Two things make this loader materially more involved than the CRISPR one,
both verified against the real files rather than assumed:

1. **Separate Figshare release.** PRISM ships as "Repurposing Public
   YYQN", NOT inside the "DepMap YYQN Public" bundle, so `depmap.py`'s
   release-title regex will never match it and this module has to resolve
   its own release.
2. **Already in (feature x sample) orientation, unlike DepMap's other
   files.** The data matrix ships as compounds-as-ROWS x cell-lines-as-
   COLUMNS, which happens to match this app's (feature x sample)
   convention directly -- so it is loaded as-is with NO transpose, while
   `depmap_crispr.py` and `depmap.py` both DO transpose their sources.
   Worth stating explicitly because the inconsistency between sibling
   loaders is surprising and invites a "fix" that would silently break
   this one.

The matrix index is a Broad compound id (`BRD:BRD-K00104122-001-01-9`),
which is meaningless to a reader on its own, so the compound list is
joined in as this layer's `features` frame to carry drug name and
mechanism of action.
"""
from __future__ import annotations

import re

import httpx
import pandas as pd

from app.config import settings
from app.models.contract import AuxLayer, AuxMatrix

FIGSHARE_SEARCH = "https://api.figshare.com/v2/articles/search"
CACHE_DIR = settings.cache_dir / "depmap_prism"

# Distinct from depmap.py's "^DepMap \d\dQ\d Public$" -- see module docstring.
_RELEASE_TITLE = re.compile(r"^Repurposing Public (\d\d)Q(\d)$")

_MATRIX_SUFFIX = "_Extended_Primary_Data_Matrix.csv"
_COMPOUND_SUFFIX = "_Extended_Primary_Compound_List.csv"


def _resolve_release() -> tuple[dict[str, str], str] | None:
    """Newest "Repurposing Public YYQN" release on Figshare, as
    {filename: download_url}. None if no matching release is published."""
    found: dict[int, dict] = {}
    for page in range(1, 4):
        resp = httpx.post(
            FIGSHARE_SEARCH,
            json={"search_for": "Repurposing Public", "page_size": 100, "page": page},
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
        return None
    latest = max(found.values(), key=lambda a: a["published_date"])

    detail = httpx.get(latest["url_public_api"], timeout=30.0)
    detail.raise_for_status()
    files = {f["name"]: f["download_url"] for f in detail.json()["files"]}
    return files, latest["title"]


def _find(files: dict[str, str], suffix: str) -> str | None:
    for name, url in files.items():
        if name.endswith(suffix):
            return url
    return None


def load_drug_sensitivity(use_cache: bool = True) -> AuxMatrix | None:
    """Drug-sensitivity matrix (compound x ModelID), or None if no
    Repurposing release is reachable. Additive layer -- never fatal to the
    dataset it attaches to."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    matrix_path = CACHE_DIR / "drug_sensitivity.parquet"
    compounds_path = CACHE_DIR / "compounds.parquet"
    release_path = CACHE_DIR / "release.txt"

    if use_cache and matrix_path.exists() and compounds_path.exists():
        matrix = pd.read_parquet(matrix_path)
        compounds = pd.read_parquet(compounds_path)
        release_title = release_path.read_text().strip() if release_path.exists() else "cached"
    else:
        resolved = _resolve_release()
        if resolved is None:
            return None
        files, release_title = resolved
        matrix_url = _find(files, _MATRIX_SUFFIX)
        compound_url = _find(files, _COMPOUND_SUFFIX)
        if matrix_url is None or compound_url is None:
            return None

        # rows=BRD compound id, cols=ACH ModelID. That is ALREADY this
        # app's (feature x sample) convention, so -- unlike depmap.py and
        # depmap_crispr.py, whose sources are sample-rows and need `.T` --
        # this one is loaded as-is. Do not "fix" the missing transpose.
        matrix = pd.read_csv(matrix_url, index_col=0, low_memory=False)
        matrix.index = matrix.index.astype(str)

        compound_raw = pd.read_csv(compound_url, low_memory=False)
        # One row per (compound, dose, screen); collapse to one row per
        # compound id so this is a usable feature table. Keeping the first
        # occurrence is enough for the name/MOA fields, which don't vary
        # within a compound id.
        compound_raw = compound_raw.drop_duplicates(subset=["IDs"], keep="first")
        compounds = pd.DataFrame(
            {
                "feature_id": compound_raw["IDs"].astype(str),
                "symbol": compound_raw["Drug.Name"].fillna(compound_raw["IDs"]).astype(str),
                "moa": compound_raw["MOA"].fillna("").astype(str),
                "target": compound_raw["repurposing_target"].fillna("").astype(str),
            }
        )

        matrix.to_parquet(matrix_path)
        compounds.to_parquet(compounds_path)
        release_path.write_text(release_title)

    return AuxMatrix(
        layer=AuxLayer.DRUG_SENSITIVITY,
        matrix=matrix,
        value_label="Viability log2 fold-change",
        value_description=(
            "PRISM Repurposing primary screen. Log2 fold-change in cell viability after "
            "compound exposure, relative to untreated control -- more negative means the "
            "compound killed more of that cell line."
        ),
        source_note=f"DepMap {release_title} (Figshare), Extended Primary Data Matrix.",
        features=compounds,
    )
