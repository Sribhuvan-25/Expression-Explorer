"""
DepMap CRISPR gene-effect scores -- an auxiliary layer on the existing
`depmap` dataset, not a dataset of its own.

Same cell lines, same ModelID join key, same Figshare release bundle as
the expression matrix `depmap.py` already downloads: this is those lines
measured a second way (knock out gene X, how much does the line care?),
which is exactly what `AuxLayer` exists for -- see METHODS.md 8.1 for why
this isn't registered as a separate dataset.

Values are Chronos gene-effect scores. The scale is anchored so that 0 is
"no effect of knocking this gene out" and -1 is roughly the median of
known common-essential genes, so *more negative = stronger dependency*.
Positive values (a line growing better without the gene) exist but are
rarer and usually small.

Inherits `depmap.py`'s release-staleness limitation wholesale: it resolves
the same Figshare release, so it can only ever reach 24Q4 for the same
reason (METHODS.md 7.5). That's deliberate -- pulling CRISPR from a
*different* release than the expression matrix it's joined against would
silently mismatch the two layers' cell-line rosters.
"""
from __future__ import annotations

import pandas as pd

from app.config import settings
from app.ingest.depmap import _resolve_release_urls
from app.models.contract import AuxLayer, AuxMatrix

CACHE_DIR = settings.cache_dir / "depmap"
CRISPR_FILENAME = "CRISPRGeneEffect.csv"


def load_crispr(use_cache: bool = True) -> AuxMatrix | None:
    """Gene-effect matrix (gene symbol x ModelID), or None if this release
    doesn't carry the file.

    Returns None rather than raising on a missing file: the CRISPR layer is
    additive, and a release without it should degrade the `depmap` dataset
    to expression-only rather than taking the whole dataset down with it.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / "crispr_gene_effect.parquet"

    if use_cache and cache_path.exists():
        matrix = pd.read_parquet(cache_path)
    else:
        urls, _release_title = _resolve_release_urls()
        if CRISPR_FILENAME not in urls:
            return None
        # rows=ModelID, cols="SYMBOL (ENTREZID)" -- same orientation and
        # same column-naming convention as the expression file, so the
        # transpose-and-strip handling mirrors depmap.py exactly.
        raw = pd.read_csv(urls[CRISPR_FILENAME], index_col=0)
        matrix = raw.T
        matrix.index = matrix.index.str.replace(r"\s*\(\d+\)$", "", regex=True)
        # A handful of symbols appear twice on this platform; keep the
        # first rather than silently averaging two different guides'
        # scores into a number that matches neither.
        matrix = matrix[~matrix.index.duplicated(keep="first")]
        matrix.to_parquet(cache_path)

    return AuxMatrix(
        layer=AuxLayer.CRISPR_GENE_EFFECT,
        matrix=matrix,
        value_label="Gene effect (Chronos)",
        value_description=(
            "CRISPR knockout gene-effect score. 0 = no effect; -1 is about the median "
            "of known common-essential genes, so more negative means the line depends "
            "more strongly on that gene."
        ),
        source_note=f"DepMap {CRISPR_FILENAME}, same Figshare release as the expression matrix.",
    )
