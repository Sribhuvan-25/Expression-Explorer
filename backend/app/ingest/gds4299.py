"""
GSE28703 / GDS4299 — ETP vs non-ETP T-ALL expression profiling.

Zhang et al., Nature 2012 (481:157-163; PMID 22237106), the landmark
ETP-ALL genomic study. 52 diagnosis tumour samples on Affymetrix HT
HG-U133+ PM (GPL13158), already log2-transformed by the submitter.

Why this dataset: TARGET-ALL-P2 classifies ETP status for only 190 of
its 469 samples, and the ETP arm itself is just 19 patients -- the
smallest and most important group in the comparison the project cares
about. Adding this cohort's 12 ETP samples takes the pooled ETP count
from 19 to 31 (+63%).

Note on combining: this is *microarray log2 intensity*, while TARGET is
*RNA-seq TPM*. Those are not on a comparable scale and must not be
pooled naively -- see METHODS.md section 7.2. This module's job is only
to expose the cohort as its own dataset; how (or whether) it gets
combined with others is a separate decision made above this layer.

Data sources, both public:
  - Expression: GEO series matrix for GSE28703
  - Platform annotation: GEO's official GPL13158 annotation file
  - ETP/non-ETP labels: curated table (see below) -- GEO's own
    !Sample_characteristics carries no ETP status whatsoever, so an
    external label source is mandatory, exactly as for TARGET.
"""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import httpx
import pandas as pd

from app.config import settings
from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

CACHE_DIR = settings.cache_dir / "gds4299"

SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28703/matrix/"
    "GSE28703_series_matrix.txt.gz"
)
PLATFORM_ANNOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz"
)

# GEO carries no ETP annotation for this series -- every sample's
# !Sample_characteristics_ch1 reads only "cell type: tumor cells" -- so
# grouping is impossible without an external label source, exactly as for
# TARGET (see liu2017_etp.py). These labels are therefore scientific
# input, not a cache: they are vendored into the repo rather than
# downloaded or read from a sibling checkout, so the dataset loads
# identically in development and in a deployed container.
#
# Verified before use: this table splits 12 ETP / 40 non-ETP, matching the
# published composition of GSE28703 (Zhang et al., Nature 2012), and joins
# to all 52 samples in the matrix with none missing or extra.
SAMPLE_LABELS_PATH = Path(__file__).parent / "data" / "gds4299_etp_labels.csv"

# Symbols that changed since GPL13158 was annotated (2016). Without these
# a user searching the current symbol gets "gene not found" even though
# the probe is present under the legacy name.
_SYMBOL_ALIASES = {
    "KIAA1524": "CIP2A",
    "PME1": "PPME1",
}

# Normalize to the hyphenated labels used everywhere else in this project
# (TARGET's Liu-2017 enrichment emits 'ETP' / 'near-ETP' / 'non-ETP').
_SUBTYPE_MAP = {"ETP": "ETP", "nonETP": "non-ETP", "non-ETP": "non-ETP"}

# GEO uses these as "no mapping", not as gene symbols.
_NULL_SYMBOLS = {"---", "NA", ""}


def _download(url: str, dest: Path) -> Path:
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    with httpx.stream("GET", url, follow_redirects=True, timeout=120.0) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as fh:
            for chunk in r.iter_bytes():
                fh.write(chunk)
        tmp.rename(dest)
    return dest


def _read_geo_table(path: Path, begin: str, end: str) -> pd.DataFrame:
    """Pull the delimited table out of a GEO file's !..._begin/!..._end block."""
    with gzip.open(path, "rt", errors="replace") as fh:
        text = fh.read()
    start, stop = text.index(begin), text.index(end)
    body = text[start:stop].split("\n", 1)[1]
    return pd.read_csv(io.StringIO(body), sep="\t", low_memory=False)


def _probe_to_symbol(annot_path: Path) -> pd.Series:
    """probe_id -> gene symbol, for probes that actually map to a gene.

    A probe may list several symbols ('A /// B'); GEO also uses '---' and
    'NA' to mean "no mapping". Both were handled in the prior verified
    analysis of this dataset and are handled the same way here: take the
    first symbol, drop the sentinels. ~44k of 54,715 probes survive.
    """
    ann = _read_geo_table(annot_path, "!platform_table_begin", "!platform_table_end")
    symbol_col = "Gene symbol" if "Gene symbol" in ann.columns else "Gene Symbol"
    ann = ann[["ID", symbol_col]].dropna()
    ann["ID"] = ann["ID"].astype(str).str.strip()
    symbols = (
        ann[symbol_col]
        .astype(str)
        .str.split(r"[;,/\s]+", regex=True)
        .str[0]
        .str.strip()
    )
    symbols = symbols.replace(_SYMBOL_ALIASES)
    keep = ~symbols.isin(_NULL_SYMBOLS)
    return pd.Series(symbols[keep].to_numpy(), index=ann["ID"][keep].to_numpy())


def _collapse_probes(expr: pd.DataFrame, probe_symbol: pd.Series) -> pd.DataFrame:
    """One row per gene, choosing the probe with the highest mean expression.

    The selection rule is deliberately independent of any sample grouping.
    The earlier exploratory analysis picked, per gene, the probe with the
    lowest ETP-vs-non-ETP p-value -- that is circular for inference
    (choose the probe that best separates the groups, then test whether
    the groups separate), so it is not used here. See METHODS.md 7.1.

    Highest mean also does the obvious right thing on this platform: MYCN
    has 10 probes, one at mean 5.36 and nine sitting at ~2.7 background.
    """
    shared = expr.index.intersection(probe_symbol.index)
    expr = expr.loc[shared]
    genes = probe_symbol.loc[shared]

    order = expr.mean(axis=1).sort_values(ascending=False)
    # First occurrence per gene in mean-descending order == highest mean.
    best_probe = genes.loc[order.index].reset_index().drop_duplicates(subset=0, keep="first")
    best_probe.columns = ["probe_id", "symbol"]

    collapsed = expr.loc[best_probe["probe_id"]]
    collapsed.index = pd.Index(best_probe["symbol"], name="feature_id")
    return collapsed.sort_index()


def _load_labels() -> pd.Series:
    if not SAMPLE_LABELS_PATH.exists():
        raise FileNotFoundError(
            f"ETP/non-ETP labels not found at {SAMPLE_LABELS_PATH}. GEO carries no ETP "
            "annotation for GSE28703, so this dataset cannot be grouped without them."
        )
    labels = pd.read_csv(SAMPLE_LABELS_PATH)
    missing = {"sample_id", "subtype"} - set(labels.columns)
    if missing:
        raise ValueError(f"Sample label file is missing columns: {sorted(missing)}")
    mapped = labels["subtype"].map(_SUBTYPE_MAP)
    unknown = labels.loc[mapped.isna(), "subtype"].unique()
    if len(unknown):
        raise ValueError(f"Unrecognized subtype value(s) in labels: {list(unknown)}")
    return pd.Series(mapped.to_numpy(), index=labels["sample_id"].to_numpy())


def load(use_cache: bool = True) -> Dataset:
    matrix_path = _download(SERIES_MATRIX_URL, CACHE_DIR / "GSE28703_series_matrix.txt.gz")
    annot_path = _download(PLATFORM_ANNOT_URL, CACHE_DIR / "GPL13158.annot.gz")

    expr = _read_geo_table(matrix_path, "!series_matrix_table_begin", "!series_matrix_table_end")
    expr = expr.set_index(expr.columns[0])
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = [str(c).strip('"') for c in expr.columns]

    matrix = _collapse_probes(expr, _probe_to_symbol(annot_path))

    labels = _load_labels()
    # Only keep samples we can actually group -- an unlabelled sample can't
    # take part in the ETP comparison this dataset exists for.
    usable = [s for s in matrix.columns if s in labels.index]
    matrix = matrix[usable]

    samples = pd.DataFrame(
        {
            "sample_id": usable,
            "dataset_id": "gds4299",
            "group_columns": [{"etp_status": labels[s]} for s in usable],
        }
    )

    features = pd.DataFrame(
        {
            "feature_id": matrix.index,
            "symbol": matrix.index,
            "aliases": [[] for _ in range(len(matrix))],
        }
    )

    source = DatasetSource(
        dataset_id="gds4299",
        display_name="GSE28703 / GDS4299 (ETP vs non-ETP T-ALL)",
        accession="GSE28703 (GDS4299)",
        repository="GEO",
        assay_type=AssayType.MICROARRAY,
        expression_unit=ExpressionUnit.LOG2_INTENSITY,
        gene_identifier="symbol",
        n_samples=matrix.shape[1],
        disease_area="ETP-ALL",
        notes=(
            "Zhang et al., Nature 2012 (PMID 22237106). Affymetrix HT HG-U133+ PM "
            "(GPL13158), log2 intensity as submitted. Multi-probe genes collapsed to "
            "the probe with the highest mean expression -- a rule independent of ETP "
            "status, see METHODS.md 7.1. ETP/non-ETP labels are curated (GEO carries "
            "no ETP annotation for this series). Microarray intensities are not on the "
            "same scale as RNA-seq TPM; do not pool with TARGET without normalisation."
        ),
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


from app.registry import DatasetDescriptor, register  # noqa: E402

register(
    DatasetDescriptor(
        dataset_id="gds4299",
        display_name="GSE28703 / GDS4299 (ETP vs non-ETP T-ALL)",
        loader=load,
        group_columns=("etp_status",),
        supports_survival=False,  # expression profiling only, no clinical follow-up
    )
)
