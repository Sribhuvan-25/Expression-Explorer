"""
GSE66586 -- Neuroblastoma cell line gene expression profiling.

Gu et al. 2016 (PMID 26844271), City of Hope. 10 untreated cell lines on
Affymetrix Human Gene 1.0 ST (GPL6244, transcript-cluster level): 5
MYCN-amplified neuroblastoma lines, 3 MYCN non-amplified neuroblastoma
lines, and 2 non-malignant controls (bone-marrow MSCs and the HCN-1A
cortical neuronal line). Values are RMA log2 as submitted (Partek
Genomics Suite 6.6) -- not re-normalised here.

Labels come from the series matrix itself, not a curated table (contrast
gds4299.py): `!Sample_source_name_ch1` encodes the MYCN group
(`MYCN_Amplified_n` / `MYCN_NotAmplified_n` / `Control_n`) and
`!Sample_characteristics_ch1` the cell type. The parse is asserted to give
exactly the published 5 / 3 / 2 split, so a re-issued matrix that changes
its labelling fails loudly instead of silently regrouping samples.

**Cell-line identity.** Two pairs here are not independent lines:
  - GSM1625938 "BE(2)c" (ATCC CRL-2268) is SK-N-BE(2)-C (CVCL_0529), and
    GSM1625941 "SK-N-BE(2)c" is -- by its own ATCC catalogue number
    (CRL-2271) and CEL file name (`..._SK-N-BE2_...`) -- the *parental*
    SK-N-BE(2) (CVCL_0528), despite a sample title that reads like the
    subclone. Either way both come from one patient.
  - SH-SY5Y (CVCL_0019) is a subclone of SK-N-SH (CVCL_0531).
So the MYCN-amplified group is 5 samples from 4 patients, and the
non-amplified group 3 samples from 2. `cell_line_identity` carries the
Cellosaurus name resolved from each sample's ATCC catalogue number, so
this is visible rather than buried. No sample is dropped or merged.

Microarray log2 intensity: never pooled with RNA-seq datasets (METHODS 7.2).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from app.config import settings
from app.ingest._geo import (
    collapse_probes,
    download,
    probe_to_symbol,
    read_series_matrix_expression,
    read_series_matrix_header,
)
from app.models.contract import AssayType, Dataset, DatasetSource, ExpressionUnit

logger = logging.getLogger(__name__)

DATASET_ID = "gse66586"
CACHE_DIR = settings.cache_dir / DATASET_ID

SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66586/matrix/"
    "GSE66586_series_matrix.txt.gz"
)
PLATFORM_ANNOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
)

# Legacy symbols GPL6244's annotation (dated Aug 2016) still carries for
# genes since renamed by HGNC. Checked against the file for the MYC / PP2A
# panel: CIP2A is only present as KIAA1524 and PABIR1 only as FAM122A;
# neither current symbol has a probe of its own, so no collision arises on
# this platform today (collapse_probes would resolve one by highest mean if
# it did). PPME1 and MYCL are already annotated under their current names.
_SYMBOL_ALIASES = {
    "KIAA1524": "CIP2A",
    "FAM122A": "PABIR1",
}

MYCN_AMPLIFIED = "MYCN-amplified"
MYCN_NON_AMPLIFIED = "MYCN non-amplified"
NON_MALIGNANT = "non-malignant"

# !Sample_source_name_ch1 -> mycn_status
_SOURCE_NAME_PATTERNS = (
    (re.compile(r"^MYCN_Amplified_\d+$"), MYCN_AMPLIFIED),
    (re.compile(r"^MYCN_NotAmplified_\d+$"), MYCN_NON_AMPLIFIED),
    (re.compile(r"^Control_\d+$"), NON_MALIGNANT),
)

# The published design (Series_summary: "5 MYCN-amplified, 3 non-MYCN-
# amplified, and 2 non-malignant cell lines").
EXPECTED_GROUP_COUNTS = {MYCN_AMPLIFIED: 5, MYCN_NON_AMPLIFIED: 3, NON_MALIGNANT: 2}

# ATCC catalogue number (!Sample_biomaterial_provider_ch1) -> Cellosaurus
# recommended name. Resolved one by one against the Cellosaurus API
# (dr:"ATCC; <cat>") on 2026-11-15. Keyed on the catalogue number rather
# than the free-text sample title because the title is what disagrees:
# GSM1625941 is titled "SK-N-BE(2)c" but is CRL-2271, which Cellosaurus
# lists as the parental SK-N-BE(2), not the -C subclone (CRL-2268).
# PCS-500-012 is primary cells and has no Cellosaurus entry.
_ATCC_TO_CELLOSAURUS = {
    "PCS-500-012": "hBM-MSC (primary cells, ATCC PCS-500-012)",
    "CRL-10442": "HCN-1A",
    "CRL-2268": "SK-N-BE(2)-C",
    "CRL-2149": "SK-N-DZ",
    "CCL-127": "IMR-32",
    "CRL-2271": "SK-N-BE(2)",
    "CRL-2273": "CHP-212",
    "CRL-2137": "SK-N-AS",
    "CRL-2266": "SH-SY5Y",
    "HTB-11": "SK-N-SH",
}


def _single_row(header: dict[str, list[list[str]]], key: str) -> list[str]:
    rows = header.get(key)
    if not rows:
        raise ValueError(f"GSE66586 series matrix has no {key} row.")
    if len(rows) != 1:
        raise ValueError(f"GSE66586 series matrix has {len(rows)} {key} rows, expected 1.")
    return rows[0]


def _mycn_status(source_name: str) -> str:
    for pattern, label in _SOURCE_NAME_PATTERNS:
        if pattern.match(source_name):
            return label
    raise ValueError(f"Unrecognised GSE66586 source name {source_name!r}.")


def _cell_type(characteristics: list[list[str]], i: int) -> str:
    for row in characteristics:
        key, _, value = row[i].partition(":")
        if key.strip().lower() == "cell type":
            return value.strip()
    raise ValueError(f"Sample column {i} has no 'cell type:' characteristic.")


def _atcc_catalogue(provider: str) -> str:
    # GEO stores "ATCC:CRL<soft hyphen>10442" for HCN-1A; treat U+00AD as '-'.
    cat = provider.replace("\u00ad", "-").strip()
    if not cat.startswith("ATCC:"):
        raise ValueError(f"Unrecognised biomaterial provider {provider!r}.")
    return cat[len("ATCC:"):].strip()


def _parse_sample_labels(matrix_path: Path) -> pd.DataFrame:
    """One row per GSM: mycn_status, cell_line, cell_type, cell_line_identity.

    Fails loudly if the parse does not reproduce the published 5 / 3 / 2
    design, or if the cell-type characteristic contradicts the group (a
    sample labelled non-malignant must not be a neuroblastoma, and vice
    versa).
    """
    header = read_series_matrix_header(matrix_path)
    gsm = _single_row(header, "!Sample_geo_accession")
    titles = _single_row(header, "!Sample_title")
    sources = _single_row(header, "!Sample_source_name_ch1")
    providers = _single_row(header, "!Sample_biomaterial_provider_ch1")
    characteristics = header.get("!Sample_characteristics_ch1", [])

    n = len(gsm)
    for key, row in [("title", titles), ("source_name", sources), ("provider", providers)] + [
        ("characteristics", r) for r in characteristics
    ]:
        if len(row) != n:
            raise ValueError(f"GSE66586 {key} row has {len(row)} values for {n} samples.")

    records = []
    for i in range(n):
        status = _mycn_status(sources[i])
        cell_type = _cell_type(characteristics, i)
        is_nb = cell_type.lower() == "neuroblastoma"
        if is_nb == (status == NON_MALIGNANT):
            raise ValueError(
                f"{gsm[i]}: source name {sources[i]!r} ({status}) contradicts "
                f"cell type {cell_type!r}."
            )
        cat = _atcc_catalogue(providers[i])
        if cat not in _ATCC_TO_CELLOSAURUS:
            raise ValueError(f"{gsm[i]}: ATCC {cat!r} has no Cellosaurus mapping.")
        records.append(
            {
                "sample_id": gsm[i],
                "mycn_status": status,
                "cell_line": titles[i],
                "cell_type": cell_type,
                "cell_line_identity": _ATCC_TO_CELLOSAURUS[cat],
            }
        )
    labels = pd.DataFrame.from_records(records)

    counts = labels["mycn_status"].value_counts().to_dict()
    if counts != EXPECTED_GROUP_COUNTS:
        raise ValueError(
            f"GSE66586 labels parsed to {counts}; expected {EXPECTED_GROUP_COUNTS} "
            "(5 MYCN-amplified / 3 non-amplified / 2 non-malignant). Refusing to "
            "load a regrouped cohort."
        )
    if labels["sample_id"].duplicated().any():
        raise ValueError("GSE66586 series matrix lists a GSM accession twice.")
    return labels


def _symbol_table(annot_path: Path) -> pd.DataFrame:
    """probe_id, legacy (annotation) symbol, current symbol."""
    raw = probe_to_symbol(annot_path, {})
    return pd.DataFrame(
        {"probe_id": raw.index, "annot_symbol": raw.to_numpy(), "symbol": raw.replace(_SYMBOL_ALIASES).to_numpy()}
    )


def _log_alias_collisions(table: pd.DataFrame, expr: pd.DataFrame) -> None:
    """Make any legacy/current collision visible. Resolution itself is
    `collapse_probes`' highest-mean rule, applied to all of a gene's
    probes regardless of which name they were annotated under."""
    present = table[table["probe_id"].isin(expr.index)]
    renamed = set(present.loc[present["annot_symbol"] != present["symbol"], "symbol"])
    for sym in sorted(renamed):
        probes = present[present["symbol"] == sym]
        if probes["annot_symbol"].nunique() > 1:
            means = expr.loc[probes["probe_id"]].mean(axis=1)
            logger.info(
                "GSE66586: %s has probes under both legacy and current symbols %s; "
                "keeping highest-mean probe %s (mean %.3f).",
                sym, sorted(probes["annot_symbol"].unique()), means.idxmax(), means.max(),
            )


def load(use_cache: bool = True) -> Dataset:
    matrix_path = download(SERIES_MATRIX_URL, CACHE_DIR / "GSE66586_series_matrix.txt.gz")
    annot_path = download(PLATFORM_ANNOT_URL, CACHE_DIR / "GPL6244.annot.gz")

    labels = _parse_sample_labels(matrix_path)
    expr = read_series_matrix_expression(matrix_path)
    if list(expr.columns) != list(labels["sample_id"]):
        raise ValueError(
            f"GSE66586 expression columns {list(expr.columns)} do not match the "
            f"parsed sample accessions {list(labels['sample_id'])}."
        )

    table = _symbol_table(annot_path)
    _log_alias_collisions(table, expr)
    matrix = collapse_probes(expr, pd.Series(table["symbol"].to_numpy(), index=table["probe_id"].to_numpy()))

    # Legacy names actually seen on this platform, recorded on the current
    # symbol's feature so a reader can see why e.g. CIP2A is present.
    legacy = (
        table[(table["annot_symbol"] != table["symbol"]) & table["probe_id"].isin(expr.index)]
        .groupby("symbol")["annot_symbol"]
        .agg(lambda s: sorted(set(s)))
    )

    labels = labels.set_index("sample_id")
    group_keys = ["mycn_status", "cell_line", "cell_type", "cell_line_identity"]
    samples = pd.DataFrame(
        {
            "sample_id": list(matrix.columns),
            "dataset_id": DATASET_ID,
            "group_columns": [{k: labels.at[s, k] for k in group_keys} for s in matrix.columns],
        }
    )

    features = pd.DataFrame(
        {
            "feature_id": matrix.index,
            "symbol": matrix.index,
            "aliases": [list(legacy.get(sym, [])) for sym in matrix.index],
        }
    )

    source = DatasetSource(
        dataset_id=DATASET_ID,
        display_name="GSE66586 (neuroblastoma cell lines, MYCN status)",
        accession="GSE66586",
        repository="GEO",
        assay_type=AssayType.MICROARRAY,
        expression_unit=ExpressionUnit.LOG2_INTENSITY,
        gene_identifier="symbol",
        n_samples=matrix.shape[1],
        disease_area="Neuroblastoma",
        notes=(
            "Gu et al. 2016 (PMID 26844271). Affymetrix Human Gene 1.0 ST (GPL6244), "
            "RMA log2 as submitted (Partek). 5 MYCN-amplified, 3 MYCN non-amplified "
            "and 2 non-malignant cell lines, labels parsed from the series matrix. "
            "Multi-probe genes collapsed to the highest-mean probe, independent of "
            "MYCN status (METHODS.md 7.1). Legacy symbols renamed to current HGNC "
            "(KIAA1524 -> CIP2A, FAM122A -> PABIR1). Not all lines are independent: "
            "BE(2)c and SK-N-BE(2)c are SK-N-BE(2)-C and its parent SK-N-BE(2) (one "
            "patient), and SH-SY5Y is a subclone of SK-N-SH -- see the "
            "cell_line_identity column. Small groups limit what any test can show: "
            "the smallest possible two-sided Mann-Whitney p is 0.036 for amplified vs "
            "non-amplified (5 vs 3), 0.095 for amplified vs non-malignant (5 vs 2) and "
            "0.2 for non-amplified vs non-malignant (3 vs 2). Microarray intensities: "
            "do not pool with RNA-seq datasets."
        ),
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


from app.registry import DatasetDescriptor, register  # noqa: E402

register(
    DatasetDescriptor(
        dataset_id=DATASET_ID,
        display_name="GSE66586 (neuroblastoma cell lines, MYCN status)",
        loader=load,
        group_columns=("mycn_status", "cell_line", "cell_type", "cell_line_identity"),
        supports_survival=False,  # cell lines, no clinical follow-up
    )
)
