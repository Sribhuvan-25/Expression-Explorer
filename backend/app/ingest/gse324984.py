"""
GSE324984 -- MYCL (L-MYC) gain and loss in prostate cancer cell lines.

"A MYC Family Switch: L-MYC Drives and Maintains Neuroendocrine Lineage
Programs in Prostate Cancer" (PMID 42001796). Bulk RNA-seq, hg38, 8
samples forming TWO independent experiments:

    C4-2B (adenocarcinoma)   empty vector x2   vs  MYCL overexpression x2
    NCI-H660 (NEPC)          shCtrl x2         vs  shMYCL knockdown x2

The design is the point of this module. The two cell lines differ far
more from each other than any perturbation changes them (replicate r =
0.99, cross-line r = 0.89 on log2 TPM), so the only valid perturbation
comparisons are *within* a cell line. They are declared here as
pre-specified contrasts (METHODS.md 11.1) and computed at ingest with the
count-based pipeline in app/analysis/contrast.py -- n=2 per group is too
small for the app's Mann-Whitney default to ever reach p < 0.05
(METHODS.md 11.2).

Two upstream files, both from the submitter:
  - TPM matrix: shown as the dataset's expression values.
  - gene counts: used for the statistics.

**The TPM file is truncated** (METHODS.md 11.3). It has exactly 65,535
gene rows -- the legacy Excel row limit -- and is a strict prefix of the
78,691-row counts file. TPM values for the genes it does contain are as
submitted and unaffected; the missing 13,156 genes are the tail of the
Ensembl id ordering (ENSG00000297233 onward, almost all recently-added
lncRNA/novel loci). So expression *display* covers 65,535 genes and the
*statistics* cover all 78,691, and both numbers are stated in the notes.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

from app.analysis.contrast import METHOD_LABEL, run_contrast
from app.config import settings
from app.models.contract import (
    AssayType,
    ContrastResult,
    ContrastSpec,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)

DATASET_ID = "gse324984"
DISPLAY_NAME = "GSE324984 (MYCL gain/loss, C4-2B & NCI-H660 prostate cancer)"
CACHE_DIR = settings.cache_dir / DATASET_ID

_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324984"
TPM_URL = f"{_BASE}/suppl/GSE324984_RNAseq_TPM_matrix.txt.gz"
COUNTS_URL = f"{_BASE}/suppl/GSE324984_RNAseq_gene_counts.txt.gz"
SERIES_MATRIX_URL = f"{_BASE}/matrix/GSE324984_series_matrix.txt.gz"

# GEO "treatment:" characteristic -> arm label shown in the UI, per cell line.
# Every value the series uses must be listed; an unknown one fails ingest
# rather than producing an ungroupable sample.
_ARMS = {
    ("C4-2B", "empty vector"): "C4-2B control (empty vector)",
    ("C4-2B", "MYCL overexpression"): "C4-2B MYCL overexpression",
    ("NCI-H660", "shCtrl"): "NCI-H660 control (shCtrl)",
    ("NCI-H660", "shMYCL"): "NCI-H660 MYCL knockdown (shMYCL)",
}

CONTRASTS = (
    ContrastSpec(
        contrast_id="c42b_mycl_oe",
        label="C4-2B: MYCL overexpression vs empty vector",
        group_column="arm",
        reference="C4-2B control (empty vector)",
        treatment="C4-2B MYCL overexpression",
        stratum_column="cell_line",
        stratum_value="C4-2B",
        method=METHOD_LABEL,
        direction_note="log2FC > 0: higher with MYCL overexpressed",
    ),
    ContrastSpec(
        contrast_id="h660_mycl_kd",
        label="NCI-H660: MYCL knockdown vs shCtrl",
        group_column="arm",
        reference="NCI-H660 control (shCtrl)",
        treatment="NCI-H660 MYCL knockdown (shMYCL)",
        stratum_column="cell_line",
        stratum_value="NCI-H660",
        method=METHOD_LABEL,
        direction_note="log2FC > 0: higher with MYCL knocked down",
    ),
)

# Current HGNC symbol -> earlier symbols a researcher may still use. The
# submitter's annotation already uses the current names; these are
# recorded as aliases so the gene can be found, and named, under either.
_LEGACY_SYMBOLS = {
    "CIP2A": ["KIAA1524"],
    "PABIR1": ["FAM122A"],
    "PPME1": ["PME1"],
}


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


def _sample_sheet(series_matrix: Path) -> pd.DataFrame:
    """One row per GSM from the series matrix header: accession, title,
    cell line, treatment, library name (which is the column header used
    in the supplementary matrices)."""
    rows: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(series_matrix, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            values = [v.strip().strip('"') for v in values]
            if key == "!Sample_characteristics_ch1":
                characteristics.append(values)
            else:
                rows.setdefault(key, values)

    sheet = pd.DataFrame(
        {
            "sample_id": rows["!Sample_geo_accession"],
            "title": rows["!Sample_title"],
            "library_name": [re.sub(r"^Library name:\s*", "", v) for v in rows["!Sample_description"]],
        }
    )
    for values in characteristics:
        key = values[0].split(":", 1)[0].strip()
        sheet[key] = [v.split(":", 1)[1].strip() if ":" in v else None for v in values]
    missing = {"cell line", "treatment"} - set(sheet.columns)
    if missing:
        raise ValueError(f"Series matrix lacks characteristics {sorted(missing)}; cannot assign arms.")

    sheet["arm"] = [_ARMS.get((c, t)) for c, t in zip(sheet["cell line"], sheet["treatment"])]
    if sheet["arm"].isna().any():
        bad = sheet.loc[sheet["arm"].isna(), ["sample_id", "cell line", "treatment"]]
        raise ValueError(f"Unrecognised cell line/treatment combination(s):\n{bad}")
    counts = sheet["arm"].value_counts()
    if len(sheet) != 8 or set(counts.index) != set(_ARMS.values()) or not (counts == 2).all():
        raise ValueError(f"Expected 2 samples in each of 4 arms, got {counts.to_dict()}")
    if sheet["library_name"].duplicated().any():
        raise ValueError("Duplicate library names in series matrix; cannot map matrix columns.")
    return sheet


def _read_supplement(path: Path, sheet: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """(values indexed by unversioned Ensembl id with GSM columns, symbol by id)."""
    frame = pd.read_csv(path, sep="\t")
    if list(frame.columns[:2]) != ["gene_id", "gene_name"]:
        raise ValueError(f"{path.name}: unexpected leading columns {list(frame.columns[:2])}")
    libraries = list(frame.columns[2:])
    if sorted(libraries) != sorted(sheet["library_name"]):
        raise ValueError(
            f"{path.name}: columns {libraries} do not match the series' library names "
            f"{list(sheet['library_name'])}"
        )
    ids = frame["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    if ids.duplicated().any():
        raise ValueError(f"{path.name}: Ensembl ids collide once versions are stripped.")
    by_library = dict(zip(sheet["library_name"], sheet["sample_id"]))
    values = frame[libraries].rename(columns=by_library)[list(sheet["sample_id"])]
    values.index = pd.Index(ids, name="feature_id")
    if values.isna().any().any() or (values < 0).any().any():
        raise ValueError(f"{path.name}: contains missing or negative values.")
    return values.astype(float), pd.Series(frame["gene_name"].astype(str).to_numpy(), index=values.index)


def _check_perturbation(values: pd.DataFrame, symbols: pd.Series, sheet: pd.DataFrame, what: str) -> dict:
    """Guard against mislabelled columns: MYCL itself must move the way
    each arm says it should. Run on BOTH files -- the TPM drives what is
    shown, the counts drive the statistics, and a swap in either alone
    would otherwise flow silently into every result. Returns the per-arm
    MYCL means so the check is visible in the notes. For counts, `values`
    should be library-size scaled (CPM) so unequal depths cannot mask or
    fake a change."""
    mycl = symbols.index[symbols == "MYCL"]
    if len(mycl) != 1:
        raise ValueError(f"Expected exactly one MYCL row, found {len(mycl)}.")
    row = values.loc[mycl[0]]
    means = {arm: float(row[sheet.loc[sheet["arm"] == arm, "sample_id"]].mean()) for arm in _ARMS.values()}
    for spec in CONTRASTS:
        ref, trt = means[spec.reference], means[spec.treatment]
        expect_up = "overexpression" in spec.treatment
        if (trt > ref) != expect_up:
            raise ValueError(
                f"MYCL {what} in '{spec.treatment}' ({trt:.1f}) vs '{spec.reference}' ({ref:.1f}) "
                "contradicts the arm labels -- columns may be mislabelled."
            )
    return means


def _tpm_direction_disagreement(tpm: pd.DataFrame, contrast: ContrastResult) -> tuple[int, int]:
    """(genes whose mean-TPM change points against log2FC, genes judged).
    Judged = tested, present in the TPM matrix, arm means differ (same rule
    as the panel's `display_direction_agrees`)."""
    table = contrast.table[contrast.table["tested"] & contrast.table["in_expression_matrix"]]
    ref = tpm.loc[table.index, contrast.reference_samples].mean(axis=1)
    trt = tpm.loc[table.index, contrast.treatment_samples].mean(axis=1)
    judged = (ref != trt) & (table["log2_fc"] != 0)
    against = (trt > ref) != (table["log2_fc"] > 0)
    return int((against & judged).sum()), int(judged.sum())


def _gene_identity(counts: pd.DataFrame, tpm: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    """Per Ensembl id (every row of the complete counts file): the symbol
    to show, whether this row is the one a symbol search should reach, and
    how many loci share its name.

    1,595 TPM rows (more in the counts) share a name with another id --
    Y_RNA/snoRNA families, and PAR X/Y copies such as GTPBP6, CD99, PLCXD1,
    whose reads the aligner splits arbitrarily between the two copies
    (opposite-signed fold-changes on the X and Y copy were seen in review).
    For each name the most-expressed row keeps it (scored as below); the
    others show their Ensembl id. The choice:
      - is independent of the sample groups (METHODS.md 7.1 reasoning);
      - is made over the *complete* counts, not the truncated TPM file --
        the expressed CAST (ENSG00000310517) exists only in the counts,
        and choosing within the TPM rows sent "CAST" to a near-silent
        locus (found in review);
      - is the same row in every view (Compare, panel, top-genes list).
    """
    # Score: mean TPM when every locus with that name has a TPM value (TPM
    # is length-normalised; raw counts favour longer loci -- C4orf36's two
    # loci have equal mean counts, 24.7 vs 25.0, but TPM 1.77 vs 0.71); mean
    # count when any of them lies past the TPM truncation, since TPM cannot
    # rank a locus it lacks (CAST). Computed per name, so the two scales are
    # never compared with each other.
    in_tpm = pd.Series(counts.index.isin(tpm.index), index=counts.index)
    all_in_tpm = in_tpm.groupby(symbols).transform("all")
    score = counts.mean(axis=1).where(~all_in_tpm, tpm.mean(axis=1).reindex(counts.index))
    order = score.sort_values(ascending=False, kind="mergesort").index
    primary = (~symbols.loc[order].duplicated(keep="first")).reindex(counts.index)
    n_loci = symbols.map(symbols.value_counts())
    return pd.DataFrame(
        {
            "gene_name": symbols,
            "symbol": np.where(primary, symbols, counts.index),
            "primary_for_name": primary.astype(bool),
            "n_loci_with_name": n_loci.astype(int),
        },
        index=counts.index,
    )


def _features(tpm: pd.DataFrame, identity: pd.DataFrame) -> pd.DataFrame:
    """One row per matrix (TPM) gene. A row whose name belongs to another
    locus shows its Ensembl id as symbol and carries NO alias for the
    shared name -- `aliases` means "earlier names of this gene", and
    putting shared names there made the legacy-name fallback resolve
    e.g. "C4orf36" to a secondary locus (found in review)."""
    ident = identity.loc[tpm.index]
    aliases = [
        list(_LEGACY_SYMBOLS.get(name, [])) if is_primary else []
        for name, is_primary in zip(ident["gene_name"], ident["primary_for_name"])
    ]
    return pd.DataFrame({"feature_id": tpm.index, "symbol": ident["symbol"].to_numpy(), "aliases": aliases})


def load(use_cache: bool = True) -> Dataset:
    series_path = _download(SERIES_MATRIX_URL, CACHE_DIR / "GSE324984_series_matrix.txt.gz")
    tpm_path = _download(TPM_URL, CACHE_DIR / "GSE324984_RNAseq_TPM_matrix.txt.gz")
    counts_path = _download(COUNTS_URL, CACHE_DIR / "GSE324984_RNAseq_gene_counts.txt.gz")

    sheet = _sample_sheet(series_path)
    tpm, tpm_symbols = _read_supplement(tpm_path, sheet)
    counts, count_symbols = _read_supplement(counts_path, sheet)

    if not tpm.index.isin(counts.index).all():
        raise ValueError("TPM file contains genes absent from the counts file.")
    # The two files must agree on gene *names* for the same Ensembl id.
    # They do for all but one row: the counts file names ENSG00000073282
    # "KDM5A", while the TPM file (and Ensembl) name it TP63 -- KDM5A is
    # ENSG00000073614, present separately in both files. The id is the
    # authoritative key and the values match TP63's (near-zero, as in the
    # TPM file), so the TPM file's name wins wherever both exist. More
    # than a handful of disagreements would mean a misaligned file rather
    # than a stray label, so that still fails ingest. METHODS.md 11.3.
    paired = count_symbols.reindex(tpm.index)
    name_conflicts = {fid: (tpm_symbols[fid], paired[fid]) for fid in tpm.index[paired != tpm_symbols]}
    if len(name_conflicts) > 5:
        raise ValueError(
            f"TPM and counts files disagree on {len(name_conflicts)} gene names -- "
            "files look misaligned."
        )
    symbols = count_symbols.copy()
    symbols.loc[tpm.index] = tpm_symbols
    n_missing_tpm = int((~counts.index.isin(tpm.index)).sum())

    mycl_means = _check_perturbation(tpm, tpm_symbols, sheet, "TPM")
    _check_perturbation(counts / counts.sum() * 1e6, symbols, sheet, "CPM (counts file)")
    identity = _gene_identity(counts, tpm, symbols)

    contrasts = {}
    for spec in CONTRASTS:
        ref = list(sheet.loc[sheet["arm"] == spec.reference, "sample_id"])
        trt = list(sheet.loc[sheet["arm"] == spec.treatment, "sample_id"])
        table = run_contrast(counts, ref, trt)
        details = dict(table.attrs)
        table = table.copy()
        table.attrs = {}
        ident = identity.reindex(table.index)
        table.insert(0, "symbol", ident["symbol"].to_numpy())
        table.insert(1, "gene_name", ident["gene_name"].to_numpy())
        table.insert(2, "n_loci_with_name", ident["n_loci_with_name"].to_numpy())
        table["in_expression_matrix"] = table.index.isin(tpm.index)
        contrasts[spec.contrast_id] = ContrastResult(spec, table, ref, trt, details)

    samples = pd.DataFrame(
        {
            "sample_id": sheet["sample_id"],
            "dataset_id": DATASET_ID,
            "group_columns": [
                {"cell_line": c, "arm": a, "treatment": t}
                for c, a, t in zip(sheet["cell line"], sheet["arm"], sheet["treatment"])
            ],
            "sample_title": sheet["title"],
            "library_name": sheet["library_name"],
        }
    )

    tested = {k: v.n_genes_tested for k, v in contrasts.items()}
    disagree = {k: _tpm_direction_disagreement(tpm, v) for k, v in contrasts.items()}
    source = DatasetSource(
        dataset_id=DATASET_ID,
        display_name=DISPLAY_NAME,
        accession="GSE324984",
        repository="GEO",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.TPM,
        gene_identifier="ensembl_gene_id",
        n_samples=int(tpm.shape[1]),
        disease_area="Neuroendocrine prostate cancer",
        notes=(
            "PMID 42001796. Two independent experiments, n=2 per arm: C4-2B MYCL "
            "overexpression vs empty vector, and NCI-H660 MYCL knockdown (shMYCL) vs "
            "shCtrl. Compare perturbations only within a cell line -- the cell lines "
            "differ far more than any perturbation (log2-TPM r = 0.89 across lines vs "
            "0.99 between replicates), so cross-line comparisons and correlations "
            "across all 8 samples mostly reflect cell line, not MYCL. Expression shown "
            "is the submitter's TPM; their TPM file is truncated at 65,535 genes "
            f"(Excel's row limit), so {n_missing_tpm:,} of {len(counts):,} genes have "
            "no TPM here. Contrast statistics use the complete counts file: "
            f"{METHOD_LABEL}, fitted per cell line "
            f"(genes tested: C4-2B {tested['c42b_mycl_oe']:,}, "
            f"NCI-H660 {tested['h660_mycl_kd']:,}). With n=2 per arm, Mann-Whitney "
            "p-values elsewhere in the app cannot fall below 0.33. The TPM and counts "
            "files come from different quantifications (TPM is not a fixed-length "
            "rescaling of the counts), so for some genes the mean TPM moves the other "
            "way from the count-based fold-change: "
            + ", ".join(
                f"{contrasts[k].spec.stratum_value} {n:,} of the {m:,} tested genes that have TPM ({n / m:.0%})"
                for k, (n, m) in disagree.items()
            )
            + ", mostly small changes; the gene panel flags them and shows the "
            "normalised counts the test used. Label check: MYCL "
            "mean TPM " + ", ".join(f"{k} {v:.1f}" for k, v in mycl_means.items()) + "."
            + (
                " Gene-name conflicts between the two files, resolved to the TPM file's "
                "name: " + "; ".join(f"{fid} TPM '{a}' vs counts '{b}'" for fid, (a, b) in name_conflicts.items()) + "."
                if name_conflicts else ""
            )
        ),
    )
    return Dataset(
        source=source,
        matrix=tpm,
        samples=samples,
        features=_features(tpm, identity),
        contrasts=contrasts,
    )


from app.registry import DatasetDescriptor, register  # noqa: E402

register(
    DatasetDescriptor(
        dataset_id=DATASET_ID,
        display_name=DISPLAY_NAME,
        loader=load,
        group_columns=("arm", "cell_line", "treatment"),
        supports_survival=False,  # cell lines, no clinical follow-up
    )
)
