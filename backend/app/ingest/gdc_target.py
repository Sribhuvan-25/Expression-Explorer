"""
Loader for TARGET-ALL-P2 (GDC). Confirmed open-access: 532 STAR-Counts
files, 1,587 cases with clinical follow-up. This is the only module that
knows about GDC's API shape — everything downstream sees a Dataset.
"""
from __future__ import annotations

import gzip
import io

import httpx
import pandas as pd

from app.config import settings
from app.ingest.liu2017_etp import load_etp_status
from app.ingest.target_mrd import load_mrd_status
from app.models.contract import (
    AssayType,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)

GDC_API = "https://api.gdc.cancer.gov"
PROJECT_ID = "TARGET-ALL-P2"
CACHE_DIR = settings.cache_dir / "gdc_target"


def _client() -> httpx.Client:
    return httpx.Client(base_url=GDC_API, timeout=60.0)


def list_open_rna_files(client: httpx.Client) -> list[dict]:
    """All open-access Gene Expression Quantification files for the project."""
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
                {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
                {"op": "in", "content": {"field": "access", "value": ["open"]}},
            ],
        },
        "fields": "file_id,file_name,cases.submitter_id,cases.case_id",
        "size": 1000,
    }
    resp = client.post("/files", json=payload)
    resp.raise_for_status()
    return resp.json()["data"]["hits"]


def fetch_clinical(client: httpx.Client) -> pd.DataFrame:
    """Vital status + follow-up times for every case, for survival analysis."""
    payload = {
        "filters": {
            "op": "in",
            "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]},
        },
        "fields": ",".join(
            [
                "submitter_id",
                "demographic.vital_status",
                "demographic.days_to_death",
                "diagnoses.days_to_last_follow_up",
            ]
        ),
        "size": 2000,
    }
    resp = client.post("/cases", json=payload)
    resp.raise_for_status()
    rows = []
    for hit in resp.json()["data"]["hits"]:
        dem = hit.get("demographic") or {}
        diag = (hit.get("diagnoses") or [{}])[0]
        rows.append(
            {
                "sample_id": hit["submitter_id"],
                "vital_status": dem.get("vital_status"),
                "days_to_death": dem.get("days_to_death"),
                "days_to_last_follow_up": diag.get("days_to_last_follow_up"),
            }
        )
    return pd.DataFrame(rows)


def _download_star_counts(client: httpx.Client, file_id: str) -> pd.DataFrame:
    """One STAR-Counts TSV -> a 2-column frame (tpm, gene_name) indexed by
    versioned Ensembl gene id.

    GDC's STAR-Counts output carries a '# gene-model:' comment line and
    N_unmapped/N_multi/N_noFeature/N_ambiguous summary rows — all dropped
    before use. gene_name is carried through here (not derived from the
    Ensembl id later) since GDC ships the authoritative GENCODE symbol
    alongside each row.
    """
    resp = client.get(f"/data/{file_id}", follow_redirects=True)
    resp.raise_for_status()
    raw = resp.content
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lines = raw.decode().splitlines()
    lines = [ln for ln in lines if not ln.startswith("# ")]
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
    df = df[~df["gene_id"].str.startswith("N_")]
    df = df.set_index("gene_id")
    return df[["gene_name", "tpm_unstranded"]]


def load(limit: int | None = None, use_cache: bool = True) -> Dataset:
    """Build the TARGET-ALL-P2 Dataset. Caches the assembled matrix to
    Parquet so repeat runs don't re-hit the GDC API."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    matrix_path = CACHE_DIR / "matrix.parquet"
    clinical_path = CACHE_DIR / "clinical.parquet"

    if use_cache and matrix_path.exists() and clinical_path.exists():
        matrix = pd.read_parquet(matrix_path)
        clinical = pd.read_parquet(clinical_path)
    else:
        with _client() as client:
            files = list_open_rna_files(client)
            if limit:
                files = files[:limit]
            columns = {}
            gene_names: pd.Series | None = None
            for f in files:
                case = (f.get("cases") or [{}])[0]
                sample_id = case.get("submitter_id", f["file_id"])
                gene_frame = _download_star_counts(client, f["file_id"])
                columns[sample_id] = gene_frame["tpm_unstranded"]
                if gene_names is None:
                    gene_names = gene_frame["gene_name"]
            matrix = pd.DataFrame(columns)
            clinical = fetch_clinical(client)
            gene_names.to_frame("gene_name").to_parquet(CACHE_DIR / "gene_names.parquet")
        matrix.to_parquet(matrix_path)
        clinical.to_parquet(clinical_path)

    # ETP/near-ETP/non-ETP status isn't in any GDC clinical file for this
    # project (checked all 9 clinical supplements) — sourced separately
    # from Liu et al. 2017, the genomic characterization paper for this
    # exact TARGET T-ALL cohort. Only ~190 of the ~530 RNA-seq samples get
    # a classification; the rest simply lacked immunophenotyping and stay
    # unlabeled rather than guessed at.
    etp_status = load_etp_status()
    # Day-29 MRD status, same source and threshold Wang et al. 2025 used
    # (TARGET's Phase II Validation clinical supplement, MRD_neg <=0.01).
    mrd_status = load_mrd_status()

    def _group_columns(sid: str) -> dict:
        cols = {}
        if sid in etp_status.index:
            cols["etp_status"] = etp_status[sid]
        if sid in mrd_status.index:
            cols["mrd_status"] = mrd_status[sid]
        return cols

    samples = clinical.copy()
    samples["dataset_id"] = "target_all_p2"
    samples["group_columns"] = samples["sample_id"].map(_group_columns)
    samples = samples[samples["sample_id"].isin(matrix.columns)].reset_index(drop=True)
    matrix = matrix[samples["sample_id"].tolist()]

    gene_names_path = CACHE_DIR / "gene_names.parquet"
    gene_names = pd.read_parquet(gene_names_path)["gene_name"] if gene_names_path.exists() else pd.Series(dtype=str)
    features = pd.DataFrame({"feature_id": matrix.index})
    features["symbol"] = features["feature_id"].map(gene_names).fillna(features["feature_id"])
    features["aliases"] = [[] for _ in range(len(features))]

    source = DatasetSource(
        dataset_id="target_all_p2",
        display_name="TARGET ALL-P2 (paediatric T-ALL)",
        accession="phs000218 / phs000463 (TARGET-ALL-P2)",
        repository="GDC",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.TPM,
        n_samples=matrix.shape[1],
        notes="Open-access STAR-Counts RNA-seq, one uniform workflow.",
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features)


from app.registry import DatasetDescriptor, register  # noqa: E402

register(
    DatasetDescriptor(
        dataset_id="target_all_p2",
        display_name="TARGET ALL-P2 (paediatric T-ALL)",
        loader=load,
        group_columns=("vital_status", "etp_status", "mrd_status"),
        supports_survival=True,
    )
)
