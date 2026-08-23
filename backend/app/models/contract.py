"""
The dataset contract. Every ingest loader normalizes its source into these
three objects; every analysis function consumes only these, never a
source-specific shape. This is what makes the tool adaptive across datasets.
"""
from __future__ import annotations

from enum import Enum
from typing import Optional

import pandas as pd
from pydantic import BaseModel, Field


class AssayType(str, Enum):
    RNA_SEQ = "rna_seq"
    MICROARRAY = "microarray"
    QPCR = "qpcr"
    SC_RNA_SEQ = "sc_rna_seq"


class ExpressionUnit(str, Enum):
    TPM = "tpm"
    FPKM = "fpkm"
    COUNTS = "counts"
    LOG2_INTENSITY = "log2_intensity"
    DELTA_CT = "delta_ct"


class DatasetSource(BaseModel):
    """Provenance metadata for one ingested dataset."""

    dataset_id: str
    display_name: str
    accession: str
    repository: str  # "GDC" | "GEO" | "DepMap"
    assay_type: AssayType
    expression_unit: ExpressionUnit
    gene_identifier: str = "ensembl_gene_id"
    n_samples: int
    disease_area: str = "ETP-ALL"
    notes: Optional[str] = None


class SampleMetadata(BaseModel):
    """One row per sample; arbitrary extra clinical/annotation columns
    are allowed via `extra` so new disease areas add columns, not code."""

    sample_id: str
    dataset_id: str
    group_columns: dict = Field(default_factory=dict)
    vital_status: Optional[str] = None
    days_to_death: Optional[float] = None
    days_to_last_follow_up: Optional[float] = None


class FeatureMetadata(BaseModel):
    """One row per gene/feature."""

    feature_id: str  # canonical Ensembl ID
    symbol: str
    aliases: list[str] = Field(default_factory=list)


class Dataset:
    """In-memory bundle: matrix (features x samples) + both metadata tables.

    An assay-type guard lives here rather than in the UI: any analysis that
    requires genome-wide expression should call `require_quantitative()`
    before running, so a qPCR-only dataset fails loudly instead of silently
    returning a meaningless answer (see PDF p.32 comment on qPCR data).
    """

    def __init__(
        self,
        source: DatasetSource,
        matrix: pd.DataFrame,
        samples: pd.DataFrame,
        features: pd.DataFrame,
    ):
        self.source = source
        self.matrix = matrix  # index=feature_id, columns=sample_id
        self.samples = samples.set_index("sample_id", drop=False)
        self.features = features.set_index("feature_id", drop=False)

    def require_quantitative(self) -> None:
        if self.source.assay_type == AssayType.QPCR:
            raise ValueError(
                f"{self.source.display_name} is qPCR (targeted assay) — "
                "cannot be queried for genes outside its designed panel."
            )
