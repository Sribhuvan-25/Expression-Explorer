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
    # log2(TPM+1) -- a transformed RNA-seq unit. Distinct from
    # LOG2_INTENSITY, which is a MICROARRAY unit: the two are not on a
    # comparable scale and do not come from comparable assays. DepMap was
    # previously labelled LOG2_INTENSITY for want of this value, so the UI
    # showed the same "log2 intensity" string for DepMap (RNA-seq) and
    # GDS4299 (microarray) -- which undercut the "never pooled across
    # assay types or units" guarantee shown right beside it, since two
    # genuinely different units read as identical (caught in QA).
    LOG2_TPM = "log2_tpm"
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


# Minimum paired samples for an aux correlation to mean anything. Defined
# once because three things must agree on it: the correlation itself, the
# usable-feature count shown in the UI, and the feature picker -- a picker
# offering features the correlation would then reject is worse than no
# picker at all.
MIN_AUX_SAMPLES = 3


class AuxLayer(str, Enum):
    """Non-expression measurements attached to an existing dataset's samples.

    These are NOT separate datasets: DepMap's CRISPR gene-effect scores are
    measured on the same cell lines, keyed by the same ModelIDs, as its
    expression matrix. Registering them as their own `dataset_id` would put
    "DepMap CRISPR" in the sidebar next to "DepMap cell line RNA-seq" as if
    they were different cohorts, which is misleading -- they're the same
    186 samples measured a second way. See METHODS.md 8.1.
    """

    CRISPR_GENE_EFFECT = "crispr_gene_effect"
    DRUG_SENSITIVITY = "drug_sensitivity"
    MUTATION_STATUS = "mutation_status"


class AuxMatrix:
    """One auxiliary measurement layer: a (feature x sample) frame plus
    enough provenance to describe what its numbers actually mean.

    `value_label`/`value_description` exist because these layers are not
    interchangeable the way expression matrices are -- a CRISPR gene-effect
    score near -1 means "knocking this gene out kills the line", while a
    drug log-fold-change near -1 means "this compound kills the line", and
    a mutation entry is a boolean. Any UI rendering one of these has to say
    which, so the meaning travels with the data rather than being hardcoded
    per-layer in the frontend.
    """

    def __init__(
        self,
        layer: AuxLayer,
        matrix: pd.DataFrame,
        value_label: str,
        value_description: str,
        source_note: str,
        features: Optional[pd.DataFrame] = None,
    ):
        self.layer = layer
        self.matrix = matrix  # index=feature_id (gene/compound), columns=sample_id
        self.value_label = value_label
        self.value_description = value_description
        self.source_note = source_note
        # Optional per-feature metadata (e.g. compound id -> drug name for
        # the PRISM layer, where the matrix index is a Broad compound id
        # that means nothing to a reader on its own).
        self.features = features.set_index("feature_id", drop=False) if features is not None else None

    @property
    def n_features(self) -> int:
        return int(self.matrix.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.matrix.shape[1])

    def n_usable_features(self, sample_ids: list[str], min_samples: int = MIN_AUX_SAMPLES) -> int:
        """How many features have enough non-null values among `sample_ids`
        to actually be analysed.

        `n_features` alone oversells a sparse layer badly: PRISM advertises
        6,790 compounds, but against DepMap's 79 covered lymphoid lines
        only ~1,518 (22%) were assayed on at least 3 of them -- the rest
        return "not enough samples" on any lookup. A user told "6,790
        compounds" and then failing on 3 of 4 random picks has been
        misled by a true-but-useless number (caught in QA). CRISPR has no
        such problem (100% usable), which is exactly why this has to be
        measured per layer rather than assumed.
        """
        present = [s for s in sample_ids if s in self.matrix.columns]
        if not present:
            return 0
        return int((self.matrix[present].notna().sum(axis=1) >= min_samples).sum())


class Dataset:
    """In-memory bundle: matrix (features x samples) + both metadata tables,
    plus any auxiliary measurement layers on the same samples.

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
        aux: Optional[dict[AuxLayer, AuxMatrix]] = None,
    ):
        self.source = source
        self.matrix = matrix  # index=feature_id, columns=sample_id
        self.samples = samples.set_index("sample_id", drop=False)
        self.features = features.set_index("feature_id", drop=False)
        # Lazily populated by loaders that have a second measurement layer;
        # every existing loader passes nothing and is unaffected.
        self.aux: dict[AuxLayer, AuxMatrix] = aux or {}

    def require_quantitative(self) -> None:
        if self.source.assay_type == AssayType.QPCR:
            raise ValueError(
                f"{self.source.display_name} is qPCR (targeted assay) — "
                "cannot be queried for genes outside its designed panel."
            )

    def require_aux(self, layer: AuxLayer) -> AuxMatrix:
        """Fetch an auxiliary layer or fail with a message naming what's
        actually available -- an analysis asking for a layer this dataset
        doesn't carry is a caller error worth surfacing clearly, not a
        KeyError."""
        if layer not in self.aux:
            available = ", ".join(sorted(a.value for a in self.aux)) or "none"
            raise ValueError(
                f"{self.source.display_name} has no '{layer.value}' layer "
                f"(available: {available})."
            )
        return self.aux[layer]

    def aux_sample_overlap(self, layer: AuxLayer) -> tuple[int, int]:
        """(n_samples_with_this_layer, n_samples_in_dataset).

        Every aux layer added so far covers only part of its dataset --
        TARGET's mutation calls reach 279 of 469 RNA-seq samples, PRISM's
        drug screen covers a subset of DepMap lines. Reporting that split
        is mandatory for the same reason it is on grouping columns
        (METHODS.md 6.1): an analysis silently running on 60% of a cohort
        reads as the whole cohort unless it says otherwise.
        """
        aux = self.require_aux(layer)
        shared = set(aux.matrix.columns) & set(self.matrix.columns)
        return len(shared), int(self.matrix.shape[1])
