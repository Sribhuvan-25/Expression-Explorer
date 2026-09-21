"""
Dataset metadata schema.

**What lives here and what deliberately does not.** These tables hold the
*description* of a dataset -- its identity, provenance, per-sample
clinical fields, and gene catalogue. The expression values themselves
stay in parquet files, and a row here only points at them
(`datasets.matrix_uri`).

That split is not a matter of taste. Measured on the current data:

| matrix           | shape        | cells  |
|------------------|--------------|--------|
| TARGET expr      | 60660 x 469  | 28.4M  |
| DepMap expr      | 19193 x 1673 | 32.1M  |
| DepMap CRISPR    | 17916 x 1178 | 21.1M  |
| PRISM drug       | 6790 x 919   | 6.2M   |
| TARGET mutations | 6250 x 717   | 4.5M   |

92M dense numeric cells today, and every analysis reads whole rows or
whole columns of them -- a genome-wide scan touches all 60,660 genes and
currently costs 0.03s against parquet. Stored relationally that becomes
92M rows now and 600M+ as datasets are added, turning each scan into a
60,000-row aggregate. Row storage buys nothing for dense numeric data
whose access pattern is "read it all"; it costs an order of magnitude.

So: small and relational goes in Postgres, large and dense stays in
files. This schema is the small half.

**Why this exists at all.** Dataset identity used to be hardcoded --
`registry.ensure_loaded()` held a literal import list, and each
`ingest/*.py` carried its own accession and group columns. Adding a
dataset therefore meant editing code and redeploying. Moving the
description into a table is what makes "add a dataset" a data operation
instead of a code change.
"""
from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import (
    Boolean,
    DateTime,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from sqlalchemy.types import JSON


class Base(DeclarativeBase):
    pass


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


class DatasetRow(Base):
    """One row per dataset -- the replacement for the hardcoded registry."""

    __tablename__ = "datasets"

    dataset_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    display_name: Mapped[str] = mapped_column(String(255))
    accession: Mapped[str] = mapped_column(String(255))
    repository: Mapped[str] = mapped_column(String(128))
    assay_type: Mapped[str] = mapped_column(String(32))
    # Stored as the enum's string value rather than a native DB enum: a new
    # unit must not require a schema migration, and this column already
    # gained `log2_tpm` once (METHODS.md 7.2) when DepMap was found sharing
    # the microarray label. Validated against ExpressionUnit on read.
    expression_unit: Mapped[str] = mapped_column(String(32))
    gene_identifier: Mapped[str] = mapped_column(String(32), default="symbol")
    n_samples: Mapped[int] = mapped_column(Integer)
    supports_survival: Mapped[bool] = mapped_column(Boolean, default=False)
    # Which columns this dataset exposes for grouping. A list, not a
    # separate table: it is read whole, never queried into, and is
    # meaningless outside its dataset.
    group_columns: Mapped[list] = mapped_column(JSON, default=list)
    # Provenance caveats shown in the UI (e.g. DepMap's release-staleness
    # disclosure, METHODS.md 7.5). Free text by nature.
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Where the expression matrix actually is. A URI, not a path, so the
    # same column addresses a local volume today ("file:///data/cache/...")
    # and object storage later ("s3://bucket/...") -- moving the files
    # becomes a config change rather than a schema change.
    matrix_uri: Mapped[str] = mapped_column(String(1024))

    # draft   -> row exists, matrix not yet written (ingest in progress)
    # ready    -> usable
    # failed   -> ingest failed; kept so the failure is visible rather than
    #             the dataset silently missing
    status: Mapped[str] = mapped_column(String(16), default="draft")
    status_detail: Mapped[str | None] = mapped_column(Text, nullable=True)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=_utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=_utcnow, onupdate=_utcnow
    )

    samples: Mapped[list["SampleRow"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )
    features: Mapped[list["FeatureRow"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )
    aux_layers: Mapped[list["AuxLayerRow"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )


class SampleRow(Base):
    """One row per sample per dataset -- the genuinely queryable part.

    `attributes` is a single JSON column rather than typed columns because
    the cohorts really do differ: TARGET carries vital_status /
    days_to_death / etp_status / mrd_status, DepMap carries lineage /
    subtype / cell_line_name, GDS4299 carries only etp_status. Modelling
    those as real columns would mean a schema migration for every new
    cohort -- which defeats the entire point of making dataset addition a
    data operation.

    The tradeoff is real and worth stating: we give up column-level type
    checking and NOT NULL constraints. Postgres can still index and filter
    inside JSONB, so query capability is retained; what is lost is the
    database enforcing shape. Validation therefore happens at ingest time
    (see scripts/ingest_dataset.py) rather than at write time.
    """

    __tablename__ = "samples"
    __table_args__ = (
        UniqueConstraint("dataset_id", "sample_id", name="uq_sample_per_dataset"),
        Index("ix_samples_dataset", "dataset_id"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[str] = mapped_column(
        String(64), ForeignKey("datasets.dataset_id", ondelete="CASCADE")
    )
    sample_id: Mapped[str] = mapped_column(String(128))
    attributes: Mapped[dict] = mapped_column(JSON, default=dict)

    dataset: Mapped[DatasetRow] = relationship(back_populates="samples")


class FeatureRow(Base):
    """Gene catalogue per dataset.

    Exists so "which datasets contain MEF2C" is answerable without loading
    a single matrix -- today that question requires reading every parquet
    file into memory. Largest table by row count (~87k rows across the
    three current datasets) and still trivial for Postgres.
    """

    __tablename__ = "features"
    __table_args__ = (
        UniqueConstraint("dataset_id", "feature_id", name="uq_feature_per_dataset"),
        Index("ix_features_symbol", "symbol"),
        Index("ix_features_dataset", "dataset_id"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[str] = mapped_column(
        String(64), ForeignKey("datasets.dataset_id", ondelete="CASCADE")
    )
    # The dataset's own key -- an Ensembl id for TARGET, a symbol for
    # DepMap/GDS4299. Deliberately not normalised to one id space: the
    # matrix is indexed by this exact value, so rewriting it here would
    # break the join to the file.
    feature_id: Mapped[str] = mapped_column(String(128))
    symbol: Mapped[str] = mapped_column(String(128))
    aliases: Mapped[list] = mapped_column(JSON, default=list)

    dataset: Mapped[DatasetRow] = relationship(back_populates="features")


class AuxLayerRow(Base):
    """An auxiliary measurement layer attached to a dataset's samples.

    Mirrors the AuxLayer/AuxMatrix contract (METHODS.md 8.1): these are
    not separate datasets, they are second measurements on the same
    samples, so they hang off a dataset rather than standing alone.

    Coverage numbers are stored rather than recomputed because they are
    what the UI states up front, and because `n_usable_features` is the
    number that stopped PRISM overselling itself (6,790 compounds
    advertised, 1,518 actually analysable -- METHODS.md 8.4).
    """

    __tablename__ = "aux_layers"
    __table_args__ = (
        UniqueConstraint("dataset_id", "layer", name="uq_layer_per_dataset"),
        Index("ix_aux_dataset", "dataset_id"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[str] = mapped_column(
        String(64), ForeignKey("datasets.dataset_id", ondelete="CASCADE")
    )
    layer: Mapped[str] = mapped_column(String(64))
    value_label: Mapped[str] = mapped_column(String(255))
    value_description: Mapped[str] = mapped_column(Text)
    source_note: Mapped[str] = mapped_column(Text)
    n_features: Mapped[int] = mapped_column(Integer)
    n_usable_features: Mapped[int] = mapped_column(Integer)
    n_samples_covered: Mapped[int] = mapped_column(Integer)
    matrix_uri: Mapped[str] = mapped_column(String(1024))
    # Optional companion table (PRISM's compound -> drug name / MOA frame).
    features_uri: Mapped[str | None] = mapped_column(String(1024), nullable=True)

    dataset: Mapped[DatasetRow] = relationship(back_populates="aux_layers")
