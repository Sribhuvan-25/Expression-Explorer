"""
Dataset registry — the single place new datasets get plugged in.

Before this module existed, adding a dataset meant editing a hardcoded
Literal type and an if/elif chain inside app/api/main.py. That's the
opposite of "adaptive": every new dataset touched API dispatch code.

Now: a loader module calls `register()` once, at import time, with a
DatasetDescriptor. The API and frontend both work purely off this
registry — `main.py` never names a dataset_id in code, and the frontend
builds its dataset picker and grouping-column options from what
`/datasets/{id}` reports, not from a hardcoded array.

To add a dataset: write an ingest/<name>.py with a `load()` function
(see ingest/gdc_target.py for the shape), then add one `register(...)`
call at the bottom of that file. Nothing else changes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from app.models.contract import AssayType, Dataset

# Analyses that make sense depend on assay type: qPCR data can't support
# arbitrary-gene queries (see Dataset.require_quantitative), and a dataset
# with no clinical follow-up columns can't support survival. Each loader
# declares which analyses it actually supports so the frontend doesn't
# have to special-case datasets by id.
ANALYSES_BY_ASSAY: dict[AssayType, tuple[str, ...]] = {
    AssayType.RNA_SEQ: ("compare", "signature_score", "survival"),
    AssayType.MICROARRAY: ("compare", "signature_score", "survival"),
    AssayType.SC_RNA_SEQ: ("compare", "signature_score"),
    AssayType.QPCR: (),  # require_quantitative() blocks genome-wide queries anyway
}


@dataclass(frozen=True)
class DatasetDescriptor:
    dataset_id: str
    display_name: str
    loader: Callable[[], Dataset]
    # Group columns this dataset is known to expose in `samples`, either as
    # top-level columns (e.g. "vital_status") or keys inside each sample's
    # group_columns dict (e.g. "lineage" for DepMap). Declared by the
    # loader author since it's the one place that actually knows what its
    # source data contains.
    group_columns: tuple[str, ...] = field(default_factory=tuple)
    supports_survival: bool = False


_REGISTRY: dict[str, DatasetDescriptor] = {}


def register(descriptor: DatasetDescriptor) -> None:
    if descriptor.dataset_id in _REGISTRY:
        raise ValueError(f"Dataset '{descriptor.dataset_id}' is already registered.")
    _REGISTRY[descriptor.dataset_id] = descriptor


def get_descriptor(dataset_id: str) -> DatasetDescriptor:
    if dataset_id not in _REGISTRY:
        raise KeyError(dataset_id)
    return _REGISTRY[dataset_id]


def list_descriptors() -> list[DatasetDescriptor]:
    return list(_REGISTRY.values())


def ensure_loaded() -> None:
    """Import every ingest module so its register() call runs. Called
    once at API startup; new ingest modules just need to be imported
    here to join the registry — no other wiring."""
    from app.ingest import depmap, gdc_target, gds4299  # noqa: F401
