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

**Datasets now come from the database first (2026-09-21).** Every row in
`datasets` with `status='ready'` becomes a descriptor whose loader reads
the parquet recorded in its `matrix_uri` -- so adding a dataset is
`python scripts/ingest_dataset.py <name>`, with no redeploy and no edit
to this file. See METHODS.md §10.

Code-registered datasets remain, and are the fallback. A dataset that is
in code but not yet ingested still works, which is what makes the
cutover safe: an empty or unreachable database degrades to exactly the
previous behaviour rather than emptying the app. Where both exist, the
database wins, since that is the row an operator most recently wrote.

To add a dataset: write an ingest/<name>.py with a `load()` function
(see ingest/gdc_target.py for the shape) and one `register(...)` call,
then run the ingest script once. The ingest module is only needed to
parse the upstream format -- after ingest, the dataset loads from its
recorded parquet without it.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Callable

from app.models.contract import AssayType, Dataset

logger = logging.getLogger(__name__)

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


def _descriptor_from_row(row) -> DatasetDescriptor:
    """Build a descriptor whose loader reads the recorded parquet.

    `dataset_id` is bound by value rather than closing over `row`: the
    row belongs to a session that is closed by the time the loader runs,
    and touching a detached ORM object then raises. Re-fetching by id
    inside the loader also means a re-ingest is picked up on the next
    load instead of serving a stale description.
    """
    from app.db.dataset_loader import load_from_row
    from app.db.models import DatasetRow
    from app.db.session import session_scope

    dataset_id = row.dataset_id

    def loader():
        with session_scope() as session:
            fresh = session.get(DatasetRow, dataset_id)
            if fresh is None:
                raise KeyError(f"Dataset '{dataset_id}' is no longer in the database.")
            # Relationships are lazy-loaded, so they must be realised
            # while the session is still open.
            _ = fresh.samples, fresh.features, fresh.aux_layers
            return load_from_row(fresh)

    return DatasetDescriptor(
        dataset_id=dataset_id,
        display_name=row.display_name,
        loader=loader,
        group_columns=tuple(row.group_columns or ()),
        supports_survival=bool(row.supports_survival),
    )


def _db_descriptors() -> dict[str, DatasetDescriptor]:
    """Ready datasets from the database, or {} if it cannot be read.

    Deliberately swallowing every error: a database that is unreachable,
    not yet migrated, or empty must degrade to the code-registered
    datasets rather than taking the API down. The previous behaviour is
    the floor, never the casualty.
    """
    try:
        from app.db.models import DatasetRow
        from app.db.session import session_scope
        from sqlalchemy import select

        with session_scope() as session:
            rows = session.scalars(
                select(DatasetRow).where(DatasetRow.status == "ready")
            ).all()
            return {row.dataset_id: _descriptor_from_row(row) for row in rows}
    except Exception as exc:  # noqa: BLE001 - see docstring
        logger.warning(
            "Dataset metadata database unavailable (%s); "
            "falling back to code-registered datasets.",
            exc,
        )
        return {}


def _merged() -> dict[str, DatasetDescriptor]:
    merged = dict(_REGISTRY)
    merged.update(_db_descriptors())  # database wins on conflict
    return merged


def get_descriptor(dataset_id: str) -> DatasetDescriptor:
    merged = _merged()
    if dataset_id not in merged:
        raise KeyError(dataset_id)
    return merged[dataset_id]


def list_descriptors() -> list[DatasetDescriptor]:
    return list(_merged().values())


def ensure_loaded() -> None:
    """Import every ingest module so its register() call runs. Called
    once at API startup; new ingest modules just need to be imported
    here to join the registry — no other wiring."""
    from app.ingest import depmap, gdc_target, gds4299  # noqa: F401
