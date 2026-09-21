"""
Ingest a dataset: run its loader, write its matrices to disk, record its
metadata in the database.

This is the command that replaces "edit registry.py and redeploy". After
it runs, the dataset is listed by the API without any code change.

    python scripts/ingest_dataset.py target_all_p2
    python scripts/ingest_dataset.py --all
    python scripts/ingest_dataset.py --all --dry-run

The loader still lives in code (`app/ingest/<name>.py`) -- writing a
parser for a new upstream format is genuinely programming, and pretending
otherwise would just move that work somewhere worse. What moves into data
is the dataset's *description*: identity, provenance, sample attributes,
gene catalogue, layer coverage. That is the part that previously forced a
redeploy.

Ordering is deliberate: matrices are written to their final location
BEFORE the row is marked `ready`. A crash mid-ingest therefore leaves a
`draft`/`failed` row rather than a `ready` row pointing at a file that
does not exist -- the API skips non-ready rows, so a half-finished ingest
is invisible rather than broken.
"""
from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.config import settings
from app.db.models import AuxLayerRow, DatasetRow, FeatureRow, SampleRow
from app.db.session import create_all, session_scope
from app.models.contract import Dataset
from app.registry import ensure_loaded, get_descriptor, list_descriptors

# Columns that live on `samples` for bookkeeping rather than as clinical
# facts about the sample -- excluded from the attributes blob so it stays
# a clean record of what was actually measured/observed.
_NON_ATTRIBUTE_COLUMNS = {"sample_id", "dataset_id", "group_columns"}


def _matrix_dir(dataset_id: str) -> Path:
    return Path(settings.cache_dir) / "datasets" / dataset_id


def _write_parquet(df, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path)
    # A URI, not a bare path: the same column addresses a local volume now
    # and object storage later without a schema change.
    return path.resolve().as_uri()


def _sample_attributes(dataset: Dataset) -> list[dict]:
    """Flatten each sample's clinical fields into one dict.

    Two shapes exist in the current data and both must survive: TARGET
    keeps vital_status/days_to_death as real DataFrame columns, while
    DepMap and GDS4299 keep everything inside a per-sample
    `group_columns` dict. Merging both into one attributes blob is what
    lets a single schema hold cohorts that genuinely differ.
    """
    rows = []
    for record in dataset.samples.to_dict(orient="records"):
        attributes = {
            key: value
            for key, value in record.items()
            if key not in _NON_ATTRIBUTE_COLUMNS and value is not None
        }
        nested = record.get("group_columns") or {}
        if isinstance(nested, dict):
            attributes.update({k: v for k, v in nested.items() if v is not None})
        # pandas NaN is not JSON-serialisable and means "absent" here.
        attributes = {
            k: (None if isinstance(v, float) and v != v else v) for k, v in attributes.items()
        }
        attributes = {k: v for k, v in attributes.items() if v is not None}
        rows.append({"sample_id": str(record["sample_id"]), "attributes": attributes})
    return rows


def ingest(dataset_id: str, dry_run: bool = False) -> None:
    descriptor = get_descriptor(dataset_id)
    print(f"[ingest] {dataset_id}: loading via {descriptor.loader.__module__}...", flush=True)
    start = time.monotonic()
    dataset = descriptor.loader()
    print(
        f"[ingest] {dataset_id}: loaded {dataset.matrix.shape[0]} features x "
        f"{dataset.matrix.shape[1]} samples in {time.monotonic() - start:.0f}s",
        flush=True,
    )

    if dry_run:
        print(f"[ingest] {dataset_id}: dry run, nothing written")
        return

    target_dir = _matrix_dir(dataset_id)
    matrix_uri = _write_parquet(dataset.matrix, target_dir / "matrix.parquet")
    print(f"[ingest] {dataset_id}: matrix -> {matrix_uri}", flush=True)

    aux_records = []
    for layer, aux in dataset.aux.items():
        shared = [c for c in dataset.matrix.columns if c in aux.matrix.columns]
        layer_uri = _write_parquet(aux.matrix, target_dir / f"aux_{layer.value}.parquet")
        features_uri = None
        if aux.features is not None:
            features_uri = _write_parquet(
                aux.features.reset_index(drop=True), target_dir / f"aux_{layer.value}_features.parquet"
            )
        aux_records.append(
            AuxLayerRow(
                dataset_id=dataset_id,
                layer=layer.value,
                value_label=aux.value_label,
                value_description=aux.value_description,
                source_note=aux.source_note,
                n_features=aux.n_features,
                n_usable_features=aux.n_usable_features(shared),
                n_samples_covered=len(shared),
                matrix_uri=layer_uri,
                features_uri=features_uri,
            )
        )
        print(f"[ingest] {dataset_id}: aux '{layer.value}' -> {layer_uri}", flush=True)

    source = dataset.source
    with session_scope() as session:
        # Replace wholesale rather than diffing: re-ingesting is how a
        # dataset gets refreshed, and a partial update could leave sample
        # rows from an older release mixed with a newer matrix.
        existing = session.get(DatasetRow, dataset_id)
        if existing is not None:
            session.delete(existing)
            session.flush()

        row = DatasetRow(
            dataset_id=dataset_id,
            display_name=descriptor.display_name,
            accession=source.accession,
            repository=source.repository,
            assay_type=source.assay_type.value,
            expression_unit=source.expression_unit.value,
            gene_identifier=getattr(source, "gene_identifier", "symbol"),
            n_samples=int(dataset.matrix.shape[1]),
            supports_survival=descriptor.supports_survival,
            group_columns=list(descriptor.group_columns),
            notes=source.notes,
            matrix_uri=matrix_uri,
            status="ready",
        )
        row.samples = [
            SampleRow(dataset_id=dataset_id, **rec) for rec in _sample_attributes(dataset)
        ]
        row.features = [
            FeatureRow(
                dataset_id=dataset_id,
                feature_id=str(rec["feature_id"]),
                symbol=str(rec.get("symbol") or rec["feature_id"]),
                aliases=list(rec.get("aliases") or []),
            )
            for rec in dataset.features.to_dict(orient="records")
        ]
        row.aux_layers = aux_records
        session.add(row)

    print(
        f"[ingest] {dataset_id}: recorded "
        f"({len(dataset.samples)} samples, {len(dataset.features)} features, "
        f"{len(aux_records)} aux layers)",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset_id", nargs="?", help="dataset to ingest")
    parser.add_argument("--all", action="store_true", help="ingest every registered dataset")
    parser.add_argument("--dry-run", action="store_true", help="load and report, write nothing")
    parser.add_argument("--list", action="store_true", help="list registered dataset ids")
    args = parser.parse_args()

    ensure_loaded()
    create_all()

    if args.list:
        for d in list_descriptors():
            print(f"  {d.dataset_id:20} {d.display_name}")
        return 0

    if args.all:
        targets = [d.dataset_id for d in list_descriptors()]
    elif args.dataset_id:
        targets = [args.dataset_id]
    else:
        parser.error("give a dataset_id, or --all, or --list")
        return 2

    failures = []
    for dataset_id in targets:
        try:
            ingest(dataset_id, dry_run=args.dry_run)
        except Exception as exc:  # noqa: BLE001 - one bad source must not stop the rest
            failures.append(dataset_id)
            print(f"[ingest] {dataset_id}: FAILED -- {exc}", file=sys.stderr, flush=True)
            traceback.print_exc()
            if not args.dry_run:
                # Record the failure so the dataset is visibly broken
                # rather than silently absent.
                with session_scope() as session:
                    row = session.get(DatasetRow, dataset_id)
                    if row is not None:
                        row.status = "failed"
                        row.status_detail = str(exc)[:1000]

    if failures:
        print(f"\n[ingest] {len(failures)} failed: {', '.join(failures)}", file=sys.stderr)
        return 1
    print(f"\n[ingest] done ({len(targets)} dataset(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
