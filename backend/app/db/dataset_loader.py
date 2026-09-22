"""
Rebuild a `Dataset` from database metadata plus its parquet files.

This is the other half of `scripts/ingest_dataset.py`: ingest writes the
description into the database and the matrices to disk, and this reads
both back into the same in-memory `Dataset` the analysis layer already
expects. Nothing downstream of the registry changes -- the contract is
identical, only its source moved.

The point of the split is that the *ingest* loader (`app/ingest/*.py`)
knows how to talk to GDC or Figshare and parse their formats, while this
one only knows how to read what ingest already produced. A dataset
recorded in the database therefore loads without its original ingest
module being importable at all, which is what makes adding a dataset a
data operation rather than a code change.
"""
from __future__ import annotations

from pathlib import Path
from urllib.parse import unquote, urlparse

import pandas as pd

from app.db.models import DatasetRow
from app.models.contract import (
    AssayType,
    AuxLayer,
    AuxMatrix,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)


def uri_to_path(uri: str) -> Path:
    """Resolve a stored matrix URI to a local path.

    Only file:// is supported today. Object storage (s3://, r2://) would
    be handled here -- the reason matrices are addressed by URI rather
    than bare path is precisely so that becomes a change in this one
    function instead of a schema migration.
    """
    parsed = urlparse(uri)
    if parsed.scheme in ("", "file"):
        return Path(unquote(parsed.path))
    raise ValueError(
        f"Unsupported matrix URI scheme {parsed.scheme!r} in {uri!r} -- "
        "only local file:// URIs are handled so far."
    )


def _read_matrix(uri: str) -> pd.DataFrame:
    path = uri_to_path(uri)
    if not path.exists():
        # A `ready` row whose file is missing is a real inconsistency, not
        # a missing dataset -- say which file, so it is fixable rather
        # than mysterious.
        raise FileNotFoundError(
            f"Matrix file for this dataset is missing: {path}. "
            "Re-run scripts/ingest_dataset.py for it."
        )
    return pd.read_parquet(path)


def _samples_frame(row: DatasetRow, matrix: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct the samples table the analysis layer expects.

    Ingest flattened each sample's clinical fields into one `attributes`
    blob (see METHODS.md 10.3). The analysis layer reads group values
    from BOTH top-level columns and a nested `group_columns` dict, so
    both shapes are restored here: attributes become real columns, and a
    copy goes into `group_columns`. Restoring only one shape would break
    whichever datasets used the other.
    """
    by_id = {s.sample_id: dict(s.attributes or {}) for s in row.samples}
    # Ordered by the matrix, not by the database: the analysis layer
    # indexes samples positionally against matrix columns in places, and
    # a row order that drifts from the matrix would misalign them.
    sample_ids = [str(c) for c in matrix.columns]
    records = []
    for sample_id in sample_ids:
        attributes = by_id.get(sample_id, {})
        record = {
            "sample_id": sample_id,
            "dataset_id": row.dataset_id,
            "group_columns": dict(attributes),
            **attributes,
        }
        records.append(record)
    frame = pd.DataFrame(records)
    # Guarantee every declared grouping column exists even if no sample
    # carries a value, so a group-by on it yields "no values" rather than
    # a KeyError.
    for column in row.group_columns or []:
        if column not in frame.columns:
            frame[column] = None
    return frame


def _features_frame(row: DatasetRow, matrix: pd.DataFrame) -> pd.DataFrame:
    if row.features:
        frame = pd.DataFrame(
            [
                {
                    "feature_id": f.feature_id,
                    "symbol": f.symbol,
                    "aliases": list(f.aliases or []),
                }
                for f in row.features
            ]
        )
    else:
        # Degrade to the matrix index rather than failing: a dataset with
        # no recorded catalogue is still analysable, it just has no alias
        # information.
        frame = pd.DataFrame(
            {
                "feature_id": [str(i) for i in matrix.index],
                "symbol": [str(i) for i in matrix.index],
                "aliases": [[] for _ in matrix.index],
            }
        )
    return frame


def _aux_layers(row: DatasetRow) -> dict[AuxLayer, AuxMatrix]:
    layers: dict[AuxLayer, AuxMatrix] = {}
    for record in row.aux_layers:
        try:
            layer = AuxLayer(record.layer)
        except ValueError:
            # A layer recorded by a newer version of the app than this one
            # understands. Skipping keeps the dataset usable instead of
            # failing the whole load over an unknown extra.
            continue
        features = None
        if record.features_uri:
            features = pd.read_parquet(uri_to_path(record.features_uri))
        layers[layer] = AuxMatrix(
            layer=layer,
            matrix=_read_matrix(record.matrix_uri),
            value_label=record.value_label,
            value_description=record.value_description,
            source_note=record.source_note,
            features=features,
        )
    return layers


def load_from_row(row: DatasetRow) -> Dataset:
    """Materialise the `Dataset` described by one database row."""
    matrix = _read_matrix(row.matrix_uri)
    source = DatasetSource(
        dataset_id=row.dataset_id,
        display_name=row.display_name,
        accession=row.accession,
        repository=row.repository,
        assay_type=AssayType(row.assay_type),
        expression_unit=ExpressionUnit(row.expression_unit),
        gene_identifier=row.gene_identifier,
        n_samples=int(matrix.shape[1]),
        notes=row.notes,
    )
    return Dataset(
        source=source,
        matrix=matrix,
        samples=_samples_frame(row, matrix),
        features=_features_frame(row, matrix),
        aux=_aux_layers(row),
    )
