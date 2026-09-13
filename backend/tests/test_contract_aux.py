import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from app.models.contract import (
    AssayType,
    AuxLayer,
    AuxMatrix,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)


def _dataset(sample_ids=("S1", "S2", "S3"), aux=None):
    matrix = pd.DataFrame(
        {s: [1.0, 2.0] for s in sample_ids}, index=["GENE1", "GENE2"]
    )
    samples = pd.DataFrame(
        {"sample_id": list(sample_ids), "dataset_id": "test", "group_columns": [{} for _ in sample_ids]}
    )
    features = pd.DataFrame(
        {"feature_id": ["GENE1", "GENE2"], "symbol": ["GENE1", "GENE2"], "aliases": [[], []]}
    )
    source = DatasetSource(
        dataset_id="test",
        display_name="Test dataset",
        accession="none",
        repository="test",
        assay_type=AssayType.RNA_SEQ,
        expression_unit=ExpressionUnit.TPM,
        n_samples=len(sample_ids),
    )
    return Dataset(source=source, matrix=matrix, samples=samples, features=features, aux=aux)


def _aux(layer=AuxLayer.CRISPR_GENE_EFFECT, sample_ids=("S1", "S2")):
    return AuxMatrix(
        layer=layer,
        matrix=pd.DataFrame({s: [-0.5, 0.1] for s in sample_ids}, index=["GENE1", "GENE2"]),
        value_label="Gene effect",
        value_description="Chronos score; more negative = stronger dependency.",
        source_note="test",
    )


def test_dataset_without_aux_defaults_to_empty_not_none():
    """Every pre-existing loader constructs Dataset with no aux argument --
    that must keep working and yield an empty dict, not None, so callers
    can iterate without a null check."""
    ds = _dataset()
    assert ds.aux == {}


def test_require_aux_returns_the_layer_when_present():
    aux_layer = _aux()
    ds = _dataset(aux={AuxLayer.CRISPR_GENE_EFFECT: aux_layer})
    assert ds.require_aux(AuxLayer.CRISPR_GENE_EFFECT) is aux_layer


def test_require_aux_missing_layer_names_what_is_available():
    """A caller asking for a layer this dataset doesn't carry should get a
    message naming what it *does* carry, not a bare KeyError."""
    ds = _dataset(aux={AuxLayer.CRISPR_GENE_EFFECT: _aux()})
    with pytest.raises(ValueError, match="crispr_gene_effect"):
        ds.require_aux(AuxLayer.DRUG_SENSITIVITY)


def test_require_aux_on_dataset_with_no_layers_says_none_available():
    ds = _dataset()
    with pytest.raises(ValueError, match="available: none"):
        ds.require_aux(AuxLayer.MUTATION_STATUS)


def test_aux_sample_overlap_reports_partial_coverage():
    """The whole reason this method exists: every aux layer covers only
    part of its dataset, and that split has to be reportable (METHODS.md
    6.1) rather than silently narrowing an analysis."""
    ds = _dataset(sample_ids=("S1", "S2", "S3"), aux={AuxLayer.CRISPR_GENE_EFFECT: _aux(sample_ids=("S1", "S2"))})
    covered, total = ds.aux_sample_overlap(AuxLayer.CRISPR_GENE_EFFECT)
    assert (covered, total) == (2, 3)


def test_aux_sample_overlap_ignores_aux_samples_not_in_the_dataset():
    """An aux source may cover samples the expression matrix doesn't have
    (DepMap's drug screen spans far more cell lines than our lymphoid
    filter keeps) -- overlap must count the intersection, not the aux
    layer's own width."""
    ds = _dataset(sample_ids=("S1", "S2"), aux={AuxLayer.CRISPR_GENE_EFFECT: _aux(sample_ids=("S1", "S9", "S8"))})
    covered, total = ds.aux_sample_overlap(AuxLayer.CRISPR_GENE_EFFECT)
    assert (covered, total) == (1, 2)


def test_aux_matrix_reports_its_own_dimensions():
    aux_layer = _aux(sample_ids=("S1", "S2", "S3"))
    assert aux_layer.n_features == 2
    assert aux_layer.n_samples == 3
