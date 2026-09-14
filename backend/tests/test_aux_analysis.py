import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from app.analysis.aux_analysis import (
    expression_by_mutation_status,
    expression_vs_aux_correlation,
    most_mutated_genes,
)
from app.models.contract import (
    AssayType,
    AuxLayer,
    AuxMatrix,
    Dataset,
    DatasetSource,
    ExpressionUnit,
)


def _dataset(expression: pd.DataFrame, aux_layers=None):
    sample_ids = list(expression.columns)
    samples = pd.DataFrame(
        {"sample_id": sample_ids, "dataset_id": "test", "group_columns": [{} for _ in sample_ids]}
    )
    features = pd.DataFrame(
        {"feature_id": list(expression.index), "symbol": list(expression.index),
         "aliases": [[] for _ in expression.index]}
    )
    source = DatasetSource(
        dataset_id="test", display_name="Test", accession="none", repository="test",
        assay_type=AssayType.RNA_SEQ, expression_unit=ExpressionUnit.TPM, n_samples=len(sample_ids),
    )
    return Dataset(source=source, matrix=expression, samples=samples, features=features, aux=aux_layers or {})


def _mutation_layer(frame: pd.DataFrame):
    return AuxMatrix(
        layer=AuxLayer.MUTATION_STATUS, matrix=frame,
        value_label="Non-silent mutation present", value_description="bool", source_note="test",
    )


def _numeric_layer(frame: pd.DataFrame, layer=AuxLayer.CRISPR_GENE_EFFECT):
    return AuxMatrix(
        layer=layer, matrix=frame, value_label="Gene effect",
        value_description="more negative = stronger dependency", source_note="test",
    )


# --- mutation-stratified expression ---------------------------------------


def test_expression_by_mutation_status_splits_and_tests():
    expression = pd.DataFrame(
        {f"S{i}": [10.0 if i < 4 else 1.0] for i in range(8)}, index=["MYCN"]
    )
    mutations = pd.DataFrame(
        {f"S{i}": [i < 4] for i in range(8)}, index=["NOTCH1"]
    )
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    result = expression_by_mutation_status(ds, "MYCN", "NOTCH1")

    assert result["n_mutated"] == 4
    assert result["n_wildtype"] == 4
    assert {p["group"] for p in result["points"]} == {"Mutated", "Wild-type"}
    # Perfectly separated groups -> the test should fire.
    assert result["test"]["p_value"] is not None
    assert result["test"]["p_value"] < 0.05


def test_expression_by_mutation_status_reports_partial_coverage():
    """The headline caveat for this layer: mutation data covers only part
    of the RNA-seq cohort, and that split must be reported, not hidden."""
    expression = pd.DataFrame({f"S{i}": [5.0] for i in range(10)}, index=["MYCN"])
    # Only 6 of the 10 expression samples have mutation calls.
    mutations = pd.DataFrame({f"S{i}": [i % 2 == 0] for i in range(6)}, index=["NOTCH1"])
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    result = expression_by_mutation_status(ds, "MYCN", "NOTCH1")

    assert result["n_dataset_total"] == 10
    assert result["n_excluded"] == 4
    assert "mutation data" in result["exclusion_reason"]


def test_expression_by_mutation_status_allows_same_gene_both_sides():
    """Unlike gene-vs-gene correlation, splitting a gene's own expression by
    its own mutation status is a real question, not a user error."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(6)}, index=["NOTCH1"])
    mutations = pd.DataFrame({f"S{i}": [i < 3] for i in range(6)}, index=["NOTCH1"])
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    result = expression_by_mutation_status(ds, "NOTCH1", "NOTCH1")
    assert result["n_mutated"] == 3 and result["n_wildtype"] == 3


def test_expression_by_mutation_status_unknown_gene_raises_value_error():
    expression = pd.DataFrame({f"S{i}": [5.0] for i in range(6)}, index=["MYCN"])
    mutations = pd.DataFrame({f"S{i}": [True] for i in range(6)}, index=["NOTCH1"])
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    with pytest.raises(ValueError, match="no recorded non-silent mutation"):
        expression_by_mutation_status(ds, "MYCN", "NOTAGENE")


def test_expression_by_mutation_status_requires_the_layer():
    expression = pd.DataFrame({f"S{i}": [5.0] for i in range(6)}, index=["MYCN"])
    ds = _dataset(expression)
    with pytest.raises(ValueError, match="mutation_status"):
        expression_by_mutation_status(ds, "MYCN", "NOTCH1")


# --- most-mutated genes ---------------------------------------------------


def test_most_mutated_genes_ranks_by_sample_count():
    expression = pd.DataFrame({f"S{i}": [1.0] for i in range(10)}, index=["MYCN"])
    mutations = pd.DataFrame(
        {f"S{i}": [i < 7, i < 3, False] for i in range(10)},
        index=["COMMON", "RARE", "NEVER"],
    )
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    result = most_mutated_genes(ds, top_n=10)
    genes = result["genes"]

    assert genes[0]["gene"] == "COMMON" and genes[0]["n_mutated"] == 7
    assert genes[1]["gene"] == "RARE" and genes[1]["n_mutated"] == 3
    # A gene mutated in nobody is not a useful thing to offer as a split.
    assert all(g["gene"] != "NEVER" for g in genes)
    assert genes[0]["fraction"] == pytest.approx(0.7)


# --- expression vs. numeric aux layer -------------------------------------


def test_expression_vs_aux_correlation_matches_scipy():
    rng = np.random.default_rng(0)
    sample_ids = [f"ACH-{i:06d}" for i in range(20)]
    expr_vals = rng.normal(size=20)
    aux_vals = expr_vals * -2 + rng.normal(scale=0.1, size=20)

    expression = pd.DataFrame([expr_vals], index=["MYCN"], columns=sample_ids)
    effects = pd.DataFrame([aux_vals], index=["MYCN"], columns=sample_ids)
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    result = expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "MYCN")

    expected_r, expected_p = stats.pearsonr(expr_vals, aux_vals)
    assert result["coefficient"] == pytest.approx(expected_r)
    assert result["p_value"] == pytest.approx(expected_p)
    assert result["n"] == 20
    # Strongly anti-correlated by construction.
    assert result["coefficient"] < -0.9


def test_expression_vs_aux_correlation_reports_partial_overlap():
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame({f"S{i}": [float(i)] for i in range(4)}, index=["MYCN"])
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    result = expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "MYCN")
    assert result["n"] == 4
    assert result["n_dataset_total"] == 10
    assert result["n_excluded"] == 6


def test_expression_vs_aux_correlation_carries_the_value_meaning():
    """A CRISPR score and a drug log-fold-change are both 'more negative is
    more effect', which is the opposite of expression -- the label has to
    reach the caller or a UI renders an unlabelled axis."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(6)}, index=["MYCN"])
    effects = pd.DataFrame({f"S{i}": [float(-i)] for i in range(6)}, index=["MYCN"])
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    result = expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "MYCN")
    assert result["value_label"] == "Gene effect"
    assert "negative" in result["value_description"]


def test_expression_vs_aux_correlation_too_few_shared_samples_raises():
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame({"S0": [1.0], "S1": [2.0]}, index=["MYCN"])
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    with pytest.raises(ValueError, match="at least 3"):
        expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "MYCN")


def test_expression_vs_aux_correlation_unknown_aux_feature_raises():
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(6)}, index=["MYCN"])
    effects = pd.DataFrame({f"S{i}": [float(i)] for i in range(6)}, index=["MYCN"])
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    with pytest.raises(KeyError):
        expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "NOT_A_COMPOUND")


# --- QA-driven regression tests -------------------------------------------


def test_expression_vs_aux_splits_the_two_exclusion_causes():
    """A sample can be absent for two genuinely different reasons -- it
    has no aux layer at all, or it has the layer but this particular
    feature was never assayed on it. Collapsing both into one count
    produced a self-contradicting UI ("covers 91 of 186" above "173
    excluded, no crispr data") that blamed missing-layer for samples that
    had the layer. Caught in QA on the CRISPR tab's own default value."""
    # 10 dataset samples; 6 carry the layer; of those, 2 are NaN for this feature.
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame(
        {f"S{i}": [float(i) if i < 4 else np.nan] for i in range(6)}, index=["SPARSE"]
    )
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    result = expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "SPARSE")

    assert result["n"] == 4
    assert result["n_missing_layer"] == 4          # S6..S9 have no layer
    assert result["n_missing_feature"] == 2        # S4, S5 have the layer but no value
    assert result["n_excluded"] == 6
    assert "no crispr gene effect data" in result["exclusion_reason"]
    assert "not assayed for 'SPARSE'" in result["exclusion_reason"]


def test_expression_vs_aux_too_few_samples_error_distinguishes_the_causes():
    """The 'not enough samples' error should say how many carry the layer
    versus how many were assayed for this feature -- otherwise a user
    picking a sparsely-screened compound can't tell whether the dataset
    or their compound choice is the problem."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame(
        {f"S{i}": [1.0 if i == 0 else np.nan] for i in range(8)}, index=["SPARSE"]
    )
    ds = _dataset(expression, {AuxLayer.CRISPR_GENE_EFFECT: _numeric_layer(effects)})

    with pytest.raises(ValueError) as exc:
        expression_vs_aux_correlation(ds, AuxLayer.CRISPR_GENE_EFFECT, "MYCN", "SPARSE")

    # Wording changed when the zero-coverage case was split out, but the
    # requirement is unchanged: both denominators must be visible, so a
    # user can tell "this dataset barely carries the layer" from "this
    # feature specifically was barely screened".
    msg = str(exc.value)
    assert "measured on only 1" in msg     # this feature, on this cohort
    assert "8 covered samples" in msg      # what the layer covers at all


def test_expression_by_mutation_status_rejects_gene_mutated_only_outside_the_cohort():
    """A gene can be present in the mutation matrix yet mutated in zero
    samples that also have expression data -- 2,701 of TARGET's 6,250
    genes are in that state. That used to return a degenerate 200
    (n_mutated=0, p_value=null) instead of explaining the problem."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(6)}, index=["MYCN"])
    # ELSEWHERE is mutated only in S9, which has no expression data.
    mutations = pd.DataFrame(
        {**{f"S{i}": [False] for i in range(6)}, "S9": [True]}, index=["ELSEWHERE"]
    )
    ds = _dataset(expression, {AuxLayer.MUTATION_STATUS: _mutation_layer(mutations)})

    with pytest.raises(ValueError, match="not mutated in any of the"):
        expression_by_mutation_status(ds, "MYCN", "ELSEWHERE")


def test_n_usable_features_counts_only_adequately_assayed_features():
    """n_features alone oversells a sparse layer: PRISM ships 6,790
    compounds but only ~22% were assayed on enough DepMap lines to be
    analysable."""
    from app.models.contract import AuxMatrix

    frame = pd.DataFrame(
        {
            "S1": [1.0, 1.0, np.nan],
            "S2": [2.0, 2.0, np.nan],
            "S3": [3.0, np.nan, np.nan],
        },
        index=["DENSE", "THIN", "EMPTY"],
    )
    layer = AuxMatrix(
        layer=AuxLayer.DRUG_SENSITIVITY, matrix=frame, value_label="x",
        value_description="x", source_note="x",
    )
    assert layer.n_features == 3
    # DENSE has 3 values, THIN has 2, EMPTY has 0 -> only DENSE clears min 3.
    assert layer.n_usable_features(["S1", "S2", "S3"], min_samples=3) == 1
    assert layer.n_usable_features(["S1", "S2", "S3"], min_samples=2) == 2
    assert layer.n_usable_features([], min_samples=3) == 0


def _aux_layer(index: list[str], layer: AuxLayer = AuxLayer.CRISPR_GENE_EFFECT) -> AuxMatrix:
    frame = pd.DataFrame(
        {"S1": [1.0] * len(index), "S2": [2.0] * len(index), "S3": [3.0] * len(index)},
        index=index,
    )
    return AuxMatrix(
        layer=layer, matrix=frame, value_label="x", value_description="x", source_note="x"
    )


def test_resolve_aux_feature_normalises_case_and_whitespace():
    """`?gene=myb` resolved while `?aux_feature=myb` 404'd, so one pane
    disagreed with itself about whether input was case-sensitive. Worse,
    the error read as a coverage claim ("not found in the layer"), which
    could talk a user out of a valid analysis."""
    from app.api.main import _resolve_aux_feature

    layer = _aux_layer(["MYB", "RPL13A"])
    for probe in ("MYB", "myb", "  myb  ", "MyB"):
        assert _resolve_aux_feature(layer, probe) == "MYB"


def test_resolve_aux_feature_preserves_exact_match_over_folding():
    """Drug compounds are keyed by Broad id, not symbol, so the exact
    identifier must win before any case folding is attempted."""
    from app.api.main import _resolve_aux_feature

    compound = "BRD:BRD-A00047421-001-01-7"
    layer = _aux_layer([compound], layer=AuxLayer.DRUG_SENSITIVITY)
    assert _resolve_aux_feature(layer, compound) == compound
    assert _resolve_aux_feature(layer, compound.lower()) == compound


def test_resolve_aux_feature_refuses_ambiguous_fold():
    """If a layer ever carries two labels differing only by case, guessing
    silently would return an arbitrary one -- say so instead."""
    from fastapi import HTTPException

    from app.api.main import _resolve_aux_feature

    layer = _aux_layer(["ABC", "abc"])
    assert _resolve_aux_feature(layer, "ABC") == "ABC"  # exact still wins
    with pytest.raises(HTTPException) as exc:
        _resolve_aux_feature(layer, "AbC")
    assert exc.value.status_code == 422


def test_resolve_aux_feature_unknown_raises_404():
    from fastapi import HTTPException

    from app.api.main import _resolve_aux_feature

    layer = _aux_layer(["MYB"])
    with pytest.raises(HTTPException) as exc:
        _resolve_aux_feature(layer, "ZZZNOPE")
    assert exc.value.status_code == 404


def test_zero_coverage_feature_says_it_was_never_assayed():
    """5,272 of PRISM's 6,790 compounds (77.6%) have zero measurements on
    DepMap's covered lymphoid lines. Reporting that as "need at least 3 to
    correlate" frames it as a threshold the user narrowly missed, when the
    truth is the compound was never run on this cohort and no gene choice
    will change it."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame({f"S{i}": [np.nan] for i in range(6)}, index=["NEVER_RUN"])
    ds = _dataset(expression, {AuxLayer.DRUG_SENSITIVITY: _numeric_layer(effects)})

    with pytest.raises(ValueError) as exc:
        expression_vs_aux_correlation(ds, AuxLayer.DRUG_SENSITIVITY, "MYCN", "NEVER_RUN")

    msg = str(exc.value)
    assert "never assayed on this cohort" in msg
    assert "no correlation is possible for any gene" in msg
    # Must NOT frame a never-measured compound as a near-miss on a threshold.
    assert "need at least" not in msg


def test_partial_coverage_feature_reports_measured_count():
    """Distinct from the zero case: the compound WAS run here, just not on
    enough lines. The user should see how many, so they can tell a sparse
    compound from an unusable one."""
    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    effects = pd.DataFrame(
        {f"S{i}": [1.0 if i < 2 else np.nan] for i in range(6)}, index=["SPARSE"]
    )
    ds = _dataset(expression, {AuxLayer.DRUG_SENSITIVITY: _numeric_layer(effects)})

    with pytest.raises(ValueError) as exc:
        expression_vs_aux_correlation(ds, AuxLayer.DRUG_SENSITIVITY, "MYCN", "SPARSE")

    msg = str(exc.value)
    assert "measured on only 2" in msg
    assert "need at least 3 paired samples" in msg
    assert "never assayed" not in msg


def test_min_aux_samples_is_shared_by_correlation_and_usable_count():
    """The picker offers features that clear `n_usable_features`, and the
    correlation then accepts or rejects them. If those two disagreed, the
    UI would list features that immediately fail -- worse than no list."""
    from app.models.contract import MIN_AUX_SAMPLES

    expression = pd.DataFrame({f"S{i}": [float(i)] for i in range(10)}, index=["MYCN"])
    # Exactly at the threshold: must be counted usable AND must correlate.
    effects = pd.DataFrame(
        {f"S{i}": [float(i) if i < MIN_AUX_SAMPLES else np.nan] for i in range(6)},
        index=["EXACTLY_MIN"],
    )
    layer = _numeric_layer(effects)
    ds = _dataset(expression, {AuxLayer.DRUG_SENSITIVITY: layer})
    shared = [c for c in ds.matrix.columns if c in layer.matrix.columns]

    assert layer.n_usable_features(shared) == 1
    result = expression_vs_aux_correlation(ds, AuxLayer.DRUG_SENSITIVITY, "MYCN", "EXACTLY_MIN")
    assert result["n"] == MIN_AUX_SAMPLES
