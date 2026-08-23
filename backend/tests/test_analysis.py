import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from app.analysis.compare import expression_by_group, kruskal_wallis, pairwise_tests
from app.analysis.signature import auc_signature_score, log2_mean_signature_score
from app.analysis.survival import binarize_by_median, build_survival_frame, kaplan_meier_curves


def _synthetic_matrix(n_genes=200, n_samples=40, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)]
    data = rng.lognormal(mean=2, sigma=1.5, size=(n_genes, n_samples))
    return pd.DataFrame(data, index=genes, columns=samples)


def test_signature_score_ranges():
    matrix = _synthetic_matrix()
    scores = auc_signature_score(matrix, ["GENE1", "GENE2", "GENE3"])
    assert (scores >= 0).all() and (scores <= 1).all()
    assert not scores.isna().any()


def test_signature_score_missing_genes_returns_nan():
    matrix = _synthetic_matrix()
    scores = auc_signature_score(matrix, ["NOT_A_GENE"])
    assert scores.isna().all()


def test_log2_mean_matches_manual_calc():
    matrix = pd.DataFrame({"S1": [3.0, 7.0], "S2": [1.0, 1.0]}, index=["G1", "G2"])
    scores = log2_mean_signature_score(matrix, ["G1", "G2"])
    expected_s1 = np.mean([np.log2(4), np.log2(8)])
    assert abs(scores["S1"] - expected_s1) < 1e-9


def test_compare_groups_and_stats():
    matrix = _synthetic_matrix()
    samples = pd.DataFrame({
        "sample_id": matrix.columns,
        "group_columns": [{"status": "A" if i % 2 == 0 else "B"} for i in range(matrix.shape[1])],
    })
    df = expression_by_group(matrix, samples, "GENE0", "status")
    assert set(df["group"].unique()) == {"A", "B"}
    tests = pairwise_tests(df)
    assert "p_value" in tests.columns and len(tests) == 1
    kw = kruskal_wallis(df)
    assert kw["p_value"] is not None


def test_survival_pipeline_runs():
    n = 60
    rng = np.random.default_rng(1)
    samples = pd.DataFrame({
        "sample_id": [f"P{i}" for i in range(n)],
        "vital_status": rng.choice(["Alive", "Dead"], n),
        "days_to_death": [rng.integers(100, 2000) if i % 3 == 0 else np.nan for i in range(n)],
        "days_to_last_follow_up": rng.integers(100, 2000, n).astype(float),
    })
    scores = pd.Series(rng.normal(size=n), index=samples["sample_id"])
    surv_df = build_survival_frame(samples, scores)
    assert len(surv_df) > 0
    group = binarize_by_median(surv_df)
    result = kaplan_meier_curves(surv_df, group)
    assert "curves" in result and len(result["curves"]) == 2
    assert "logrank_p_value" in result
