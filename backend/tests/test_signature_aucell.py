"""
The AUC signature score is AUCell's normalised AUC (METHODS.md 2.1-2.3).

It used to integrate the recovery curve with a trapezoid (interpolating
between hits, which over-counts area) and divide by aucThreshold x n
rather than AUCell's maxAUC. Scores ran up to 0.10 high and moved samples
across the median split that the survival pane uses (16 of 186 DepMap
lines and 8 of 469 TARGET patients for one five-gene set).
"""
import itertools

import numpy as np
import pandas as pd
import pytest

from app.analysis.signature import auc_signature_score


def aucell_auc(ranks, threshold, n_genes):
    """AUCell's `.auc` with normAUC, transcribed from R/02_calcAUC.R:
    x <- sort(x[x < aucThreshold]); sum(diff(c(x, aucThreshold)) * seq_along(x)) / maxAUC."""
    x = np.sort([r for r in ranks if r < threshold])
    auc = float(np.sum(np.diff(np.append(x, threshold)) * np.arange(1, len(x) + 1)))
    x_th = np.arange(1, n_genes + 1)
    x_th = x_th[x_th < threshold]
    max_auc = float(np.sum(np.diff(np.append(x_th, threshold)) * np.arange(1, len(x_th) + 1)))
    return auc / max_auc


def test_matches_aucell_on_untied_ranks():
    # 8 genes, top_fraction 0.5 -> aucThreshold 4. S1 and S2 rank 1 and 3.
    matrix = pd.DataFrame({"s": [80, 70, 60, 50, 40, 30, 20, 10]}, index=["S1", "A", "S2", "B", "C", "D", "E", "F"])
    score = auc_signature_score(matrix, ["S1", "S2"], top_fraction=0.5)["s"]
    # (4-1) + (4-3) = 4 over the best case (4-1) + (4-2) = 5.
    assert score == pytest.approx(0.8)
    assert score == pytest.approx(aucell_auc([1, 3], 4, 2))


def test_signature_at_the_top_scores_one_and_duplicates_count_once():
    matrix = pd.DataFrame({"s": [9.0, 8.0, 7.0, 1.0, 0.5, 0.2, 0.1, 0.0]}, index=list("XYZabcde"))
    assert auc_signature_score(matrix, ["X", "Y", "Y", "Z"], top_fraction=0.5)["s"] == pytest.approx(1.0)


def test_ties_give_aucells_score_averaged_over_its_random_tie_breaks():
    # AUCell ranks ties at random. Enumerate every ordering of each tie
    # block, score each with AUCell's formula, and average.
    rng = np.random.default_rng(3)
    values = rng.integers(0, 6, size=(14, 5)).astype(float)  # heavy ties, some straddling the cut
    genes = [f"G{i}" for i in range(14)]
    matrix = pd.DataFrame(values, index=genes, columns=[f"s{j}" for j in range(5)])
    signature = ["G0", "G3", "G7", "G11"]
    threshold = round(0.5 * 14)
    got = auc_signature_score(matrix, signature, top_fraction=0.5)

    for col in matrix.columns:
        v = matrix[col]
        blocks = [list(v.index[v == level]) for level in sorted(v.unique(), reverse=True)]
        scores = []
        for orders in itertools.product(*(itertools.permutations(b) for b in blocks)):
            rank = {g: i + 1 for i, g in enumerate(itertools.chain.from_iterable(orders))}
            scores.append(aucell_auc([rank[g] for g in signature], threshold, len(signature)))
        assert got[col] == pytest.approx(np.mean(scores), abs=1e-12), col


def test_matches_r_aucell():
    """Pinned to R AUCell 1.20.2 (AUCell_buildRankings + AUCell_calcAUC,
    aucMaxRank = 50, normAUC = TRUE) on this matrix, which has no ties, so
    AUCell's random tie-breaking does not enter."""
    rng = np.random.default_rng(11)
    matrix = pd.DataFrame(
        rng.lognormal(2, 1.5, size=(200, 4)),
        index=[f"G{i}" for i in range(200)],
        columns=[f"s{j}" for j in range(4)],
    )
    got = auc_signature_score(matrix, [f"G{i}" for i in range(0, 200, 17)], top_fraction=0.25)
    assert got.to_numpy() == pytest.approx(R_AUCELL_PINNED, abs=1e-12)


R_AUCELL_PINNED = [0.22413793103448276, 0.11494252873563218, 0.31417624521072796, 0.23563218390804597]
