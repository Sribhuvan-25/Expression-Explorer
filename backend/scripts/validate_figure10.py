"""
Reproduces Wang et al. 2025, Figure 10A/B against real TARGET-ALL-P2 data:
averaged log2 expression of the ZMIZ1-5 and ZMIZ1-ETP-8 signatures, split
by day-29 MRD status, restricted to ETP + near-ETP patients.

Paper's reported numbers (Fig. 10 caption, page 22 of the manuscript):
  Cohort: ETP + near-ETP TARGET T-ALL patients, n = 16 MRD-neg / 26 MRD-pos
  Panel A (ZMIZ1-5): MEF2C, BCL2, MYB, MYCN, ZMIZ1 -- p < 0.05 (*)
  Panel B (ZMIZ1-ETP-8): ZMIZ1-5 genes + ETP-TF5 (MEF2C, LYL1, HHEX, LMO2,
    MYCN) merged, deduplicated -- p = 0.017
  Method: averaged LOG2 expression, two-sided Wilcoxon (Mann-Whitney) test

This script does not aim to match the paper's p-value exactly -- TARGET's
open RNA-seq (GDC STAR-Counts TPM) is a different processing pipeline than
whatever the paper used, and n may differ slightly depending on exact
patient-ID matching across three separate data sources (GDC, Liu et al.
2017, TARGET clinical). The goal is to confirm the *direction* and rough
*significance* reproduce -- if MRD-pos shows higher signature expression
than MRD-neg with a similarly small p-value, that's a real validation of
the pipeline. If it doesn't reproduce at all, that's a real finding too,
and should be reported as such, not hidden.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy import stats

from app.ingest.gdc_target import load

ZMIZ1_5 = ["MEF2C", "BCL2", "MYB", "MYCN", "ZMIZ1"]
ETP_TF5 = ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"]
ZMIZ1_ETP_8 = sorted(set(ZMIZ1_5) | set(ETP_TF5))


def log2_mean_score(matrix, gene_ids, symbol_to_id):
    present_ids = [symbol_to_id[g] for g in gene_ids if g in symbol_to_id]
    missing = [g for g in gene_ids if g not in symbol_to_id]
    if missing:
        print(f"  WARNING: genes not found in matrix: {missing}")
    log_vals = np.log2(matrix.loc[present_ids] + 1)
    return log_vals.mean(axis=0), present_ids


def run_panel(name, gene_list, matrix, symbol_to_id, cohort_samples, mrd_by_sample):
    print(f"\n=== Panel: {name} ===")
    print(f"Genes ({len(gene_list)}): {', '.join(sorted(gene_list))}")

    scores, present_ids = log2_mean_score(matrix, gene_list, symbol_to_id)
    print(f"Genes resolved in matrix: {len(present_ids)} / {len(gene_list)}")

    df = cohort_samples.copy()
    df["score"] = df["sample_id"].map(scores)
    df["mrd_status"] = df["sample_id"].map(mrd_by_sample)
    df = df.dropna(subset=["score", "mrd_status"])

    neg = df.loc[df["mrd_status"] == "MRD-neg", "score"]
    pos = df.loc[df["mrd_status"] == "MRD-pos", "score"]
    print(f"n = {len(neg)} MRD-neg / {len(pos)} MRD-pos (paper: n = 16 / 26)")

    if len(neg) < 2 or len(pos) < 2:
        print("Too few samples in one group to test -- skipping.")
        return

    stat, p = stats.mannwhitneyu(neg, pos, alternative="two-sided")
    direction = "MRD-pos HIGHER" if pos.mean() > neg.mean() else "MRD-neg HIGHER"
    print(f"MRD-neg mean log2 expr: {neg.mean():.3f} (n={len(neg)})")
    print(f"MRD-pos mean log2 expr: {pos.mean():.3f} (n={len(pos)})")
    print(f"Direction: {direction} (paper expects MRD-pos higher)")
    print(f"Mann-Whitney p = {p:.4g}")
    print(f"Paper's reported p: ~0.05 (panel A, marked *) / 0.017 (panel B)")


def main():
    print("Loading TARGET-ALL-P2 (using cache if present)...")
    ds = load(use_cache=True)
    print(f"Matrix: {ds.matrix.shape[0]} genes x {ds.matrix.shape[1]} samples")

    symbol_to_id = ds.features.drop_duplicates("symbol").set_index("symbol")["feature_id"].to_dict()

    etp_by_sample = ds.samples["group_columns"].apply(lambda d: d.get("etp_status"))
    mrd_by_sample = ds.samples.set_index("sample_id")["group_columns"].apply(lambda d: d.get("mrd_status"))

    cohort_mask = etp_by_sample.isin(["ETP", "near-ETP"])
    cohort_samples = ds.samples.loc[cohort_mask, ["sample_id"]].copy()
    print(f"\nETP + near-ETP patients with RNA-seq in this cohort: {len(cohort_samples)}")
    print("(Paper's Fig. 10A/B cohort: n=42 total, 16 MRD-neg / 26 MRD-pos)")

    run_panel("ZMIZ1-5 (Fig. 10A)", ZMIZ1_5, ds.matrix, symbol_to_id, cohort_samples, mrd_by_sample)
    run_panel("ZMIZ1-ETP-8 (Fig. 10B)", ZMIZ1_ETP_8, ds.matrix, symbol_to_id, cohort_samples, mrd_by_sample)


if __name__ == "__main__":
    main()
