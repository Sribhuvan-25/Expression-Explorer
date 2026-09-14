"""
Analyses over auxiliary measurement layers (CRISPR gene effect, drug
sensitivity, mutation status) attached to an existing dataset's samples.

Every function here follows the same rule the rest of this app follows for
partial data (METHODS.md 6.1): an aux layer never covers a whole dataset,
so each result reports how many samples it actually ran on versus how many
the dataset holds. A mutation-stratified comparison silently running on 60%
of TARGET would read as the full cohort otherwise -- exactly the confusion
the domain-expert review already caught once on survival.
"""
from __future__ import annotations

import pandas as pd
from scipy import stats

from app.models.contract import MIN_AUX_SAMPLES, AuxLayer, Dataset


def expression_by_mutation_status(
    ds: Dataset,
    expression_feature_id: str,
    mutated_gene: str,
) -> dict:
    """Expression of one gene, split by whether each sample carries a
    non-silent mutation in another (possibly the same) gene.

    Unlike the correlation endpoint, `mutated_gene == expression gene` is
    legitimate here and a real use case ("do NOTCH1-mutant samples express
    NOTCH1 differently?"), so it is deliberately NOT blocked.
    """
    aux = ds.require_aux(AuxLayer.MUTATION_STATUS)
    if mutated_gene not in aux.matrix.index:
        raise ValueError(
            f"'{mutated_gene}' has no recorded non-silent mutation in this cohort -- "
            "nothing to split on."
        )
    if expression_feature_id not in ds.matrix.index:
        raise KeyError(f"{expression_feature_id} not found in this dataset's matrix.")

    shared = [c for c in ds.matrix.columns if c in aux.matrix.columns]
    if len(shared) < 4:
        raise ValueError(
            f"Only {len(shared)} samples have both expression and mutation data -- "
            "too few to compare."
        )

    values = ds.matrix.loc[expression_feature_id, shared]
    mutated = aux.matrix.loc[mutated_gene, shared].astype(bool)

    # Being present in the mutation matrix isn't enough -- a gene can be
    # mutated only in samples that have no RNA-seq, leaving zero mutants
    # among the samples actually comparable. 2,701 of TARGET's 6,250
    # mutated genes are in exactly that state, and they previously
    # returned a degenerate 200 (n_mutated=0, p_value=null) rather than
    # saying why. Fail with the same message shape as the not-in-cohort
    # case, since to a caller it is the same problem: nothing to split on.
    if not bool(mutated.any()):
        raise ValueError(
            f"'{mutated_gene}' is not mutated in any of the {len(shared)} samples that have "
            "both expression and mutation data -- nothing to split on."
        )

    points = [
        {
            "sample_id": sid,
            "group": "Mutated" if bool(mutated[sid]) else "Wild-type",
            "value": float(values[sid]),
        }
        for sid in shared
        if pd.notna(values[sid])
    ]

    mut_vals = [p["value"] for p in points if p["group"] == "Mutated"]
    wt_vals = [p["value"] for p in points if p["group"] == "Wild-type"]

    test: dict = {"u_stat": None, "p_value": None}
    if len(mut_vals) >= 2 and len(wt_vals) >= 2:
        u, p = stats.mannwhitneyu(mut_vals, wt_vals, alternative="two-sided")
        test = {"u_stat": float(u), "p_value": float(p)}

    n_total = int(ds.matrix.shape[1])
    return {
        "mutated_gene": mutated_gene,
        "n_mutated": len(mut_vals),
        "n_wildtype": len(wt_vals),
        "n_dataset_total": n_total,
        "n_excluded": n_total - len(points),
        "exclusion_reason": (
            "no somatic mutation data for this sample" if n_total - len(points) > 0 else None
        ),
        "points": points,
        "test": test,
    }


def most_mutated_genes(ds: Dataset, top_n: int = 25) -> dict:
    """Genes carrying a non-silent mutation in the most samples -- the
    natural entry point for "what should I even split on?", since a user
    has no way to guess which genes are recurrently mutated in this cohort.
    """
    aux = ds.require_aux(AuxLayer.MUTATION_STATUS)
    shared = [c for c in ds.matrix.columns if c in aux.matrix.columns]
    if not shared:
        raise ValueError("No samples have both expression and mutation data.")

    counts = aux.matrix[shared].sum(axis=1).sort_values(ascending=False)
    counts = counts[counts > 0].head(top_n)
    return {
        "n_samples_with_mutation_data": len(shared),
        "n_dataset_total": int(ds.matrix.shape[1]),
        "genes": [
            {"gene": gene, "n_mutated": int(n), "fraction": float(n) / len(shared)}
            for gene, n in counts.items()
        ],
    }


def expression_vs_aux_correlation(
    ds: Dataset,
    layer: AuxLayer,
    expression_feature_id: str,
    aux_feature_id: str,
    method: str = "pearson",
) -> dict:
    """Correlate a gene's expression against an auxiliary measurement on
    the same samples -- e.g. does MYCN expression track with sensitivity to
    a given compound, or with dependency on a given gene.

    This is the payoff the drug/CRISPR layers exist for: both are keyed by
    the same ModelIDs as the expression matrix, so the join is direct.
    """
    aux = ds.require_aux(layer)
    if expression_feature_id not in ds.matrix.index:
        raise KeyError(f"{expression_feature_id} not found in this dataset's matrix.")
    if aux_feature_id not in aux.matrix.index:
        raise KeyError(f"{aux_feature_id} not found in the {layer.value} layer.")

    methods = {"pearson": stats.pearsonr, "spearman": stats.spearmanr}
    if method not in methods:
        raise ValueError(f"Unknown correlation method: {method!r}")

    shared = [c for c in ds.matrix.columns if c in aux.matrix.columns]
    paired = pd.DataFrame(
        {
            "expression": ds.matrix.loc[expression_feature_id, shared],
            "aux": aux.matrix.loc[aux_feature_id, shared],
        }
    ).dropna()

    if len(paired) < MIN_AUX_SAMPLES:
        n_measured = int(aux.matrix.loc[aux_feature_id, shared].notna().sum())
        label = layer.value.replace("_", " ")
        if n_measured == 0:
            # By far the common case for PRISM -- 5,272 of 6,790 compounds
            # (77.6%) were never run on any of this cohort's lines. Saying
            # "need at least 3 to correlate" frames that as a threshold the
            # user just missed, when the real answer is that this compound
            # was never tested here and no gene choice will change it. Say
            # so, and point at the picker rather than leaving them to guess
            # again (caught in QA).
            raise ValueError(
                f"'{aux_feature_id}' has no {label} measurements on any of this dataset's "
                f"{len(shared)} covered samples -- it was never assayed on this cohort, so no "
                f"correlation is possible for any gene. Pick another feature from the list."
            )
        raise ValueError(
            f"'{aux_feature_id}' was measured on only {n_measured} of this dataset's "
            f"{len(shared)} covered samples, and {len(paired)} of those also have expression "
            f"for this gene -- need at least {MIN_AUX_SAMPLES} paired samples to correlate."
        )

    stat, p_value = methods[method](paired["expression"], paired["aux"])
    n_total = int(ds.matrix.shape[1])

    # Two genuinely different reasons a sample isn't in this result, and
    # collapsing them into one count produced a self-contradicting UI: the
    # layer's own coverage note would say "covers 91 of 186 samples" while
    # the result note said "173 excluded (no crispr gene effect data)" --
    # blaming all 173 on a missing layer when 78 of them DO have the layer
    # and simply weren't assayed for this particular feature. RPL13A is
    # scored on 349 of 1178 DepMap lines upstream, so this is the common
    # case, not an edge case (caught in QA on the CRISPR tab's own default
    # placeholder). Report the split.
    n_missing_layer = n_total - len(shared)
    n_missing_feature = len(shared) - len(paired)
    reasons = []
    if n_missing_layer:
        reasons.append(f"{n_missing_layer} have no {layer.value.replace('_', ' ')} data")
    if n_missing_feature:
        reasons.append(f"{n_missing_feature} were not assayed for '{aux_feature_id}'")

    return {
        "layer": layer.value,
        "value_label": aux.value_label,
        "value_description": aux.value_description,
        "method": method,
        "n": len(paired),
        "n_dataset_total": n_total,
        "n_excluded": n_total - len(paired),
        "n_missing_layer": n_missing_layer,
        "n_missing_feature": n_missing_feature,
        "exclusion_reason": "; ".join(reasons) if reasons else None,
        "coefficient": float(stat),
        "p_value": float(p_value),
        "points": [
            {"sample_id": sid, "x": float(row["expression"]), "y": float(row["aux"])}
            for sid, row in paired.iterrows()
        ],
    }
