"""
Sample ranking — sort every sample in a dataset by a gene's expression or
a signature's score, with sample metadata attached. This is a distinct
operation from group comparison: there's no statistical test here, just
an ordered list a user can scan to pick candidates.

Motivated directly by the brief behind this project (PDF p.3): "Check
MYCN and PP2A signature genes in this Cell Line RNA-Seq data! Goal is to
purchase/get some Leukemia cell lines that have high/Low MYCN
expression." That's exactly this: rank DepMap cell lines by a gene or
signature, read off the top and bottom of the list.
"""
from __future__ import annotations

import pandas as pd

from app.analysis.signature import auc_signature_score, log2_mean_signature_score


def rank_by_gene(matrix: pd.DataFrame, samples: pd.DataFrame, gene_id: str) -> pd.DataFrame:
    """One row per sample, sorted high to low by a single gene's expression,
    with every group_columns key flattened into its own column."""
    if gene_id not in matrix.index:
        raise KeyError(f"{gene_id} not found in this dataset's matrix.")
    values = matrix.loc[gene_id]
    return _assemble_ranking(samples, values.rename("value"))


def rank_by_signature(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    gene_ids: list[str],
    method: str = "auc",
) -> pd.DataFrame:
    """Same idea, but the ranking value is a multi-gene signature score
    rather than one gene's expression."""
    scorer = auc_signature_score if method == "auc" else log2_mean_signature_score
    scores = scorer(matrix, gene_ids)
    return _assemble_ranking(samples, scores.rename("value"))


# Frame-level bookkeeping columns that are never sample annotation, so
# they should not be echoed back as metadata on a ranking row. The
# survival fields are excluded deliberately: they are inputs to the
# Kaplan-Meier/Cox analysis, not annotation a reader scans while picking
# candidates, and surfacing them here just widens the table with numbers
# that mean nothing out of that context.
_NON_META_COLUMNS = {
    "sample_id",
    "dataset_id",
    "group_columns",
    "value",
    "days_to_death",
    "days_to_last_follow_up",
}


def _assemble_ranking(samples: pd.DataFrame, values: pd.Series) -> pd.DataFrame:
    df = samples[["sample_id"]].copy()
    df["value"] = values.reindex(df["sample_id"]).to_numpy()
    df = df.dropna(subset=["value"])

    def attach(col: str, series: pd.Series) -> None:
        # Not every sample carries every metadata key (e.g. TARGET-ALL-P2's
        # etp_status/mrd_status are only populated for the subset of
        # samples with that classification) -- pandas leaves those as
        # NaN, which Starlette's default JSONResponse rejects outright
        # (allow_nan=False), turning any ranking request that includes
        # an incompletely-classified sample into an unhandled 500. None
        # serializes to JSON null and the frontend already renders a
        # missing value as "—".
        mapped = df["sample_id"].map(series)
        df[col] = mapped.astype(object).where(mapped.notna(), None)

    # A dataset may expose a group column either as a key inside each
    # sample's group_columns dict (DepMap: lineage/subtype/cell_line_name)
    # or as a top-level column on the samples frame (TARGET: vital_status,
    # which is a declared field on SampleMetadata). expression_by_group
    # already handles both; ranking previously flattened only the dict, so
    # top-level columns were silently absent from every ranking row --
    # vital_status never reached the TARGET ranking table even though the
    # dataset advertises it as a group column.
    indexed = samples.set_index("sample_id")
    for col in samples.columns:
        if col in _NON_META_COLUMNS:
            continue
        attach(col, indexed[col])

    # Flatten every group_columns key into its own column so callers get
    # e.g. cell_line_name/lineage/subtype as plain fields, not a nested dict.
    if "group_columns" in samples.columns:
        meta = pd.json_normalize(indexed.loc[df["sample_id"], "group_columns"])
        meta.index = df["sample_id"].to_numpy()
        for col in meta.columns:
            attach(col, meta[col])

    return df.sort_values("value", ascending=False).reset_index(drop=True)
