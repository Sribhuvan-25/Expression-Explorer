"""
Shared helpers for GEO microarray loaders (series matrix + platform annot).

Extracted from `gds4299.py` so a second GEO loader (`gse66586.py`) reuses
the exact same probe-mapping and probe-collapse rules instead of a copy
that could drift. Everything here is platform-agnostic; anything that is
specific to one platform (its legacy-symbol alias map) is passed in by the
calling loader. The extraction was verified bit-identical on GDS4299's
matrix, samples and features before and after.
"""
from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
from typing import Mapping

import httpx
import pandas as pd

# GEO uses these as "no mapping", not as gene symbols.
NULL_SYMBOLS = frozenset({"---", "NA", ""})


def download(url: str, dest: Path) -> Path:
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    with httpx.stream("GET", url, follow_redirects=True, timeout=120.0) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as fh:
            for chunk in r.iter_bytes():
                fh.write(chunk)
        tmp.rename(dest)
    return dest


def read_geo_table(path: Path, begin: str, end: str) -> pd.DataFrame:
    """Pull the delimited table out of a GEO file's !..._begin/!..._end block."""
    with gzip.open(path, "rt", errors="replace") as fh:
        text = fh.read()
    start, stop = text.index(begin), text.index(end)
    body = text[start:stop].split("\n", 1)[1]
    return pd.read_csv(io.StringIO(body), sep="\t", low_memory=False)


def read_series_matrix_expression(path: Path) -> pd.DataFrame:
    """The probe x sample table of a series matrix, ids as clean strings."""
    expr = read_geo_table(path, "!series_matrix_table_begin", "!series_matrix_table_end")
    expr = expr.set_index(expr.columns[0])
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = [str(c).strip('"') for c in expr.columns]
    return expr


def read_series_matrix_header(path: Path) -> dict[str, list[list[str]]]:
    """The `!Sample_*` header rows of a series matrix, one list per row.

    A key can repeat (`!Sample_characteristics_ch1` appears once per
    characteristic), so each key maps to a list of rows, in file order.
    Each row holds one value per sample, in the same column order as the
    expression table.
    """
    rows: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            fields = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
            rows.setdefault(fields[0], []).append([f.strip() for f in fields[1:]])
    return rows


def probe_to_symbol(
    annot_path: Path,
    aliases: Mapping[str, str],
    null_symbols: frozenset[str] = NULL_SYMBOLS,
) -> pd.Series:
    """probe_id -> gene symbol, for probes that actually map to a gene.

    A probe may list several symbols ('A /// B'); GEO also uses '---' and
    'NA' to mean "no mapping". Both were handled in the prior verified
    analysis of GDS4299 and are handled the same way for every GEO loader:
    take the first symbol, drop the sentinels. `aliases` maps legacy
    symbols the platform annotation still carries to their current HGNC
    symbol; pass `{}` to get the annotation's own (raw) first symbols.
    """
    ann = read_geo_table(annot_path, "!platform_table_begin", "!platform_table_end")
    symbol_col = "Gene symbol" if "Gene symbol" in ann.columns else "Gene Symbol"
    ann = ann[["ID", symbol_col]].dropna()
    ann["ID"] = ann["ID"].astype(str).str.strip()
    symbols = (
        ann[symbol_col]
        .astype(str)
        .str.split(r"[;,/\s]+", regex=True)
        .str[0]
        .str.strip()
    )
    symbols = symbols.replace(dict(aliases))
    keep = ~symbols.isin(null_symbols)
    return pd.Series(symbols[keep].to_numpy(), index=ann["ID"][keep].to_numpy())


def collapse_probes(expr: pd.DataFrame, probe_symbol: pd.Series) -> pd.DataFrame:
    """One row per gene, choosing the probe with the highest mean expression.

    The selection rule is deliberately independent of any sample grouping.
    The earlier exploratory analysis of GDS4299 picked, per gene, the probe
    with the lowest ETP-vs-non-ETP p-value -- that is circular for inference
    (choose the probe that best separates the groups, then test whether
    the groups separate), so it is not used here. See METHODS.md 7.1.

    Highest mean also does the obvious right thing on GPL13158: MYCN has
    10 probes, one at mean 5.36 and nine sitting at ~2.7 background.

    Probes whose legacy symbol was aliased onto a current one compete here
    on equal terms with any probe already annotated under the current
    symbol -- a legacy/current collision is resolved by the same rule.
    """
    shared = expr.index.intersection(probe_symbol.index)
    expr = expr.loc[shared]
    genes = probe_symbol.loc[shared]

    order = expr.mean(axis=1).sort_values(ascending=False)
    # First occurrence per gene in mean-descending order == highest mean.
    best_probe = genes.loc[order.index].reset_index().drop_duplicates(subset=0, keep="first")
    best_probe.columns = ["probe_id", "symbol"]

    collapsed = expr.loc[best_probe["probe_id"]]
    collapsed.index = pd.Index(best_probe["symbol"], name="feature_id")
    return collapsed.sort_index()
