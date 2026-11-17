"""
ETP-status enrichment for TARGET T-ALL patients, from the supplementary
tables of Liu et al. 2017 (Nat Genet 49:1211-1218) — the genomic
characterization paper for this exact TARGET T-ALL cohort.

Why this exists: GDC's own clinical files for TARGET-ALL-P2 do not carry
ETP/near-ETP/non-ETP classification anywhere, even though the CDE
dictionary defines those as permissible "Cell of Origin" values (checked
all 9 GDC clinical supplement files — none populate it). Wang et al. 2025
(the seed paper for this project) filters TARGET patients by this exact
status, so without it, ETP-status grouping — the central axis in most of
that paper's figures — is unavailable. This module is the fix.

Not a standalone Dataset: this returns a lookup table joined onto
TARGET-ALL-P2's `samples` frame by gdc_target.py, not registered as its
own dataset.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from app.config import settings
from app.ingest._http import fetch_with_retry, validate_excel_bytes

SUPP_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fng.3909/MediaObjects/41588_2017_BFng3909_MOESM2_ESM.xlsx"
)
CACHE_DIR = settings.cache_dir / "liu2017_etp"

# The same labels, committed: 190 cases (146 non-ETP, 25 near-ETP, 19 ETP),
# built from SUPP_URL with exactly the parsing below and verified identical
# to a fresh download (source sha256 282a3ef4bed82ecd44df9da30b9e0c6130a1f6ae0df501ae1b39faee9a821d36,
# 2026-11-17). Used when the download fails: in production Springer Nature
# answers the server's request with a 3 KB HTML page, not the workbook, so
# TARGET loaded with no etp_status at all (METHODS.md 7.10).
SNAPSHOT_PATH = Path(__file__).parent / "data" / "liu2017_target_etp_status.csv"

logger = logging.getLogger(__name__)
SHEET_NAME = "Table S1 cohort"

# Liu et al.'s values are unhyphenated ('nearETP', 'notETP'); normalize to
# the hyphenated form used elsewhere in this project's group labels and UI.
_STATUS_MAP = {
    "ETP": "ETP",
    "nearETP": "near-ETP",
    "notETP": "non-ETP",
}


def load_etp_status(use_cache: bool = True) -> pd.Series:
    """Returns a Series of normalized ETP status, indexed by GDC-style
    sample_id (e.g. 'TARGET-10-PARNJB'). Only ETP/near-ETP/non-ETP rows
    are included — 'Unevaluable' and blank (no immunophenotyping done)
    rows are dropped since they're not a usable classification."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / "etp_status.parquet"

    if use_cache and cache_path.exists():
        df = pd.read_parquet(cache_path)
    else:
        # Retries transient failures AND a 200 whose body isn't actually
        # an Excel file -- see app/ingest/_http.py. This is the sibling
        # of target_mrd.py's fetch: same single-shot-with-no-defense
        # shape, and ETP status is used more widely, so it's fixed the
        # same way rather than waiting for it to fail in production too.
        try:
            content = fetch_with_retry(SUPP_URL, timeout=120.0, validate=validate_excel_bytes)
        except Exception as exc:  # noqa: BLE001 -- any failure falls back to the verified snapshot
            logger.warning("Liu 2017 supplement unavailable (%s); using the committed ETP label snapshot", exc)
            return pd.read_csv(SNAPSHOT_PATH).set_index("sample_id")["etp_status"]
        raw = pd.read_excel(pd.io.common.BytesIO(content), sheet_name=SHEET_NAME)
        df = raw[["USI", "ETP status"]].copy()
        df["sample_id"] = "TARGET-10-" + df["USI"].astype(str)
        df["etp_status"] = df["ETP status"].map(_STATUS_MAP)
        df = df.dropna(subset=["etp_status"])[["sample_id", "etp_status"]]
        df.to_parquet(cache_path)

    return df.set_index("sample_id")["etp_status"]
