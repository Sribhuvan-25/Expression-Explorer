"""
Day-29 minimal residual disease (MRD) status for TARGET T-ALL patients.

Source: TARGET's own clinical supplement "ClinicalData_Phase_II_Validation"
file, hosted on GDC under TARGET-ALL-P2's Clinical Supplement file set —
NOT the "Phase_II_Discovery" file (that one is B-ALL only, despite the
generic name; checked and confirmed: 220 B-Precursor / 10 B-Precursor-other,
zero T-ALL rows). Phase_II_Validation has 766 rows across B-ALL and T-ALL;
filtering to Cell of Origin == 'T Cell ALL' gives exactly 265 patients — the
same cohort Liu et al. 2017 profiles — with MRD Day 29 populated for all of
them.

This is the same source and threshold Wang et al. 2025 describes in its
methods: "MRD values for all patients were available for day 29... Broad:
MRD_neg <=0.01, MRD_pos >0.01."
"""
from __future__ import annotations

import httpx
import pandas as pd

from app.config import settings

GDC_API = "https://api.gdc.cancer.gov"
# TARGET_ALL_ClinicalData_Phase_II_Validation_20230727.xlsx
CLINICAL_FILE_ID = "c6a63e40-f4e5-4002-9b87-3e6f3b3dfab4"
CACHE_DIR = settings.cache_dir / "target_mrd"

MRD_NEG_THRESHOLD = 0.01  # <=0.01 -> MRD-neg, >0.01 -> MRD-pos, per paper methods


def load_mrd_status(use_cache: bool = True) -> pd.Series:
    """Returns Series of {'MRD-neg', 'MRD-pos'}, indexed by GDC-style
    sample_id (e.g. 'TARGET-10-PARASZ'). T-ALL patients only — this file
    also has B-ALL rows, which are dropped."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / "mrd_status.parquet"

    if use_cache and cache_path.exists():
        df = pd.read_parquet(cache_path)
    else:
        resp = httpx.get(
            f"{GDC_API}/data/{CLINICAL_FILE_ID}", follow_redirects=True, timeout=60.0
        )
        resp.raise_for_status()
        raw = pd.read_excel(pd.io.common.BytesIO(resp.content))
        tall = raw[raw["Cell of Origin"] == "T Cell ALL"].copy()
        df = pd.DataFrame(
            {
                "sample_id": tall["TARGET USI"],
                "mrd_day29": tall["MRD Day 29"],
            }
        ).dropna(subset=["mrd_day29"])
        df["mrd_status"] = df["mrd_day29"].apply(
            lambda v: "MRD-neg" if v <= MRD_NEG_THRESHOLD else "MRD-pos"
        )
        df.to_parquet(cache_path)

    return df.set_index("sample_id")["mrd_status"]
