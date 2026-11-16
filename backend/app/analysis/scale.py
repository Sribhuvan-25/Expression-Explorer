"""
The scale correlation, co-expression and PCA run on (METHODS.md 7.6a).

Linear units (TPM) are analysed as log2(x + 1); data already on a log2
scale is used as is. Expression is roughly log-normal, so Pearson r and
PCA on linear TPM are driven by a few very high values: TARGET DTX1 vs
NOTCH3 gave r = 0.113, p = 0.014 on TPM against r = 0.037, p = 0.42 on
log2(TPM + 1) (Spearman 0.064, p = 0.17). Rank-based coefficients
(Spearman, Kendall) do not change. Same transform as the log2-mean score
(METHODS.md 2.4).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

LINEAR_UNITS = {"tpm", "fpkm", "cpm", "counts"}


def unit_name(unit) -> str:
    return str(getattr(unit, "value", unit))


def on_analysis_scale(values: pd.DataFrame | pd.Series, unit) -> tuple[pd.DataFrame | pd.Series, str]:
    """(values on the analysis scale, the unit they are now in)."""
    name = unit_name(unit)
    if name in LINEAR_UNITS:
        return np.log2(values + 1), f"log2_{name}"
    return values, name
