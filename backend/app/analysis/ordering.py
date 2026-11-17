"""Deterministic ordering for "top N by count" lists.

pandas' default sort is not stable, so equal counts came out in storage
order -- which differs between machines -- and so did which items made a
top-N cut. Ties are broken by the item's name, the same everywhere.
"""
from __future__ import annotations

import pandas as pd


def top_by_count(counts: pd.Series, n: int) -> pd.Series:
    """The n largest counts, equal counts ordered by index label."""
    return counts.loc[sorted(counts.index, key=lambda k: (-counts[k], str(k)))[:n]]
