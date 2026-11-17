"""top_by_count backs the CRISPR/drug feature picker (/aux-features) and
the most-mutated-genes list. Equal counts used to come out in storage
order, so which items made the cut differed between machines."""
import pandas as pd

from app.analysis.ordering import top_by_count


def test_ties_at_the_cut_are_broken_by_name_whatever_the_storage_order():
    counts = {"TOP": 9, "ZETA": 4, "ALPHA": 4, "MID": 4, "LOW": 1}
    for order in (list(counts), list(reversed(counts))):
        out = top_by_count(pd.Series({k: counts[k] for k in order}), 3)
        assert list(out.index) == ["TOP", "ALPHA", "MID"]
        assert list(out) == [9, 4, 4]
