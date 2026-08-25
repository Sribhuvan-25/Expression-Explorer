import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from app.ingest.gds4299 import _collapse_probes, _probe_to_symbol, _SUBTYPE_MAP


def test_collapse_probes_picks_highest_mean_not_most_significant():
    """The probe-collapse rule must be independent of any sample grouping.

    Constructed so the two candidate rules disagree: probe A has by far
    the highest mean (real signal), while probe B is low-expressed noise
    that happens to separate the two sample halves perfectly. Selecting
    on group separation -- as the earlier exploratory analysis did -- would
    pick B, which is circular for inference. Highest-mean must pick A.
    """
    expr = pd.DataFrame(
        {
            "s1": [9.0, 1.0],
            "s2": [9.1, 1.1],
            "s3": [8.9, 5.0],
            "s4": [9.0, 5.1],
        },
        index=["probeA", "probeB"],
    )
    probe_symbol = pd.Series({"probeA": "GENE1", "probeB": "GENE1"})

    collapsed = _collapse_probes(expr, probe_symbol)

    assert list(collapsed.index) == ["GENE1"]
    # probeA's values, not probeB's.
    assert collapsed.loc["GENE1", "s1"] == 9.0
    assert collapsed.loc["GENE1", "s3"] == 8.9


def test_collapse_probes_yields_one_row_per_gene():
    expr = pd.DataFrame(
        {"s1": [5.0, 3.0, 7.0], "s2": [5.2, 3.1, 7.1]},
        index=["p1", "p2", "p3"],
    )
    probe_symbol = pd.Series({"p1": "A", "p2": "A", "p3": "B"})

    collapsed = _collapse_probes(expr, probe_symbol)

    assert sorted(collapsed.index) == ["A", "B"]
    assert not collapsed.index.duplicated().any()
    # Gene A keeps p1 (mean 5.1), not p2 (mean 3.05).
    assert collapsed.loc["A", "s1"] == 5.0


def test_collapse_probes_ignores_probes_with_no_gene_mapping():
    expr = pd.DataFrame({"s1": [5.0, 9.0]}, index=["mapped", "unmapped"])
    probe_symbol = pd.Series({"mapped": "GENE1"})

    collapsed = _collapse_probes(expr, probe_symbol)

    # The unmapped probe has a higher mean but isn't a gene, so it's gone.
    assert list(collapsed.index) == ["GENE1"]
    assert collapsed.loc["GENE1", "s1"] == 5.0


def test_subtype_map_normalizes_to_project_labels():
    """Labels must match the hyphenated form the rest of the project emits
    (TARGET's enrichment produces 'ETP' / 'near-ETP' / 'non-ETP'), or the
    same biological group appears as two different groups across datasets."""
    assert _SUBTYPE_MAP["nonETP"] == "non-ETP"
    assert _SUBTYPE_MAP["ETP"] == "ETP"


@pytest.mark.skipif(
    not (Path(__file__).resolve().parents[2] / "data/cache/gds4299/GPL13158.annot.gz").exists(),
    reason="GPL13158 annotation not downloaded",
)
def test_probe_symbol_mapping_handles_real_annotation():
    """Guards the sentinel/alias handling against the real GEO file: '---'
    and 'NA' are 'no mapping', not gene names, and legacy symbols need
    remapping or a current-symbol lookup silently finds nothing."""
    annot = Path(__file__).resolve().parents[2] / "data/cache/gds4299/GPL13158.annot.gz"
    mapping = _probe_to_symbol(annot)

    assert len(mapping) > 40000
    assert not mapping.isin({"---", "NA", ""}).any()
    # KIAA1524 is the platform's legacy symbol for CIP2A; it must surface
    # under the current name or a CIP2A query returns "gene not found".
    assert (mapping == "CIP2A").any()
    assert not (mapping == "KIAA1524").any()
    assert (mapping == "MYCN").any()
