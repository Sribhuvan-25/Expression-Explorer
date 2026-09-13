import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from app.ingest import depmap_crispr, depmap_prism, target_mutations
from app.models.contract import AuxLayer


# --- TARGET somatic mutations ---------------------------------------------


def test_non_silent_set_excludes_silent_and_utr_variants():
    """The whole point of the allowlist: a synonymous or untranslated-region
    variant doesn't change the protein, so counting it as "this gene is
    mutated" inflates every mutation rate. Verified against a real MAF where
    7 of 21 calls in one aliquot were Silent and 2 were 3'UTR."""
    assert "Silent" not in target_mutations.NON_SILENT
    assert "3'UTR" not in target_mutations.NON_SILENT
    assert "5'UTR" not in target_mutations.NON_SILENT
    assert "Intron" not in target_mutations.NON_SILENT


def test_non_silent_set_includes_the_protein_changing_classes():
    for cls in (
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Splice_Site",
        "In_Frame_Del",
    ):
        assert cls in target_mutations.NON_SILENT


def test_mutation_matrix_is_boolean_and_dense(monkeypatch, tmp_path):
    """Two samples with partially overlapping mutated genes must produce a
    dense boolean frame -- every (gene, sample) cell present and False where
    that sample has no call in that gene, not NaN."""
    monkeypatch.setattr(target_mutations, "CACHE_DIR", tmp_path)

    fake_files = [
        {"file_id": "f1", "cases": [{"submitter_id": "S1"}]},
        {"file_id": "f2", "cases": [{"submitter_id": "S2"}]},
    ]
    calls = {
        "f1": pd.DataFrame(
            {"Hugo_Symbol": ["NOTCH1", "TP53", "QUIET"],
             "Variant_Classification": ["Missense_Mutation", "Nonsense_Mutation", "Silent"]}
        ),
        "f2": pd.DataFrame(
            {"Hugo_Symbol": ["TP53"], "Variant_Classification": ["Frame_Shift_Del"]}
        ),
    }

    monkeypatch.setattr(target_mutations, "_list_maf_files", lambda client: fake_files)
    monkeypatch.setattr(target_mutations, "_download_maf", lambda client, fid: calls[fid])

    layer = target_mutations.load_mutations(use_cache=False)

    assert layer.layer is AuxLayer.MUTATION_STATUS
    assert layer.matrix.dtypes.unique().tolist() == [bool]
    # QUIET was Silent-only -> must not appear as a mutated gene at all.
    assert "QUIET" not in layer.matrix.index
    assert layer.matrix.loc["NOTCH1", "S1"] is True or bool(layer.matrix.loc["NOTCH1", "S1"])
    assert not bool(layer.matrix.loc["NOTCH1", "S2"])  # dense False, not missing
    assert bool(layer.matrix.loc["TP53", "S1"]) and bool(layer.matrix.loc["TP53", "S2"])


def test_mutation_loader_skips_a_file_that_fails_to_download(monkeypatch, tmp_path):
    """One unreachable MAF out of hundreds must not abort the whole run."""
    monkeypatch.setattr(target_mutations, "CACHE_DIR", tmp_path)
    fake_files = [
        {"file_id": "good", "cases": [{"submitter_id": "S1"}]},
        {"file_id": "bad", "cases": [{"submitter_id": "S2"}]},
    ]

    def _download(client, fid):
        if fid == "bad":
            raise RuntimeError("simulated GDC failure")
        return pd.DataFrame(
            {"Hugo_Symbol": ["NOTCH1"], "Variant_Classification": ["Missense_Mutation"]}
        )

    monkeypatch.setattr(target_mutations, "_list_maf_files", lambda client: fake_files)
    monkeypatch.setattr(target_mutations, "_download_maf", _download)

    layer = target_mutations.load_mutations(use_cache=False)
    assert list(layer.matrix.columns) == ["S1"]


def test_mutation_loader_returns_none_when_nothing_downloaded(monkeypatch, tmp_path):
    monkeypatch.setattr(target_mutations, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(target_mutations, "_list_maf_files", lambda client: [])
    assert target_mutations.load_mutations(use_cache=False) is None


# --- DepMap PRISM drug sensitivity ----------------------------------------


def test_prism_release_title_pattern_is_distinct_from_depmap_bundle():
    """PRISM ships as its own Figshare release ("Repurposing Public YYQN"),
    NOT inside "DepMap YYQN Public" -- a single shared regex would silently
    never match it."""
    assert depmap_prism._RELEASE_TITLE.match("Repurposing Public 24Q2")
    assert not depmap_prism._RELEASE_TITLE.match("DepMap 24Q4 Public")


def test_prism_find_matches_by_suffix_not_exact_name():
    """Filenames carry the release version ("Repurposing_Public_24Q2_..."),
    so a next release renames every file -- matching on suffix keeps this
    working across versions instead of pinning to one quarter."""
    files = {
        "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv": "u1",
        "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv": "u2",
        "README.txt": "u3",
    }
    assert depmap_prism._find(files, depmap_prism._MATRIX_SUFFIX) == "u1"
    assert depmap_prism._find(files, depmap_prism._COMPOUND_SUFFIX) == "u2"
    assert depmap_prism._find(files, "_NotPresent.csv") is None


def test_prism_returns_none_when_no_release_found(monkeypatch, tmp_path):
    monkeypatch.setattr(depmap_prism, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(depmap_prism, "_resolve_release", lambda: None)
    assert depmap_prism.load_drug_sensitivity(use_cache=False) is None


# --- DepMap CRISPR --------------------------------------------------------


def test_crispr_returns_none_when_release_lacks_the_file(monkeypatch, tmp_path):
    """A release without CRISPRGeneEffect.csv must degrade the depmap
    dataset to expression-only, not raise."""
    monkeypatch.setattr(depmap_crispr, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(
        depmap_crispr, "_resolve_release_urls", lambda: ({"SomethingElse.csv": "u"}, "DepMap 24Q4 Public")
    )
    assert depmap_crispr.load_crispr(use_cache=False) is None


def test_crispr_layer_metadata_describes_the_score_direction():
    """A CRISPR gene-effect score near -1 means 'essential', which is the
    opposite of the intuition for an expression value -- the meaning has to
    travel with the data so a UI can't render it unlabelled."""
    layer = depmap_crispr.AuxMatrix(
        layer=AuxLayer.CRISPR_GENE_EFFECT,
        matrix=pd.DataFrame({"ACH-1": [-1.0]}, index=["RPL13A"]),
        value_label="Gene effect (Chronos)",
        value_description="more negative means the line depends more strongly on that gene",
        source_note="test",
    )
    assert "negative" in layer.value_description.lower()
