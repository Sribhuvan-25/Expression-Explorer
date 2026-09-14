import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.services import gene_info


def test_shape_extracts_expected_fields():
    hit = {
        "symbol": "MYCN",
        "name": "MYCN proto-oncogene, bHLH transcription factor",
        "summary": "This gene is a member of the MYC family...",
        "alias": ["N-myc", "NMYC"],
        "ensembl": {"gene": "ENSG00000134323"},
        "entrezgene": "4613",
    }
    shaped = gene_info._shape(hit)
    assert shaped["symbol"] == "MYCN"
    assert shaped["ensembl_gene_id"] == "ENSG00000134323"
    assert "N-myc" in shaped["aliases"]


def test_shape_handles_missing_optional_fields():
    shaped = gene_info._shape({"symbol": "FOO"})
    assert shaped["aliases"] == []
    assert shaped["ensembl_gene_id"] is None


def test_cache_path_sanitizes_symbol():
    path = gene_info._cache_path("weird/../symbol!!")
    assert ".." not in str(path)
    assert path.parent == gene_info.CACHE_DIR


def test_cache_round_trip(tmp_path, monkeypatch):
    monkeypatch.setattr(gene_info, "CACHE_DIR", tmp_path)
    path = gene_info._cache_path("MYCN")
    path.write_text(json.dumps({"symbol": "MYCN"}))
    assert gene_info._read_cache(path) == {"symbol": "MYCN"}


def test_cache_expires_after_ttl(tmp_path, monkeypatch):
    monkeypatch.setattr(gene_info, "CACHE_DIR", tmp_path)
    path = gene_info._cache_path("MYCN")
    path.write_text(json.dumps({"symbol": "MYCN"}))
    old_time = time.time() - gene_info.CACHE_TTL_SECONDS - 1
    import os
    os.utime(path, (old_time, old_time))
    assert gene_info._read_cache(path) is None


def test_cache_treats_empty_dict_as_known_negative_result(tmp_path, monkeypatch):
    monkeypatch.setattr(gene_info, "CACHE_DIR", tmp_path)
    path = gene_info._cache_path("NOTAGENE")
    path.write_text(json.dumps({}))
    assert gene_info._read_cache(path) == {}


def test_shape_handles_list_shaped_ensembl_field():
    """mygene.info returns `ensembl` as a LIST when a symbol maps to
    several Ensembl genes. ZAP70 -- a core T-cell receptor signalling gene
    a T-ALL researcher would plausibly look up -- is one such symbol, and
    it crashed /genes/{symbol} with a 500 (AttributeError: 'list' object
    has no attribute 'get').

    Only 1 of 30 surveyed T-ALL genes has this shape, which is exactly why
    sampling missed it.
    """
    hit = {
        "symbol": "ZAP70",
        "name": "zeta chain of T cell receptor associated protein kinase 70",
        "alias": ["SRK", "STD", "TZK", "ZAP-70"],
        "ensembl": [
            {"gene": "ENSG00000115085", "protein": ["ENSP00000264972"]},
            {"gene": "ENSG00000288859", "protein": ["ENSP00000513759"]},
        ],
        "entrezgene": "7535",
    }
    shaped = gene_info._shape(hit)
    assert shaped["symbol"] == "ZAP70"
    assert shaped["ensembl_gene_id"] == "ENSG00000115085"
    assert shaped["entrez_gene_id"] == "7535"


def test_shape_tolerates_missing_or_odd_ensembl_shapes():
    """Absent, empty, and unexpected scalar shapes must degrade to None
    rather than raising -- an annotation lookup should never be able to
    500 the endpoint over a metadata field."""
    for ensembl in (None, {}, [], [{}], "ENSG00000115085", 42):
        shaped = gene_info._shape({"symbol": "X", "ensembl": ensembl})
        assert shaped["ensembl_gene_id"] is None, ensembl
