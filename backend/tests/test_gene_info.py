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
