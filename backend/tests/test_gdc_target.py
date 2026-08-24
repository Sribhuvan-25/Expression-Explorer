import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import httpx
import pandas as pd
import pytest

from app.ingest import gdc_target


def _tsv_response(gene_ids, gene_names, tpms):
    lines = ["# gene-model: fake"]
    lines.append("gene_id\tgene_name\ttpm_unstranded")
    lines.append("N_unmapped\t\t0")
    for gid, name, tpm in zip(gene_ids, gene_names, tpms):
        lines.append(f"{gid}\t{name}\t{tpm}")
    return "\n".join(lines).encode()


def test_download_star_counts_retries_transport_error(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda _: None)

    calls = {"n": 0}

    class FakeClient:
        def get(self, url, follow_redirects=True):
            calls["n"] += 1
            if calls["n"] < 3:
                raise httpx.RemoteProtocolError("Server disconnected without sending a response.")
            body = _tsv_response(["ENSG1"], ["GENE1"], [4.2])
            return httpx.Response(200, content=body, request=httpx.Request("GET", url))

    df = gdc_target._download_star_counts(FakeClient(), "fake-file-id")
    assert calls["n"] == 3
    assert df.loc["ENSG1", "gene_name"] == "GENE1"
    assert df.loc["ENSG1", "tpm_unstranded"] == 4.2


def test_download_star_counts_gives_up_after_max_retries(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda _: None)

    class AlwaysFailsClient:
        def get(self, url, follow_redirects=True):
            raise httpx.RemoteProtocolError("Server disconnected without sending a response.")

    with pytest.raises(httpx.RemoteProtocolError):
        gdc_target._download_star_counts(AlwaysFailsClient(), "fake-file-id")


def test_download_star_counts_does_not_retry_http_error(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda _: None)
    calls = {"n": 0}

    class HttpErrorClient:
        def get(self, url, follow_redirects=True):
            calls["n"] += 1
            req = httpx.Request("GET", url)
            return httpx.Response(404, content=b"not found", request=req)

    with pytest.raises(httpx.HTTPStatusError):
        gdc_target._download_star_counts(HttpErrorClient(), "fake-file-id")
    assert calls["n"] == 1


def test_load_recovers_from_corrupt_checkpoint(tmp_path, monkeypatch):
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    # Simulates a checkpoint write interrupted mid-write (e.g. an OOM
    # kill) -- not valid parquet at all.
    (tmp_path / "matrix.partial.parquet").write_bytes(b"not a parquet file")

    fetch_calls = []

    def fake_list_open_rna_files(client):
        return [{"file_id": "f-a", "cases": [{"submitter_id": "SAMPLE_A"}]}]

    def fake_download_star_counts(client, file_id):
        fetch_calls.append(file_id)
        return pd.DataFrame(
            {"gene_name": ["GENE1"], "tpm_unstranded": [9.0]},
            index=["ENSG1"],
        )

    def fake_fetch_clinical(client):
        return pd.DataFrame(
            {
                "sample_id": ["SAMPLE_A"],
                "vital_status": ["Alive"],
                "days_to_death": [None],
                "days_to_last_follow_up": [200],
            }
        )

    monkeypatch.setattr(gdc_target, "list_open_rna_files", fake_list_open_rna_files)
    monkeypatch.setattr(gdc_target, "_download_star_counts", fake_download_star_counts)
    monkeypatch.setattr(gdc_target, "fetch_clinical", fake_fetch_clinical)
    monkeypatch.setattr(gdc_target, "load_etp_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "load_mrd_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "_client", lambda: _NullContextClient())

    # Must not raise -- corrupt checkpoint should be discarded, not fatal.
    dataset = gdc_target.load(use_cache=True)

    assert fetch_calls == ["f-a"]
    assert set(dataset.matrix.columns) == {"SAMPLE_A"}


def test_load_resumes_from_partial_checkpoint(tmp_path, monkeypatch):
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    partial = pd.DataFrame({"SAMPLE_A": pd.Series([1.0, 2.0], index=["ENSG1", "ENSG2"])})
    partial.to_parquet(tmp_path / "matrix.partial.parquet")
    gene_names = pd.Series(["GENE1", "GENE2"], index=["ENSG1", "ENSG2"])
    gene_names.to_frame("gene_name").to_parquet(tmp_path / "gene_names.parquet")

    fetch_calls = []

    def fake_list_open_rna_files(client):
        return [
            {"file_id": "f-a", "cases": [{"submitter_id": "SAMPLE_A"}]},
            {"file_id": "f-b", "cases": [{"submitter_id": "SAMPLE_B"}]},
        ]

    def fake_download_star_counts(client, file_id):
        fetch_calls.append(file_id)
        return pd.DataFrame(
            {"gene_name": ["GENE1", "GENE2"], "tpm_unstranded": [9.0, 9.0]},
            index=["ENSG1", "ENSG2"],
        )

    def fake_fetch_clinical(client):
        return pd.DataFrame(
            {
                "sample_id": ["SAMPLE_A", "SAMPLE_B"],
                "vital_status": ["Alive", "Dead"],
                "days_to_death": [None, 100],
                "days_to_last_follow_up": [200, None],
            }
        )

    monkeypatch.setattr(gdc_target, "list_open_rna_files", fake_list_open_rna_files)
    monkeypatch.setattr(gdc_target, "_download_star_counts", fake_download_star_counts)
    monkeypatch.setattr(gdc_target, "fetch_clinical", fake_fetch_clinical)
    monkeypatch.setattr(gdc_target, "load_etp_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "load_mrd_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "_client", lambda: _NullContextClient())

    dataset = gdc_target.load(use_cache=True)

    # SAMPLE_A was already in the checkpoint -- only SAMPLE_B should be fetched.
    assert fetch_calls == ["f-b"]
    assert set(dataset.matrix.columns) == {"SAMPLE_A", "SAMPLE_B"}
    assert not (tmp_path / "matrix.partial.parquet").exists()


class _NullContextClient:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False
