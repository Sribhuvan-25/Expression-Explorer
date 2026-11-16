import sys
import threading
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


def _file(file_id, case_id, sample_type):
    return {
        "file_id": file_id,
        "cases": [{"submitter_id": case_id, "samples": [{"sample_type": sample_type}]}],
    }


def test_select_one_file_per_case_prefers_primary_over_recurrent():
    """Reproduces the data-accuracy bug found in production: a case with
    both a Primary and a Recurrent RNA-seq file must resolve to the
    Primary sample deterministically, not whichever file GDC happens to
    list first (56/469 real TARGET-ALL-P2 cases hit this)."""
    files = [
        _file("f-recurrent", "CASE_A", "Recurrent Blood Derived Cancer - Bone Marrow"),
        _file("f-primary", "CASE_A", "Primary Blood Derived Cancer - Bone Marrow"),
    ]
    selected = gdc_target.select_one_file_per_case(files)
    assert len(selected) == 1
    assert selected[0]["file_id"] == "f-primary"

    # Order in the input list must not matter.
    reversed_files = list(reversed(files))
    selected2 = gdc_target.select_one_file_per_case(reversed_files)
    assert selected2[0]["file_id"] == "f-primary"


def test_select_one_file_per_case_prefers_bone_marrow_over_peripheral_blood():
    files = [
        _file("f-blood", "CASE_B", "Primary Blood Derived Cancer - Peripheral Blood"),
        _file("f-marrow", "CASE_B", "Primary Blood Derived Cancer - Bone Marrow"),
    ]
    selected = gdc_target.select_one_file_per_case(files)
    assert selected[0]["file_id"] == "f-marrow"


def test_select_one_file_per_case_falls_back_to_recurrent_when_no_primary_exists():
    """1 real TARGET-ALL-P2 case has only Recurrent files -- must still
    resolve to a file, not disappear from the cohort."""
    files = [
        _file("f-r2", "CASE_C", "Recurrent Blood Derived Cancer - Bone Marrow"),
        _file("f-r1", "CASE_C", "Recurrent Blood Derived Cancer - Peripheral Blood"),
    ]
    selected = gdc_target.select_one_file_per_case(files)
    assert len(selected) == 1
    assert selected[0]["file_id"] == "f-r2"  # bone marrow preferred among same-tier options


def test_select_one_file_per_case_breaks_true_duplicate_ties_deterministically():
    """8 real TARGET-ALL-P2 cases have two files of the identical
    sample_type -- a genuine duplicate upload, not a second timepoint.
    Must resolve to the same file on every run (lowest file_id), not an
    arbitrary one."""
    files = [
        _file("f-zzz", "CASE_D", "Primary Blood Derived Cancer - Bone Marrow"),
        _file("f-aaa", "CASE_D", "Primary Blood Derived Cancer - Bone Marrow"),
    ]
    selected = gdc_target.select_one_file_per_case(files)
    assert len(selected) == 1
    assert selected[0]["file_id"] == "f-aaa"

    # Stable regardless of input order.
    selected2 = gdc_target.select_one_file_per_case(list(reversed(files)))
    assert selected2[0]["file_id"] == "f-aaa"


def test_select_one_file_per_case_leaves_single_file_cases_untouched():
    files = [_file("f-only", "CASE_E", "Primary Blood Derived Cancer - Bone Marrow")]
    selected = gdc_target.select_one_file_per_case(files)
    assert selected == files


def test_load_uses_one_file_per_case_selection(tmp_path, monkeypatch):
    """End-to-end: load() must never assemble two columns for the same
    case. Before the fix, whichever file GDC listed first silently won
    (see the removed `if sample_path.exists(): continue` skip) --
    reintroducing that behaviour by returning both files unfiltered from
    list_open_rna_files and skipping select_one_file_per_case would make
    this test's `fetch_calls` assertion fail (both files get downloaded,
    not just the Primary one)."""
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    def fake_list_open_rna_files(client):
        return [
            _file("f-recurrent", "CASE_A", "Recurrent Blood Derived Cancer - Bone Marrow"),
            _file("f-primary", "CASE_A", "Primary Blood Derived Cancer - Bone Marrow"),
        ]

    fetch_calls = []

    def fake_download_star_counts(client, file_id):
        fetch_calls.append(file_id)
        tpm = 1.0 if file_id == "f-primary" else 99.0
        return pd.DataFrame(
            {"gene_name": ["GENE1"], "tpm_unstranded": [tpm]},
            index=["ENSG1"],
        )

    def fake_fetch_clinical(client):
        return pd.DataFrame(
            {
                "sample_id": ["CASE_A"],
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

    dataset = gdc_target.load(use_cache=True, with_aux=False)

    # Exactly one column for CASE_A -- no patient appears twice.
    assert list(dataset.matrix.columns) == ["CASE_A"]
    assert fetch_calls == ["f-primary"]
    assert dataset.matrix.loc["ENSG1", "CASE_A"] == 1.0


def test_load_skips_samples_already_downloaded(tmp_path, monkeypatch):
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    # Simulates a prior run that crashed after SAMPLE_A was already
    # written to disk -- each sample is its own file now, so a crash
    # mid-run only costs the samples not yet written, and resuming is
    # just "does this sample's file already exist," no separate
    # checkpoint step to go stale or get corrupted.
    samples_dir = tmp_path / "samples"
    samples_dir.mkdir(parents=True)
    pd.DataFrame({"tpm_unstranded": [1.0, 2.0]}, index=["ENSG1", "ENSG2"]).to_parquet(
        samples_dir / "SAMPLE_A.parquet"
    )

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

    dataset = gdc_target.load(use_cache=True, with_aux=False)

    # SAMPLE_A's file already existed -- only SAMPLE_B should be fetched.
    assert fetch_calls == ["f-b"]
    assert set(dataset.matrix.columns) == {"SAMPLE_A", "SAMPLE_B"}
    assert dataset.matrix.loc["ENSG1", "SAMPLE_A"] == 1.0
    assert dataset.matrix.loc["ENSG1", "SAMPLE_B"] == 9.0
    # Per-sample files and the scratch dir are cleaned up once the
    # combined matrix is written.
    assert not samples_dir.exists()


def test_load_assembles_matrix_from_downloaded_samples(tmp_path, monkeypatch):
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    def fake_list_open_rna_files(client):
        return [
            {"file_id": "f-a", "cases": [{"submitter_id": "SAMPLE_A"}]},
            {"file_id": "f-b", "cases": [{"submitter_id": "SAMPLE_B"}]},
        ]

    def fake_download_star_counts(client, file_id):
        tpm = 1.0 if file_id == "f-a" else 2.0
        return pd.DataFrame(
            {"gene_name": ["GENE1", "GENE2"], "tpm_unstranded": [tpm, tpm]},
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

    dataset = gdc_target.load(use_cache=True, with_aux=False)

    assert set(dataset.matrix.columns) == {"SAMPLE_A", "SAMPLE_B"}
    assert dataset.matrix.loc["ENSG1", "SAMPLE_A"] == 1.0
    assert dataset.matrix.loc["ENSG1", "SAMPLE_B"] == 2.0
    assert (tmp_path / "matrix.parquet").exists()
    assert not (tmp_path / "samples").exists()

    # A second call with the assembled matrix.parquet already present
    # should read the cache, not touch the network layer at all.
    monkeypatch.setattr(
        gdc_target,
        "list_open_rna_files",
        lambda client: (_ for _ in ()).throw(AssertionError("should not re-list files")),
    )
    dataset2 = gdc_target.load(use_cache=True, with_aux=False)
    assert set(dataset2.matrix.columns) == {"SAMPLE_A", "SAMPLE_B"}


def test_load_survives_concurrent_calls(tmp_path, monkeypatch):
    """Reproduces the production crash: two overlapping requests for the
    same uncached dataset both entered the download-and-assemble section
    at once, and one's cleanup unlinking a sample file the other had
    already deleted raised FileNotFoundError. Runs load() from two
    threads simultaneously and asserts both succeed with a consistent
    result instead of one of them crashing."""
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    def fake_list_open_rna_files(client):
        # A small artificial delay per file widens the race window so
        # two threads are actually likely to overlap inside the
        # download loop, not just at entry.
        return [
            {"file_id": f"f-{i}", "cases": [{"submitter_id": f"SAMPLE_{i}"}]} for i in range(5)
        ]

    def fake_download_star_counts(client, file_id):
        time.sleep(0.01)
        return pd.DataFrame(
            {"gene_name": ["GENE1", "GENE2"], "tpm_unstranded": [1.0, 2.0]},
            index=["ENSG1", "ENSG2"],
        )

    def fake_fetch_clinical(client):
        return pd.DataFrame(
            {
                "sample_id": [f"SAMPLE_{i}" for i in range(5)],
                "vital_status": ["Alive"] * 5,
                "days_to_death": [None] * 5,
                "days_to_last_follow_up": [200] * 5,
            }
        )

    monkeypatch.setattr(gdc_target, "list_open_rna_files", fake_list_open_rna_files)
    monkeypatch.setattr(gdc_target, "_download_star_counts", fake_download_star_counts)
    monkeypatch.setattr(gdc_target, "fetch_clinical", fake_fetch_clinical)
    monkeypatch.setattr(gdc_target, "load_etp_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "load_mrd_status", lambda: pd.Series(dtype=str))
    monkeypatch.setattr(gdc_target, "_client", lambda: _NullContextClient())

    results: list = [None, None]
    errors: list = [None, None]

    def run(i):
        try:
            results[i] = gdc_target.load(use_cache=True, with_aux=False)
        except Exception as exc:  # noqa: BLE001
            errors[i] = exc

    t1 = threading.Thread(target=run, args=(0,))
    t2 = threading.Thread(target=run, args=(1,))
    t1.start()
    t2.start()
    t1.join(timeout=10)
    t2.join(timeout=10)

    assert errors == [None, None], f"concurrent load() raised: {errors}"
    assert results[0] is not None and results[1] is not None
    assert set(results[0].matrix.columns) == {f"SAMPLE_{i}" for i in range(5)}
    assert set(results[1].matrix.columns) == {f"SAMPLE_{i}" for i in range(5)}


class _NullContextClient:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False
