"""
Tests for app/ingest/_http.py, and the grouping-column degradation it
enables in gdc_target.py.

Production incident (2026-09-21): target_mrd.py's single-shot fetch of
TARGET's clinical supplement got a 200 response with a body pandas
couldn't parse as Excel -- not a connection failure `httpx.TransportError`
would catch, not an HTTP error `raise_for_status()` would catch. There
was no retry and no validation, so it raised straight up out of
`load_mrd_status()`, which was called with no isolation from
`gdc_target.load()` -- killing the entire 469-sample target_all_p2
dataset over one supplementary grouping column.

These tests pin both halves of the fix: the fetch itself now retries a
bad response body the same way it retries a transport error, and even if
it still fails after retries, the dataset degrades to that one column
missing rather than disappearing.
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import httpx
import pandas as pd
import pytest

from app.ingest import gdc_target
from app.ingest._http import BadResponseBody, fetch_with_retry, validate_excel_bytes

_XLSX_BYTES = b"PK\x03\x04" + b"\x00" * 20  # minimal valid-looking signature
_XLS_BYTES = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" + b"\x00" * 20
_HTML_BYTES = b"<!DOCTYPE html>\n<html>not the file you wanted</html>"


def test_validate_excel_bytes_accepts_xlsx_and_xls():
    validate_excel_bytes(_XLSX_BYTES)  # does not raise
    validate_excel_bytes(_XLS_BYTES)  # does not raise


def test_validate_excel_bytes_rejects_html():
    """The exact production failure: a 200 response whose body is an
    HTML page (an error/redirect page, or a login wall) rather than the
    requested spreadsheet."""
    with pytest.raises(BadResponseBody, match="doesn't start with an Excel"):
        validate_excel_bytes(_HTML_BYTES)


def test_validate_excel_bytes_rejects_empty():
    with pytest.raises(BadResponseBody, match="empty"):
        validate_excel_bytes(b"")


def test_fetch_with_retry_retries_a_bad_response_body(monkeypatch):
    """A 200 with an unparseable body must be retried exactly like a
    transport error -- both are 'the network gave us something
    unusable', and retrying is cheap relative to failing a whole
    dataset load over one bad attempt."""
    monkeypatch.setattr(time, "sleep", lambda _: None)
    calls = {"n": 0}

    def fake_get(url, follow_redirects, timeout):
        calls["n"] += 1
        body = _HTML_BYTES if calls["n"] < 3 else _XLSX_BYTES

        class _Resp:
            content = body

            def raise_for_status(self):
                pass

        return _Resp()

    monkeypatch.setattr(httpx, "get", fake_get)
    result = fetch_with_retry("https://example.test/file.xlsx", validate=validate_excel_bytes)
    assert result == _XLSX_BYTES
    assert calls["n"] == 3


def test_fetch_with_retry_gives_up_after_max_attempts(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda _: None)

    def fake_get(url, follow_redirects, timeout):
        class _Resp:
            content = _HTML_BYTES

            def raise_for_status(self):
                pass

        return _Resp()

    monkeypatch.setattr(httpx, "get", fake_get)
    with pytest.raises(BadResponseBody):
        fetch_with_retry("https://example.test/file.xlsx", attempts=4, validate=validate_excel_bytes)


def test_fetch_with_retry_does_not_retry_a_real_http_error(monkeypatch):
    """A genuine 404/500 is not transient -- retrying it just delays a
    failure that will not change, so it must raise immediately."""
    monkeypatch.setattr(time, "sleep", lambda _: None)
    calls = {"n": 0}

    def fake_get(url, follow_redirects, timeout):
        calls["n"] += 1

        class _Resp:
            def raise_for_status(self):
                raise httpx.HTTPStatusError("404", request=None, response=None)

        return _Resp()

    monkeypatch.setattr(httpx, "get", fake_get)
    with pytest.raises(httpx.HTTPStatusError):
        fetch_with_retry("https://example.test/file.xlsx", validate=validate_excel_bytes)
    assert calls["n"] == 1


def test_fetch_with_retry_retries_transport_errors_too(monkeypatch):
    """The pre-existing case (connection drop), still covered after the
    refactor into a shared helper."""
    monkeypatch.setattr(time, "sleep", lambda _: None)
    calls = {"n": 0}

    def fake_get(url, follow_redirects, timeout):
        calls["n"] += 1
        if calls["n"] < 2:
            raise httpx.ConnectError("boom")

        class _Resp:
            content = _XLSX_BYTES

            def raise_for_status(self):
                pass

        return _Resp()

    monkeypatch.setattr(httpx, "get", fake_get)
    result = fetch_with_retry("https://example.test/file.xlsx", validate=validate_excel_bytes)
    assert result == _XLSX_BYTES
    assert calls["n"] == 2


class _NullContextClient:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _full_mock_load(monkeypatch, tmp_path, *, etp, mrd):
    """Mirrors test_gdc_target.py's test_load_assembles_matrix_from_downloaded_samples
    setup -- every network-touching function stubbed, only etp_status/
    mrd_status left variable, so these tests exercise gdc_target.load()'s
    own try/except around those two calls without hitting the network at
    all (both real endpoints have been observed slow/flaky from this
    exact test run, which is a second, independent reason a real fetch
    doesn't belong in this suite)."""
    monkeypatch.setattr(gdc_target, "CACHE_DIR", tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        gdc_target,
        "list_open_rna_files",
        lambda client: [
            {"file_id": "f-a", "cases": [{"submitter_id": "SAMPLE_A"}]},
            {"file_id": "f-b", "cases": [{"submitter_id": "SAMPLE_B"}]},
        ],
    )
    monkeypatch.setattr(
        gdc_target,
        "_download_star_counts",
        lambda client, file_id: pd.DataFrame(
            {"gene_name": ["GENE1", "GENE2"], "tpm_unstranded": [1.0, 1.0]},
            index=["ENSG1", "ENSG2"],
        ),
    )
    monkeypatch.setattr(
        gdc_target,
        "fetch_clinical",
        lambda client: pd.DataFrame(
            {
                "sample_id": ["SAMPLE_A", "SAMPLE_B"],
                "vital_status": ["Alive", "Dead"],
                "days_to_death": [None, 100],
                "days_to_last_follow_up": [200, None],
            }
        ),
    )
    monkeypatch.setattr(gdc_target, "load_etp_status", etp)
    monkeypatch.setattr(gdc_target, "load_mrd_status", mrd)
    monkeypatch.setattr(gdc_target, "_client", lambda: _NullContextClient())


def test_target_all_p2_loads_when_both_grouping_sources_fail(tmp_path, monkeypatch):
    """The actual regression: with BOTH third-party grouping sources
    (Springer's ETP table, GDC's MRD clinical supplement) completely
    unreachable, target_all_p2 must still load with its full expression
    matrix -- only etp_status/mrd_status are absent, not the dataset."""
    def _boom_etp(use_cache=True):
        raise RuntimeError("springer is down")

    def _boom_mrd(use_cache=True):
        raise RuntimeError("gdc clinical supplement is down")

    _full_mock_load(monkeypatch, tmp_path, etp=_boom_etp, mrd=_boom_mrd)

    ds = gdc_target.load(use_cache=True, with_aux=False)

    assert set(ds.matrix.columns) == {"SAMPLE_A", "SAMPLE_B"}
    assert ds.matrix.shape[1] == ds.source.n_samples
    group_columns = ds.samples["group_columns"]
    assert not group_columns.apply(lambda d: "etp_status" in d).any()
    assert not group_columns.apply(lambda d: "mrd_status" in d).any()
    # A column sourced from GDC's own clinical file (not the two flaky
    # third-party ones) must be unaffected.
    assert "vital_status" in ds.samples.columns


def test_target_all_p2_loads_when_only_mrd_fails(tmp_path, monkeypatch):
    """Partial degradation: etp_status survives independently of
    mrd_status failing, since they're two independent additive sources."""
    def _boom_mrd(use_cache=True):
        raise RuntimeError("gdc clinical supplement is down")

    def _ok_etp(use_cache=True):
        return pd.Series({"SAMPLE_A": "ETP", "SAMPLE_B": "non-ETP"})

    _full_mock_load(monkeypatch, tmp_path, etp=_ok_etp, mrd=_boom_mrd)

    ds = gdc_target.load(use_cache=True, with_aux=False)

    group_columns = ds.samples["group_columns"]
    assert not group_columns.apply(lambda d: "mrd_status" in d).any()
    # etp_status should still be present -- its (mocked, but successful)
    # loader ran normally and was not affected by mrd failing.
    assert group_columns.apply(lambda d: "etp_status" in d).any()
    assert group_columns.loc[group_columns.apply(lambda d: d.get("etp_status") == "ETP")].any()
