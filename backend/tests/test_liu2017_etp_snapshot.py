"""
In production Springer Nature answers the server's request for the Liu
2017 supplement with a 3 KB HTML page, not the workbook; validate_excel_bytes
rightly rejects it, and TARGET-ALL-P2 then loaded with no etp_status at all
(0 of 469 samples labelled, vs 190 wherever the download works), so Compare
by etp_status silently dropped TARGET. The loader now falls back to the
committed snapshot of the same labels (METHODS.md 7.10). No network: the
fetch is replaced by one that fails the way production's did.
"""
import pandas as pd

from app.ingest import liu2017_etp


def _html_instead_of_excel(*_args, **_kwargs):
    raise ValueError("response body (3038 bytes) doesn't start with an Excel file signature -- got b'<!DOCTYPE html>\\n'")


def test_falls_back_to_committed_snapshot_when_supplement_is_blocked(tmp_path, monkeypatch):
    monkeypatch.setattr(liu2017_etp, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(liu2017_etp, "fetch_with_retry", _html_instead_of_excel)
    labels = liu2017_etp.load_etp_status(use_cache=False)
    assert len(labels) == 190
    assert labels.value_counts().to_dict() == {"non-ETP": 146, "near-ETP": 25, "ETP": 19}
    assert labels.index.str.startswith("TARGET-10-").all()
    assert labels["TARGET-10-PARNJB"] == "ETP"


def test_snapshot_uses_the_project_label_spelling():
    snap = pd.read_csv(liu2017_etp.SNAPSHOT_PATH)
    assert list(snap.columns) == ["sample_id", "etp_status"]
    assert set(snap["etp_status"]) == {"ETP", "near-ETP", "non-ETP"}
    assert snap["sample_id"].is_unique
