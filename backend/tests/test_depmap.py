import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.ingest.depmap import _release_age_warning


def test_release_age_warning_none_when_current():
    # 26Q1 evaluated as of early 2026 is at most one quarter behind.
    assert _release_age_warning("DepMap 26Q1 Public", today=date(2026, 3, 1)) is None


def test_release_age_warning_none_at_threshold():
    # Exactly 2 quarters behind is still within the loose threshold.
    assert _release_age_warning("DepMap 25Q3 Public", today=date(2026, 3, 15)) is None


def test_release_age_warning_fires_when_stale():
    # 24Q4 evaluated in mid-2026 is far more than 2 quarters behind.
    warning = _release_age_warning("DepMap 24Q4 Public", today=date(2026, 8, 29))
    assert warning is not None
    assert "24Q4" in warning
    assert "quarters behind" in warning


def test_release_age_warning_ignores_unparseable_title():
    assert _release_age_warning("cached", today=date(2026, 8, 29)) is None
