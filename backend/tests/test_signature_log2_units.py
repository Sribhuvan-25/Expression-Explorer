"""
"log2 mean" signature scores applied log2(x+1) to every matrix, so on the
three datasets already stored as log2 (DepMap log2(TPM+1), GDS4299 and
GSE66586 log2 intensity) the score was a mean of log2(log2(x)+1) -- a
GSE66586 sample at ~8 log2 units scored ~3.2 (METHODS.md 2.4).
"""
import numpy as np
import pytest
from fastapi.testclient import TestClient

from app import registry
from app.analysis.signature import is_log2_unit, log2_mean_signature_score
from app.api import main as api_main
from app.models.contract import ExpressionUnit
from tests.test_api import _build_dataset

LOGID = "apitest_log2"


def _log2_dataset():
    ds = _build_dataset()
    ds.source = ds.source.model_copy(update={"dataset_id": LOGID, "expression_unit": ExpressionUnit.LOG2_INTENSITY})
    return ds


@pytest.fixture(scope="module")
def log_client():
    if LOGID not in {d.dataset_id for d in registry.list_descriptors()}:
        registry.register(registry.DatasetDescriptor(dataset_id=LOGID, display_name="log2 test", loader=_log2_dataset))
    return TestClient(api_main.app)


def test_log2_units_are_recognised():
    assert is_log2_unit(ExpressionUnit.LOG2_INTENSITY) and is_log2_unit("log2_tpm")
    assert not is_log2_unit(ExpressionUnit.TPM) and not is_log2_unit("counts")


def test_log2_mean_averages_log2_data_as_is():
    ds = _log2_dataset()
    got = log2_mean_signature_score(ds.matrix, ["MYCN", "DTX1"], already_log2=True)
    np.testing.assert_allclose(got.to_numpy(), ds.matrix.loc[["MYCN", "DTX1"]].mean().to_numpy())


@pytest.mark.parametrize("path", ["signature-score", "rank-signature"])
def test_endpoints_do_not_relog_log2_datasets(log_client, path):
    ds = _log2_dataset()
    expected = ds.matrix.loc[["MYCN", "DTX1"]].mean()
    body = log_client.post(f"/datasets/{LOGID}/{path}", json={"genes": ["MYCN", "DTX1"], "method": "log2_mean"}).json()
    got = body["scores"] if path == "signature-score" else {r["sample_id"]: r["value"] for r in body["rows"]}
    for sid, v in got.items():
        assert v == pytest.approx(expected[sid])


def test_linear_units_are_still_log_transformed(log_client):
    from tests.test_api import DATASET_ID, client  # noqa: F401
    ds = _build_dataset()
    expected = np.log2(ds.matrix.loc[["MYCN", "DTX1"]] + 1).mean()
    if DATASET_ID not in {d.dataset_id for d in registry.list_descriptors()}:
        registry.register(registry.DatasetDescriptor(dataset_id=DATASET_ID, display_name="t", loader=_build_dataset))
    body = log_client.post(f"/datasets/{DATASET_ID}/signature-score", json={"genes": ["MYCN", "DTX1"], "method": "log2_mean"}).json()
    for sid, v in body["scores"].items():
        assert v == pytest.approx(expected[sid])
