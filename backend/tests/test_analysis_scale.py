"""
Correlation, co-expression and PCA on linear units run on log2(x + 1)
(METHODS.md 7.6a). On linear TPM a few very high values drove Pearson r:
TARGET DTX1 vs NOTCH3 r = 0.113, p = 0.014 against 0.037, p = 0.42 on
log2(TPM + 1).
"""
import numpy as np
import pytest

from app.api import main as api_main
from tests.test_api import DATASET_ID, client  # noqa: F401  (fixture; TPM unit)


def test_pearson_and_points_use_log2_tpm(client):  # noqa: F811
    from scipy import stats

    ds = api_main._get_dataset(DATASET_ID)
    body = client.get(f"/datasets/{DATASET_ID}/correlation", params={"gene_a": "MYCN", "gene_b": "DTX1"}).json()
    a, b = np.log2(ds.matrix.loc["MYCN"] + 1), np.log2(ds.matrix.loc["DTX1"] + 1)
    assert body["coefficient"] == pytest.approx(stats.pearsonr(a, b)[0])
    assert body["value_unit"] == "log2_tpm"
    assert sorted(p["x"] for p in body["points"]) == pytest.approx(sorted(a))


def test_rank_coefficients_do_not_change(client):  # noqa: F811
    from scipy import stats

    ds = api_main._get_dataset(DATASET_ID)
    body = client.get(
        f"/datasets/{DATASET_ID}/correlation", params={"gene_a": "MYCN", "gene_b": "DTX1", "method": "spearman"}
    ).json()
    assert body["coefficient"] == pytest.approx(stats.spearmanr(ds.matrix.loc["MYCN"], ds.matrix.loc["DTX1"])[0])


def test_pca_runs_on_log2_tpm(client):  # noqa: F811
    from app.analysis.dimensionality import pca

    ds = api_main._get_dataset(DATASET_ID)
    genes = ["MYCN", "DTX1", "NOTCH1"]
    body = client.post(f"/datasets/{DATASET_ID}/pca", json={"genes": genes}).json()
    expected = pca(np.log2(ds.matrix.loc[genes] + 1), genes, ds.samples)
    assert body["variance_explained"] == pytest.approx(expected["variance_explained"])
    assert body["value_unit"] == "log2_tpm"
