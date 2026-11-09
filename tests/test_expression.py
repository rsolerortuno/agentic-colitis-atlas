"""M2: positive controls change in the expected direction in at least two cohorts; statistics match scipy."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from ucatlas import expression
from ucatlas.paths import PROCESSED

UP = ["S100A8", "LCN2", "REG4", "DUOX2", "CXCL8", "MMP3", "IL1B", "OSM"]
DOWN = ["AQP8", "SLC26A2", "HMGCS2", "ABCG2"]
have = (PROCESSED / "de_results.parquet").exists()


@pytest.mark.skipif(not have, reason="run ucatlas.expression first")
@pytest.mark.parametrize("gene,sign", [(g, 1) for g in UP] + [(g, -1) for g in DOWN])
def test_positive_controls(gene, sign):
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    r = de[(de.gene == gene) & de.contrast.isin(["inflamed_vs_control", "uc_vs_control"]) & (de.padj < 0.05)]
    assert (np.sign(r.log2fc) == sign).sum() >= 2


def test_welch_matches_scipy():
    rng = np.random.default_rng(0)
    x = pd.DataFrame(rng.normal(size=(5, 12)), columns=[f"s{i}" for i in range(12)])
    r = expression.welch(x, [f"s{i}" for i in range(6)], [f"s{i}" for i in range(6, 12)])
    p = stats.ttest_ind(x.iloc[:, :6], x.iloc[:, 6:], axis=1, equal_var=False).pvalue
    assert np.allclose(r["p"], p)


def test_random_effects_equals_fixed_when_homogeneous():
    est, se, p, i2, k = expression.random_effects(np.array([1.0, 1.0, 1.0]), np.array([0.1, 0.2, 0.3]))
    assert est == pytest.approx(1.0) and i2 == 0 and k == 3
