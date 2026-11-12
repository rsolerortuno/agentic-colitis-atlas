"""M13: the discovery screen never calls a known UC drug target or a cell-influx gene 'first-in-class'."""
import pandas as pd
import pytest

from ucatlas.paths import PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "discovery.parquet").exists(), reason="run ucatlas.discovery")


def test_first_in_class_rules():
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    f = d[d.status == "first-in-class opportunity"]
    assert len(f) > 0
    assert not f["clinical_drug"].any() and (f["uc_papers"] <= 15).all() and (f["vital_max_ntpm"] <= 50).all()
    # tissue change, or neutrophil-directed (over half of the RNA in neutrophils: the composition filter is skipped)
    assert (f["composition"].isin(["intrinsic", "partly_composition"]) | (f["neutrophil_share"] > 0.5)).all()
    known = {"TNF", "IL23A", "IL12B", "ITGB7", "ITGA4", "JAK1", "S1PR1", "TNFSF15", "CXCR2", "CXCL8"}
    assert not (set(f["gene"]) & known)
