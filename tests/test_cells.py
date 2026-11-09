"""M11/M12: cell attribution of control genes, replication, neutrophils and the composition adjustment."""
import pandas as pd
import pytest

from ucatlas.mcp_servers import cell as cell_mcp
from ucatlas.paths import GRAPH_DB, PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "sc3_celltype_share.parquet").exists() or not GRAPH_DB.exists(),
                                reason="single-cell layer not built")


@pytest.mark.parametrize("gene,lineage", [("LCN2", "epithelial"), ("DUOX2", "epithelial"), ("IL11", "stromal"),
                                          ("CHI3L1", "stromal"), ("CXCL8", "myeloid")])
def test_attribution_and_replication(gene, lineage):
    r = cell_mcp.celltype_expression(gene)
    assert r["dominant_lineage"] == lineage and r["replicated"] and r["evidence_id"]


def test_neutrophil_genes():
    rows = {r["gene"]: r for r in cell_mcp.neutrophil_genes(["CXCR2", "S100A8", "AQP8"])}
    assert rows["CXCR2"]["neutrophil_share"] > 0.8 and rows["S100A8"]["neutrophil_share"] > 0.5
    assert rows["AQP8"]["neutrophil_share"] < 0.1


def test_composition_separates_influx_from_tissue_change():
    assert cell_mcp.composition_effect("CXCL8")["verdict"] in ("composition", "partly_composition")
    assert cell_mcp.composition_effect("AQP8")["verdict"] == "intrinsic"


def test_replication_share():
    r = pd.read_parquet(PROCESSED / "sc_replication.parquet")
    assert r["same_lineage"].mean() > 0.75
