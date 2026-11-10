"""M4/M5 checks: graph integrity, one test per MCP tool (real data, external calls from the disk cache), stdio smoke test."""
import asyncio
import sys

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from ucatlas.graph import GraphStore
from ucatlas.mcp_servers import expression as ex_mcp
from ucatlas.mcp_servers import external as ext_mcp
from ucatlas.mcp_servers import kg as kg_mcp
from ucatlas.paths import GRAPH_DB, ROOT

pytestmark = pytest.mark.skipif(not GRAPH_DB.exists(), reason="build the graph first")


def test_graph_has_no_claim_without_evidence():
    s = GraphStore().stats()
    assert s["claims_without_resolvable_evidence"] == 0 and s["claims_total"] > 100_000


# ---------------------------------------------------------------- expression-mcp
def test_get_de_results_accepts_legacy_symbol():
    rows = ex_mcp.get_de_results("IL8")  # previous symbol of CXCL8
    assert rows and all(r["gene"] == "CXCL8" for r in rows)
    infl = next(r for r in rows if r["contrast"] == "inflamed_vs_control")
    assert infl["log2fc"] > 1 and infl["evidence_id"].startswith("EV:omics_de:")


def test_meta_analysis():
    r = ex_mcp.meta_analysis("AQP8")
    assert r["log2fc"] < -2 and r["cohorts_significant_down"] >= 4 and r["evidence_id"]


def test_top_genes():
    up = ex_mcp.top_genes("up", 10)
    assert len(up) == 10 and all(r["log2fc"] > 0 for r in up)


def test_replication():
    r = ex_mcp.replication("LCN2")
    assert r["significant_up"] >= 4 and all(c["evidence_id"] for c in r["cohorts"])


# ---------------------------------------------------------------- kg-mcp
def test_neighbors():
    r = kg_mcp.neighbors("TNF", predicate="targets")
    assert r["total"] > 0 and all(n["evidence_ids"] for n in r["neighbors"])


def test_evidence_for_and_get_evidence():
    r = kg_mcp.evidence_for("LCN2", evidence_type="omics_de")
    ev = r["evidence"][0]["evidence_id"]
    assert kg_mcp.check_evidence_ids([ev, "EV:fake:000"]) == {"valid": [ev], "invalid": ["EV:fake:000"]}
    assert kg_mcp.get_evidence([ev])[0]["evidence_id"] == ev


def test_search_literature():
    hits = kg_mcp.search_literature("calprotectin mucosal healing", 5)
    assert hits and hits[0]["evidence_id"] and "pubmed" in hits[0]["url"]


def test_convergence_ranking_tool():
    assert isinstance(kg_mcp.convergence_ranking(5), list)


# ---------------------------------------------------------------- external-mcp (served from the cache)
def test_open_targets_association():
    r = ext_mcp.open_targets_association("IL23R")
    assert r["rank"] <= 50 and r["evidence_id"]


def test_open_targets_target():
    r = ext_mcp.open_targets_target("ITGB7")
    assert any("VEDOLIZUMAB" in d["name"] for d in r["drugs"])


def test_string_and_uniprot():
    assert ext_mcp.string_interactions("ITGB7")
    assert ext_mcp.uniprot_entry("TNF")["accession"] == "P01375"


def test_disease_top_targets():
    top = [r["gene"] for r in ext_mcp.disease_top_targets(25)]
    assert {"TNF", "IL23A", "ITGB7", "JAK1"} <= set(top)


@pytest.mark.parametrize("module,expected", [
    ("ucatlas.mcp_servers.expression", {"get_de_results", "meta_analysis", "top_genes", "replication", "protein_level"}),
    ("ucatlas.mcp_servers.kg", {"neighbors", "evidence_for", "search_literature", "check_evidence_ids", "numeric_annotation", "pubtator_relations"}),
    ("ucatlas.mcp_servers.external", {"open_targets_target", "string_interactions", "uniprot_entry",
                                      "mendelian_randomization", "mr_hits", "mr_gut_hits", "mouse_phenotype"}),
    ("ucatlas.mcp_servers.cell", {"celltype_expression", "composition_effect", "spatial_context", "novel_targets"}),
    ("ucatlas.mcp_servers.pathway", {"pathway_activity", "tf_activity", "communication", "treatment_response"}),
])
def test_stdio_server_lists_tools(module, expected):
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", module], cwd=str(ROOT))
        async with stdio_client(params) as (r, w):
            async with ClientSession(r, w) as s:
                await s.initialize()
                return {t.name for t in (await s.list_tools()).tools}
    assert expected <= asyncio.run(run())
