"""M15: pathway layer — TNFa / NFkB / JAK-STAT active in all five cohorts; the neutrophil OSM -> OSMR circuit to
inflammatory fibroblasts is found; tools return citable records."""
import pandas as pd
import pytest

from ucatlas.mcp_servers import pathway as p
from ucatlas.paths import PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "communication.parquet").exists(), reason="run ucatlas.pathways")


@pytest.mark.parametrize("name", ["TNFa", "NFkB", "JAK-STAT"])
def test_inflammatory_pathways_up_everywhere(name):
    r = p.pathway_activity(name)[0]
    assert r["cohorts_up"] == 5 and r["evidence_id"].startswith("EV:pathway_activity:")


def test_osm_osmr_circuit():
    c = pd.read_parquet(PROCESSED / "communication.parquet")
    x = c[(c.ligand == "OSM") & c.receptor.str.contains("OSMR")]
    assert x.iloc[0]["sender"].startswith("Neutrophil") and "fibroblast" in x.iloc[0]["receiver"].lower()


def test_communication_tool_cites():
    rows = p.communication("CXCL8", 3)
    assert rows and all(r["evidence_id"].startswith("EV:communication:") for r in rows)
