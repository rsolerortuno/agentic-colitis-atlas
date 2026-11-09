"""cell-mcp: which cell makes each gene in ulcerative colitis mucosa (three single-cell datasets) and, from M12, how
much of a bulk change is explained by cell composition. Every record carries the evidence_id to cite."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import PROCESSED
from .common import canonical_gene, clean, store

mcp = FastMCP("cell-mcp")


def _record(ev_type: str, gene: str) -> dict | None:
    eid = evidence_id(ev_type, gene)
    rec = store().get_evidence([eid])
    if not rec:
        return None
    out = {k: v for k, v in rec[0].items() if k not in ("source", "date", "type")}
    return clean({**out, "evidence_id": eid})


@mcp.tool()
def celltype_expression(gene: str) -> dict:
    """Which cell types make a gene's transcripts in UC colon. Smillie 2019 (SCP259; 51 clusters, inflamed,
    non-inflamed and healthy mucosa; no neutrophils captured): top clusters by share of the gene's transcripts and
    % of cells expressing; dominant lineage in inflamed mucosa and whether a second dataset (Garrido-Trigo 2023)
    agrees. Gubatan 2026 (SCP3755; inflamed and uninflamed UC, includes five neutrophil states): share of transcripts
    made by neutrophils and the top cell states. Shares compare cell types, not tissue abundance (dissociation
    protocols enrich some fractions)."""
    g = canonical_gene(gene)
    return _record("sc_summary", g) or {"gene": g, "message": "not measured in the single-cell datasets"}


@lru_cache(maxsize=1)
def _neu() -> pd.DataFrame:
    s = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    s = s[(s["health"] == "inflamed")]
    return s.groupby("gene").apply(lambda d: d.loc[d["compartment"] == "Neutrophil", "share"].sum()).rename("neutrophil_share")


@mcp.tool()
def neutrophil_genes(genes: list[str]) -> list[dict]:
    """For a list of genes: the share of each gene's transcripts made by neutrophils in inflamed UC mucosa
    (Gubatan 2026). A high share means a bulk increase may simply reflect neutrophil influx. Cite the gene's
    sc_summary evidence_id."""
    n = _neu()
    out = []
    for gene in genes:
        g = canonical_gene(gene)
        out.append({"gene": g, "neutrophil_share": round(float(n[g]), 3) if g in n.index else None,
                    "evidence_id": evidence_id("sc_summary", g) if _record("sc_summary", g) else None})
    return out


@mcp.tool()
def composition_effect(gene: str) -> dict:
    """M12: does the bulk change survive adjustment for cell composition (neutrophil, myeloid, plasma, T cell,
    epithelial, stromal proportions estimated from the biopsies)? Returns the unadjusted and adjusted log2FC of
    inflamed UC vs control in the discovery cohort, and the verdict ('intrinsic' when the adjusted change stays
    significant in the same direction, 'composition' when it does not)."""
    g = canonical_gene(gene)
    return _record("composition_effect", g) or {"gene": g, "message": "composition effect not computed for this gene"}


@mcp.tool()
def spatial_context(gene: str) -> dict:
    """In-situ view from Xenium 5K sections of UC and healthy rectum (GSE312420): share of the gene's transcripts per
    cell lineage by state, and lineages enriched among the 10 nearest neighbours of expressing cells (permutation
    p < 0.05, enrichment > 1.2). Sections with few typed cells are left out and named. Panel sensitivity is low:
    absence of enrichment is not evidence of absence."""
    g = canonical_gene(gene)
    return _record("spatial_summary", g) or {"gene": g, "message": "not on the Xenium 5K panel or too few counts"}


@mcp.tool()
def discovery_profile(gene: str) -> dict:
    """M13 discovery screen record of a gene up in UC: composition verdict and retained share, dominant lineage,
    neutrophil share, UC papers in Europe PMC, clinical drug, tractability, highest vital-organ RNA (HPA) and status
    ('first-in-class opportunity' when it passes every filter)."""
    g = canonical_gene(gene)
    return _record("discovery_profile", g) or {"gene": g, "message": "not in the discovery pool (not up in UC)"}


@mcp.tool()
def novel_targets(n: int = 20) -> list[dict]:
    """The first-in-class opportunities of the discovery screen, strongest meta-analysis increase first."""
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    d = d[d["status"] == "first-in-class opportunity"].sort_values("meta_log2fc", ascending=False).head(n)
    return [r for g in d["gene"] if (r := _record("discovery_profile", g))]


if __name__ == "__main__":
    mcp.run()
