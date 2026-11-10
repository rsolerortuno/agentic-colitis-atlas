"""kg-mcp: query the knowledge graph (GraphStore) and the literature chunk index."""
from __future__ import annotations

import json
import os

from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import RUNS
from ..search import search_chunks
from .common import canonical_gene, clean, store

mcp = FastMCP("kg-mcp")


def _resolve(entity: str) -> dict | None:
    s = store()
    return s.resolve(entity) or s.resolve(canonical_gene(entity))


@mcp.tool()
def neighbors(entity: str, predicate: str | None = None, limit: int = 50) -> dict:
    """Neighbors of a graph entity (gene symbol, node id such as 'MONDO:0005005', 'PMID:123', drug name).
    Each neighbor lists the predicate, direction and the evidence_ids supporting the link.
    Optional predicate filter, e.g. 'targets', 'interacts_with', 'in_pathway', 'upregulated_in'."""
    n = _resolve(entity)
    if not n:
        return {"error": f"entity '{entity}' not found"}
    nb = store().neighbors(n["id"], limit=500)
    if predicate:
        nb = [x for x in nb if x["predicate"] == predicate]
    for x in nb:
        x["evidence_ids"] = x["evidence_ids"][:5]
    return clean({"entity": {k: n[k] for k in ("id", "type", "name")}, "neighbors": nb[:limit],
                  "total": len(nb)})


@mcp.tool()
def evidence_for(entity: str, evidence_type: str | None = None, limit: int = 50) -> dict:
    """All evidence records attached to an entity's relations (omics statistics, literature sentences with
    PMID, Open Targets scores, drugs, STRING, UniProt). Optional evidence_type filter, e.g.
    'literature_relation', 'omics_de', 'open_targets_drug', 'open_targets_association'."""
    n = _resolve(entity)
    if not n:
        return {"error": f"entity '{entity}' not found"}
    ev = store().evidence_for(n["id"], limit=2000)
    ev = [e for e in ev if e.get("type") not in ("provenance", "pubtator_annotation")]
    if evidence_type:
        ev = [e for e in ev if e.get("type") == evidence_type]
    seen, out = set(), []
    for e in ev:
        if e["evidence_id"] not in seen:
            seen.add(e["evidence_id"])
            out.append(e)
    return clean({"entity": {k: n[k] for k in ("id", "type", "name")}, "evidence": out[:limit], "total": len(out)})


@mcp.tool()
def search_literature(query: str, k: int = 10) -> list[dict]:
    """BM25 search over sentences of the ingested ulcerative colitis abstracts. Returns sentence text, PMID and the
    evidence_id of the sentence chunk (cite it, or cite literature_relation evidence from evidence_for)."""
    hits = search_chunks(query, k=min(k, 30))
    for h in hits:
        ev = evidence_id("chunk", h["chunk_id"])
        h["evidence_id"] = ev if store().check_evidence_ids([ev])["valid"] else None
        h["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{h['pmid']}/"
    return hits


@mcp.tool()
def check_evidence_ids(ids: list[str]) -> dict:
    """Check which evidence ids exist in the graph. Use before returning findings."""
    return store().check_evidence_ids(ids)


@mcp.tool()
def get_evidence(ids: list[str]) -> list[dict]:
    """Full records for a list of evidence ids (max 50)."""
    return clean(store().get_evidence(ids[:50]))


@mcp.tool()
def numeric_annotation(evidence_id: str, field: list, display: str, unit: str = "native") -> dict:
    """Build the exact annotation a number needs to pass verification (docs/numeric-claims.md). field = path of keys
    into the cited record, e.g. ["log2fc"] or ["arms", "sunitinib", "log_hr_per_sd"]; display = the number as you will
    print it, e.g. "2.27" (sign and decimals must match the record); unit = "native", "percent" (registered fractions)
    or "native_percent" (fields already 0-100). Returns `annotation` (paste it once, unchanged, into the statement where
    the number goes) and `numeric_claim` (put it in the finding's "numeric_claims" list), or `error` explaining why the
    number does not match the record. Use it for every number you write."""
    from ..semantic import numeric_claim_text
    recs = store().get_evidence([evidence_id])
    if not recs:
        return {"error": f"evidence id {evidence_id} does not exist"}
    claim = {"evidence_id": evidence_id, "field": field, "display": display, "unit": unit}
    try:
        return {"annotation": numeric_claim_text(claim, recs[0]), "numeric_claim": claim}
    except (ValueError, TypeError) as exc:
        return {"error": str(exc), "record": clean(recs[0])}


def _latest_convergence() -> dict | None:
    from ..runs import run_dir
    if os.environ.get("ATLAS_RUN_DIR") and (run_dir() / "convergence.json").exists():
        return json.loads((run_dir() / "convergence.json").read_text())
    from ..runs import current_run
    p = current_run() / "convergence.json"
    return json.loads(p.read_text()) if p.exists() else None


@mcp.tool()
def convergence_ranking(n: int = 25, profile: str = "convergence") -> list[dict]:
    """Deterministic ranking from the latest scripts/integrate.py run.

    profile: 'convergence' (all evidence streams: expression meta-analysis, replication across cohorts, inflamed vs
    non-inflamed mucosa, genetics, literature, tractability). 'target' (expression up, composition-retained tissue change, genetics, predicted tractability) and 'biomarker'
    (any change, replication, disease activity, inflamed vs non-inflamed) leave out literature and drug history. Each row carries its components, the meta-analysis and paired values, Open Targets scores and the
    evidence ids per stream. When present, `rank_interval_90` gives the 5th-95th percentile rank and probability of
    staying in the top 10/25/100 when the weights are randomly perturbed (robustness)."""
    conv = _latest_convergence()
    if not conv:
        return [{"message": "no convergence run yet; the orchestrator runs scripts/integrate.py"}]
    if profile == "convergence":
        return conv["ranking"][:n]
    rows = (conv.get("profiles") or {}).get(profile)
    if rows is None:
        return [{"message": f"unknown profile '{profile}'; use convergence, target or biomarker"}]
    return rows[:n]


@mcp.tool()
def gene_scores(gene: str) -> dict:
    """Where one gene sits in the convergence ranking and in the target / biomarker profiles (rank, score,
    components, robustness interval). Returns only what the latest run computed; absent = outside the stored top."""
    conv = _latest_convergence()
    if not conv:
        return {"message": "no convergence run yet"}
    g = canonical_gene(gene)
    out = {"gene": g}
    for name, rows in [("convergence", conv["ranking"])] + list((conv.get("profiles") or {}).items()):
        hit = next((r for r in rows if r["gene"] == g), None)
        out[name] = hit if hit else "not in the stored top of this ranking"
    return clean(out)


@mcp.tool()
def pubtator_relations(gene: str) -> dict:
    """PubTator3's machine-extracted relations for a gene in the local literature corpus (585 abstracts): number of
    relations and publications, relation types to ulcerative colitis / IBD / colitis with their PMIDs, and the most
    frequent partners. Machine extraction, not verbatim-checked: cite it as literature context, and prefer the
    verbatim-checked relations (evidence_for with evidence_type literature_relation) for specific claims."""
    g = canonical_gene(gene)
    eid = evidence_id("pubtator_relations", g)
    r = store().get_evidence([eid])
    if not r:
        return {"gene": g, "message": "no PubTator3 relation for this gene in the local corpus"}
    return clean({**{k: v for k, v in r[0].items() if k not in ("source", "date", "type")}, "evidence_id": eid})


if __name__ == "__main__":
    mcp.run()
