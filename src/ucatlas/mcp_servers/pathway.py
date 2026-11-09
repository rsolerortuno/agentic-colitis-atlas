"""pathway-mcp: pathway and transcription-factor activity in UC mucosa (five cohorts) and ligand-receptor
communication between cell states (single cell with neutrophils). Every record carries the evidence_id to cite."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import PROCESSED
from .common import canonical_gene, clean, store

mcp = FastMCP("pathway-mcp")


def _rec(ev_id: str) -> dict | None:
    r = store().get_evidence([ev_id])
    return clean({**{k: v for k, v in r[0].items() if k not in ("source", "date", "type")}, "evidence_id": ev_id}) if r else None


@mcp.tool()
def pathway_activity(pathway: str | None = None) -> list[dict]:
    """PROGENy pathway activity (Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, PI3K, TGFb, TNFa, Trail,
    VEGF, WNT, p53), inflamed or UC vs control in each of five cohorts (t, BH padj) and the number of cohorts with a
    significant increase or decrease. Without an argument, all 14 pathways."""
    d = pd.read_parquet(PROCESSED / "pathway_activity.parquet")
    names = [pathway] if pathway else sorted(d["source"].unique())
    return [r for n in names if (r := _rec(evidence_id("pathway_activity", n)))]


@lru_cache(maxsize=1)
def _tf() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "tf_activity.parquet").drop_duplicates("source")


@mcp.tool()
def tf_activity(tf: str) -> dict:
    """Transcription-factor activity (CollecTRI regulons, decoupler ulm) of one TF, inflamed or UC vs control in
    five cohorts. Activity is inferred from the TF's target genes, not from its own mRNA."""
    g = canonical_gene(tf)
    return _rec(evidence_id("tf_activity", g)) or {"tf": g, "message": "no CollecTRI regulon for this TF"}


@mcp.tool()
def top_tfs(direction: str = "up", n: int = 20) -> list[dict]:
    """TFs whose inferred activity changes in the most cohorts (then by mean t), up or down."""
    d = _tf()
    col = "cohorts_up" if direction == "up" else "cohorts_down"
    d = d.sort_values([col, "mean_t" if "mean_t" in d else col], ascending=[False, direction != "up"]).head(n)
    return [r for s in d["source"] if (r := _rec(evidence_id("tf_activity", s)))]


@lru_cache(maxsize=1)
def _comm() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "communication.parquet")


@mcp.tool()
def communication(gene: str, n: int = 10) -> list[dict]:
    """Ligand-receptor pairs that involve a gene (as ligand or receptor subunit), in inflamed UC mucosa (Gubatan 2026
    single cell, includes neutrophils): sender and receiver cell states, score (ligand mean CP10K in sender x receptor
    mean CP10K in receiver), and log2 change vs uninflamed mucosa. A score is co-expression, not proof of signalling."""
    g = canonical_gene(gene)
    c = _comm()
    hit = c[c["ligand"].str.split("_").map(lambda x: g in x) | c["receptor"].str.split("_").map(lambda x: g in x)]
    out = []
    for r in hit.head(n).itertuples():
        eid = evidence_id("communication", r.ligand, r.receptor, r.sender, r.receiver)
        rec = _rec(eid)
        if rec:
            out.append(rec)
    return out or [{"gene": g, "message": "no ligand-receptor pair with this gene above the expression thresholds"}]


if __name__ == "__main__":
    mcp.run()
