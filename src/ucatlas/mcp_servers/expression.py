"""expression-mcp: differential expression of ulcerative colitis biopsies (eight cohorts) and the meta-analysis.

Every record carries the evidence_id to cite. Values are those computed by ucatlas.expression (M2)."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ..graph import evidence_id
from ..paths import PROCESSED
from .common import canonical_gene, clean, existing

mcp = FastMCP("expression-mcp")
COHORT_TEXT = {"GSE193677": "MSCCR (Mount Sinai), RNA-seq, discovery", "GSE206285": "UNIFI (ustekinumab trial)",
               "GSE92415": "PURSUIT-SC (golimumab trial)", "GSE16879": "Leuven infliximab cohort",
               "GSE73661": "Leuven infliximab/vedolizumab cohort",
               "GSE109142": "PROTECT, paediatric rectum at diagnosis, RNA-seq",
               "GSE117993": "RISK, paediatric rectum at diagnosis, RNA-seq", "GSE165512": "adult colon, RNA-seq"}


@lru_cache(maxsize=1)
def _de() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "de_results.parquet")


@lru_cache(maxsize=1)
def _meta() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")


def _de_row(r) -> dict:
    return clean({"gene": r.gene, "cohort": r.cohort, "cohort_text": COHORT_TEXT[r.cohort], "contrast": r.contrast,
                  "log2fc": round(r.log2fc, 4), "se": round(r.se, 4), "padj": float(f"{r.padj:.3g}"),
                  "n_disease": int(r.n1), "n_control": int(r.n2),
                  "evidence_id": existing(evidence_id("omics_de", r.cohort, r.contrast, r.gene))})


@mcp.tool()
def get_de_results(gene: str) -> list[dict]:
    """Differential expression of one gene in every cohort and contrast: MSCCR inflamed vs control, non-inflamed vs
    control and inflamed vs non-inflamed in the same patient (paired); UC vs healthy in UNIFI, PURSUIT, GSE16879,
    GSE73661, PROTECT, RISK and GSE165512. log2FC is disease minus control on each cohort's log2 scale (arrays compress fold changes)."""
    g = canonical_gene(gene)
    rows = _de()[_de()["gene"] == g]
    if rows.empty:
        return [{"gene": g, "message": "not measured in the cohorts"}]
    return [_de_row(r) for r in rows.itertuples()]


@mcp.tool()
def meta_analysis(gene: str) -> dict:
    """Random-effects meta-analysis (DerSimonian-Laird) of UC vs control across the eight cohorts: pooled log2FC, BH
    padj, heterogeneity I2, cohorts with a significant change up / down. Cite the returned evidence_id."""
    g = canonical_gene(gene)
    if g not in _meta().index:
        return {"gene": g, "message": "not in the meta-analysis"}
    r = _meta().loc[g]
    return clean({"gene": g, "log2fc": round(r.meta_log2fc, 4), "se": round(r.meta_se, 4),
                  "padj": float(f"{r.meta_padj:.3g}"), "i2": round(r.i2, 3), "cohorts": int(r.k),
                  "cohorts_significant_up": int(r.n_sig_up), "cohorts_significant_down": int(r.n_sig_down),
                  "evidence_id": existing(evidence_id("expression_meta", g))})


@mcp.tool()
def top_genes(direction: str = "up", n: int = 25, min_cohorts: int = 2) -> list[dict]:
    """Strongest meta-analysis changes (direction 'up' or 'down'), significant in the meta-analysis and in at least
    `min_cohorts` cohorts in the same direction."""
    m = _meta().reset_index()
    col = "n_sig_up" if direction == "up" else "n_sig_down"
    m = m[(m["meta_padj"] < 0.05) & (m[col] >= min_cohorts)]
    m = m.nlargest(n, "meta_log2fc") if direction == "up" else m.nsmallest(n, "meta_log2fc")
    return [meta_analysis(g) for g in m["gene"]]


@mcp.tool()
def replication(gene: str) -> dict:
    """Is the change reproduced across the eight cohorts? Per cohort: direction, significance (padj < 0.05) and
    log2FC, with the evidence_id of each cohort's record."""
    g = canonical_gene(gene)
    d = _de()[(_de()["gene"] == g) & _de()["contrast"].isin(["inflamed_vs_control", "uc_vs_control"])]
    rows = [{"cohort": r.cohort, "cohort_text": COHORT_TEXT[r.cohort], "log2fc": round(r.log2fc, 4),
             "padj": float(f"{r.padj:.3g}"), "significant": bool(r.padj < 0.05),
             "evidence_id": existing(evidence_id("omics_de", r.cohort, r.contrast, g))} for r in d.itertuples()]
    return clean({"gene": g, "cohorts": rows, "significant_up": sum(r["significant"] and r["log2fc"] > 0 for r in rows),
                  "significant_down": sum(r["significant"] and r["log2fc"] < 0 for r in rows)})


if __name__ == "__main__":
    mcp.run()
