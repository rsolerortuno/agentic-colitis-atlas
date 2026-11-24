"""Citable ranks: where a gene stands among the others, so comparisons can be cited (added 2026-11-23).

A reading check of run 2026-11-14 (claim_check.py) found about two dozen comparisons across genes ("the largest
pooled decrease", "the only intrinsic epithelial hit", "the leading plasma MR hit", "the highest vital-organ RNA of the
shortlist") that cited only the gene's own record, which cannot show a rank. These records give the rank within a
named set, with the set size and the leading genes, for:

- meta_up / meta_down: pooled UC vs control log2FC among genes with meta-analysis BH < 0.05 (meta.parquet), increases
  and decreases separately; records for the 150 strongest each way.
- first_in_class: the first-in-class candidates of the discovery screen (discovery.parquet), ordered by pooled log2FC
  as in the report; per gene its rank by log2FC, vital-organ RNA (highest first), retained share after composition
  adjustment (highest first, overall and within its lineage), and cohorts up (fewest first), and the genes that share
  its lineage and composition verdict.
- mr_plasma / mr_gut: rank by p among the significant (q < 0.05) Mendelian randomisation genes of the layer, and among
  those that colocalise (gut: colon tissue and gut cell types together, best p per gene).
- ot_association: Open Targets overall association score with UC among the genes the graph holds.

Ties share the better rank. Evidence type rank_context, one record per gene and set, linked gene -> ranked_in ->
ANALYSIS:rank_context.
"""
from __future__ import annotations

import hashlib
import json

import pandas as pd

from .paths import PROCESSED

SOURCE = "ANALYSIS:rank_context"
TOP_META = 150
LEADERS = 5


def _rank(s: pd.Series, ascending: bool) -> pd.Series:
    return s.rank(ascending=ascending, method="min").astype(int)


def _leaders(d: pd.DataFrame, col: str, ascending: bool, n: int = LEADERS) -> list[dict]:
    t = d.sort_values(col, ascending=ascending).head(n)
    return [{"gene": g, col: round(float(v), 4)} for g, v in zip(t["gene"], t[col])]


def records() -> list[dict]:
    out = []
    # ---------------------------------------------------------------- meta-analysis
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[m["meta_padj"] < 0.05]
    for name, part, asc in (("meta_up", m[m.meta_log2fc > 0], False), ("meta_down", m[m.meta_log2fc < 0], True)):
        part = part.assign(rank=_rank(part["meta_log2fc"], asc)).sort_values("rank")
        lead = _leaders(part, "meta_log2fc", asc)
        for r in part.head(TOP_META).itertuples():
            out.append({"gene": r.gene, "set": name,
                        "set_description": f"genes {'up' if name == 'meta_up' else 'down'} in UC vs control in the "
                                           "eight-cohort meta-analysis (BH < 0.05), by pooled log2FC",
                        "rank": int(r.rank), "of": len(part), "meta_log2fc": round(float(r.meta_log2fc), 4),
                        "leaders": lead})
    # ---------------------------------------------------------------- discovery shortlist
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    f = d[d.status == "first-in-class opportunity"].copy()
    f["r_lfc"] = _rank(f["meta_log2fc"], False)
    f["r_vital"] = _rank(f["vital_max_ntpm"].fillna(-1), False)
    f["r_retained"] = _rank(f["retained"].fillna(-1), False)
    f["r_cohorts_fewest"] = _rank(f["cohorts_up"], True)
    f["r_retained_in_lineage"] = f.groupby("lineage")["retained"].rank(ascending=False, method="min")
    for r in f.itertuples():
        same = f[(f.lineage == r.lineage) & (f.composition == r.composition)].sort_values("r_lfc")
        out.append({"gene": r.gene, "set": "first_in_class",
                    "set_description": "first-in-class candidates of the discovery screen, ordered by pooled log2FC",
                    "of": len(f), "rank_by_log2fc": int(r.r_lfc),
                    "rank_by_vital_organ_rna_highest_first": int(r.r_vital),
                    "rank_by_retained_share_highest_first": int(r.r_retained),
                    "rank_by_cohorts_up_fewest_first": int(r.r_cohorts_fewest),
                    "lineage": r.lineage, "composition": r.composition,
                    "rank_by_retained_share_within_lineage": None if pd.isna(r.r_retained_in_lineage)
                    else int(r.r_retained_in_lineage),
                    "candidates_in_lineage": int((f.lineage == r.lineage).sum()),
                    "same_lineage_and_composition": [{"gene": g, "rank_by_log2fc": int(k)}
                                                     for g, k in zip(same["gene"], same["r_lfc"])]})
    # ---------------------------------------------------------------- Mendelian randomisation
    for name, fn in (("mr_plasma", "mr_pqtl.parquet"), ("mr_gut", "mr_colon.parquet")):
        t = pd.read_parquet(PROCESSED / fn)
        t = t[t["q"] < 0.05].sort_values("p").drop_duplicates("gene")
        col = t["colocalises"].fillna(False).astype(bool) if "colocalises" in t else None
        if col is None:   # gut layers: colocalisation stored in the graph record, read from mr_colon's own column
            col = t.get("coloc_h4", pd.Series(0, index=t.index)).fillna(0).ge(0.8) | \
                t.get("coloc_clpp", pd.Series(0, index=t.index)).fillna(0).ge(0.01)
        t = t.assign(rank=_rank(t["p"], True), coloc=col.values)
        tc = t[t.coloc].assign(rank_c=lambda x: _rank(x["p"], True))
        lead = [{"gene": g, "p": float(p)} for g, p in zip(t["gene"].head(LEADERS), t["p"].head(LEADERS))]
        for r in t.itertuples():
            out.append({"gene": r.gene, "set": name,
                        "set_description": f"genes significant (q < 0.05) in {'plasma pQTL' if name == 'mr_plasma' else 'gut eQTL (colon tissue and gut cell types)'} "
                                           "Mendelian randomisation for UC, by p (best per gene)",
                        "rank": int(r.rank), "of": len(t), "p": float(r.p), "colocalises": bool(r.coloc),
                        "rank_among_colocalising": int(tc.loc[r.Index, "rank_c"]) if r.Index in tc.index else None,
                        "colocalising": len(tc), "leaders": lead})
    # ---------------------------------------------------------------- Open Targets association
    from .graph import GraphStore
    rows = [json.loads(p) for (p,) in GraphStore().conn.execute(
        "select props from evidence where type = 'open_targets_association'")]
    o = pd.DataFrame(rows).dropna(subset=["gene", "score"]).sort_values("score", ascending=False).drop_duplicates("gene")
    o["rank"] = _rank(o["score"], False)
    lead = _leaders(o, "score", False)
    for r in o.itertuples():
        out.append({"gene": r.gene, "set": "ot_association",
                    "set_description": "Open Targets overall association score with UC, among the genes in the graph",
                    "rank": int(r.rank), "of": len(o), "score": float(r.score), "leaders": lead})
    return out


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    n = 0
    with store.batch():
        store.add_node(SOURCE, "Dataset", "Ranks of genes within named sets (citable comparisons)")
        for rec in records():
            node = store.resolve(rec["gene"])
            if not node or node["type"] != "Gene":
                continue
            digest = hashlib.sha256(json.dumps(rec, sort_keys=True, default=str).encode()).hexdigest()
            ev = store.add_evidence(evidence_id("rank_context", rec["set"], rec["gene"]), "rank_context", SOURCE,
                                    date.today().isoformat(), **rec, result_sha256=digest)
            store.add_claim(node["id"], "ranked_in", SOURCE, ev)
            n += 1
    return n


if __name__ == "__main__":
    recs = records()
    print(pd.Series([r["set"] for r in recs]).value_counts().to_dict(), "loaded:", load())
