"""M11: load the single-cell attribution (Smillie 2019, replication in Garrido-Trigo 2023, neutrophils from Gubatan 2026) into the graph (one sc_summary record per gene, with its replication).

Record fields: the three clusters with the largest share of the gene's transcripts in inflamed, non-inflamed and
healthy mucosa, the dominant lineage in inflamed mucosa with its share, the second dataset's dominant lineage and
whether they agree. Genes: those in the meta-analysis that are measured in Smillie 2019."""
from __future__ import annotations

from datetime import date

import pandas as pd

from .graph import GraphStore, evidence_id
from .paths import PROCESSED

DATASET = "DATASET:SCP259"
DATASET_2 = "DATASET:GARRIDO2023"


def summaries(genes: set[str]) -> dict[str, dict]:
    s = pd.read_parquet(PROCESSED / "sc_celltype_share.parquet")
    s = s[s["gene"].isin(genes)]
    rep_p = PROCESSED / "sc_replication.parquet"
    rep = pd.read_parquet(rep_p).set_index("gene") if rep_p.exists() else pd.DataFrame()
    from .sc_replication import shares_2019
    lin = shares_2019("inflamed")
    lin = lin[lin["gene"].isin(genes)].sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    neu = neutrophil_view(genes)
    out: dict[str, dict] = {}
    for (gene, health), d in s.groupby(["gene", "health"]):
        top = d.nlargest(3, "share")
        out.setdefault(gene, {})[health] = [{"cluster": r.cluster, "compartment": r.compartment,
                                             "share": round(float(r.share), 3),
                                             "pct_expressing": round(float(r.pct_expressing), 1)} for r in top.itertuples()]
    recs = {}
    for gene, h in out.items():
        rec = {"gene": gene, "top_clusters_inflamed": h.get("inflamed", []), "top_clusters_noninflamed": h.get("noninflamed", []),
               "top_clusters_healthy": h.get("healthy", [])}
        if gene in lin.index:
            rec.update(dominant_lineage=lin.loc[gene, "lineage"], dominant_lineage_share=round(float(lin.loc[gene, "share"]), 3))
        if len(rep) and gene in rep.index:
            rec.update(lineage_second_dataset=rep.loc[gene, "lineage_2023"], replicated=bool(rep.loc[gene, "same_lineage"]))
        if gene in neu:
            rec.update(neutrophil_share_inflamed=neu[gene]["neutrophil_share"],
                       top_clusters_inflamed_with_neutrophils=neu[gene]["top"])
        recs[gene] = rec
    return recs


def neutrophil_view(genes: set[str]) -> dict[str, dict]:
    """Gubatan 2026 (SCP3755, captures neutrophils): share of the gene's transcripts in neutrophil states (inflamed
    UC mucosa) and the top three cell states."""
    p = PROCESSED / "sc3_celltype_share.parquet"
    if not p.exists():
        return {}
    s = pd.read_parquet(p)
    s = s[s["gene"].isin(genes) & (s["health"] == "inflamed")]
    out = {}
    for gene, d in s.groupby("gene"):
        top = d.nlargest(3, "share")
        out[gene] = {"neutrophil_share": round(float(d.loc[d["compartment"] == "Neutrophil", "share"].sum()), 3),
                     "top": [{"cell_state": r.cluster, "group": r.compartment, "share": round(float(r.share), 3)}
                             for r in top.itertuples()]}
    return out


def load(store: GraphStore | None = None) -> dict:
    store = store or GraphStore()
    today = date.today().isoformat()
    genes = set(pd.read_parquet(PROCESSED / "meta.parquet")["gene"])
    n = 0
    with store.batch():
        store.add_node(DATASET, "Dataset", "Smillie et al. 2019 colon scRNA-seq (SCP259): UC inflamed, non-inflamed, healthy")
        store.add_node(DATASET_2, "Dataset", "Garrido-Trigo et al. 2023 colon scRNA-seq (CZ CELLxGENE): UC, CD, healthy")
        store.add_node("DATASET:SCP3755", "Dataset", "Gubatan et al. 2026 UC colon scRNA-seq with neutrophils (SCP3755)")
        for gene, rec in summaries(genes).items():
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            ev = store.add_evidence(evidence_id("sc_summary", gene), "sc_summary", DATASET, today, **rec)
            for c in rec["top_clusters_inflamed"][:1]:
                cid = f"CELLTYPE:{c['cluster']}"
                store.add_node(cid, "CellType", c["cluster"], compartment=c["compartment"])
                store.add_claim(node["id"], "expressed_mainly_by", cid, ev)
            if "replicated" in rec:
                store.add_claim(node["id"], "cell_attribution_checked_in", DATASET_2, ev)
            n += 1
    return {"sc_summary_records": n}


if __name__ == "__main__":
    print(load())
