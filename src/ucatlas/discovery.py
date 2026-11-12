"""M13: screen for new candidates (deterministic, no LLM), the same idea as the kidney atlas adapted to UC.

Pool: genes up in UC mucosa (meta-analysis BH < 0.05, log2FC >= 1, significant up in >= 2 cohorts).
Filters, each recorded per gene:
- tissue change: the composition-adjusted change survives ('intrinsic' or 'partly_composition', M12);
- cell source: dominant lineage in inflamed mucosa (Smillie 2019) and neutrophil share (Gubatan 2026);
- novelty: <= 15 Europe PMC records mentioning the gene with ulcerative colitis, and no clinical-stage drug
  (Open Targets, any indication);
- tractability: Open Targets small-molecule or antibody tractability;
- safety: highest RNA level in vital organs (heart muscle, brain regions, liver, lung, kidney; HPA consensus)
  <= 50 nTPM.
Status: 'first-in-class opportunity' when all pass; otherwise the first failing filter names the status.
Neutrophil-made genes (added 2026-11-12): when neutrophils make more than half of a gene's transcripts in inflamed
mucosa (Gubatan 2026), the tissue-change filter is not applied. The composition adjustment removes the neutrophil
influx these genes are about; they go on through the other filters as neutrophil-directed candidates (blocking them
acts on neutrophil arrival or activity, like CXCR2), and the validation lists measure their change per neutrophil.
The neutrophil_directed column is kept out of the graph record (the agents cite these records).
Output: data/processed/discovery.parquet.
"""
from __future__ import annotations

import io
import zipfile

import pandas as pd

from .external import DISEASE_TERMS, europepmc_count, ot_target
from .paths import PROCESSED, RAW

VITAL = ("heart muscle", "cerebral cortex", "cerebellum", "liver", "lung", "kidney")
MAX_PAPERS, MAX_VITAL_NTPM = 15, 50
NEUTROPHIL_MADE = 0.5   # as in validation_shortlist


def vital_expression() -> pd.Series:
    with zipfile.ZipFile(RAW / "extra" / "hpa_rna_tissue_consensus.tsv.zip") as z:
        t = pd.read_csv(io.BytesIO(z.read(z.namelist()[0])), sep="\t")
    t = t[t["Tissue"].isin(VITAL)]
    return t.groupby("Gene name")["nTPM"].max()


def pool() -> pd.DataFrame:
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    return m[(m["meta_padj"] < 0.05) & (m["meta_log2fc"] >= 1) & (m["n_sig_up"] >= 2)].copy()


def run() -> pd.DataFrame:
    from .graph import GraphStore, evidence_id
    p = pool()
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    vital = vital_expression()
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    ens = dict(zip(hg["symbol"], hg["ensembl_gene_id"]))
    store = GraphStore()
    rows = []
    for r in p.itertuples():
        g = r.gene
        sc = store.get_evidence([evidence_id("sc_summary", g)])
        sc = sc[0] if sc else {}
        t = ot_target(ens[g]) if g in ens else None
        papers = europepmc_count(f'"{g}" AND {DISEASE_TERMS}')
        rows.append({"gene": g, "meta_log2fc": r.meta_log2fc, "cohorts_up": int(r.n_sig_up),
                     "composition": comp["verdict"].get(g), "retained": comp["retained"].get(g),
                     "lineage": sc.get("dominant_lineage"), "neutrophil_share": sc.get("neutrophil_share_inflamed"),
                     "uc_papers": papers,
                     "clinical_drug": bool(t and t["has_approved_or_clinical_drug"]),
                     "tractable": sorted(set(t["tractable_modalities"]) & {"SM", "AB"}) if t else [],
                     "vital_max_ntpm": float(vital.get(g)) if g in vital.index else None})
    d = pd.DataFrame(rows)

    d["neutrophil_directed"] = d["neutrophil_share"].fillna(0) > NEUTROPHIL_MADE

    def status(x):
        if x["composition"] not in ("intrinsic", "partly_composition") and not x["neutrophil_directed"]:
            return "cell influx (composition)"
        if x["clinical_drug"]:
            return "known druggable target"
        if x["uc_papers"] > MAX_PAPERS:
            return "well studied in UC"
        if not x["tractable"]:
            return "novel, not tractable"
        if x["vital_max_ntpm"] is None or x["vital_max_ntpm"] > MAX_VITAL_NTPM:
            return "novel, vital-organ expression"
        return "first-in-class opportunity"

    d["status"] = d.apply(status, axis=1)
    d.to_parquet(PROCESSED / "discovery.parquet", index=False)
    return d


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    n = 0
    with store.batch():
        store.add_node("ANALYSIS:discovery", "Dataset", "Discovery screen for new UC candidates (M13)")
        for r in d.to_dict("records"):
            node = store.resolve(r["gene"])
            if not node or node["type"] != "Gene":
                continue
            rec = {k: (None if isinstance(v, float) and pd.isna(v) else v) for k, v in r.items() if k != "neutrophil_directed"}
            rec["tractable"] = list(rec["tractable"])
            ev = store.add_evidence(evidence_id("discovery_profile", r["gene"]), "discovery_profile", "ANALYSIS:discovery",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "screened_in", "ANALYSIS:discovery", ev)
            n += 1
    return n


if __name__ == "__main__":
    d = run()
    print(d["status"].value_counts().to_dict(), "records:", load())
    print(d[d["status"] == "first-in-class opportunity"].sort_values("meta_log2fc", ascending=False)
          [["gene", "meta_log2fc", "lineage", "uc_papers", "tractable", "vital_max_ntpm", "retained"]].head(25).to_string())
