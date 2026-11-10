"""Protein-level replication in UC colon mucosa (added 2026-11-10, third round; kidney atlas: CPTAC confirmatory cohort).

Source: Schniers et al. 2019, Clin Proteomics 16:4 (additional file 2): label-free bottom-up proteomics (MaxQuant,
Perseus) of colon mucosa biopsies, inflamed mucosa of untreated UC patients at diagnosis against healthy controls;
8,562 proteins identified, 6,818 quantified in > 70 % of samples. The authors' table gives, per protein group, the
UC/healthy ratio, -log10 p of a two-sample test on z-scored intensities, and their significance call (S0 = 2,
FDR 0.01, proteins quantified in > 70 % of samples).

Here: the first gene name of each protein group is mapped to an approved HGNC symbol; log2 ratio and p are taken as
published; Benjamini-Hochberg over the proteins quantified in > 70 % of samples. Replication of the RNA change: a gene
replicates at protein level when its meta-analysis change (padj < 0.05, |log2FC| >= 1) and the protein change
(padj < 0.05) have the same sign. The study's patients are not those of the RNA cohorts, so this is an independent
check, with a different molecule and method.
Outputs: data/processed/proteomics.parquet, proteomics_summary.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED, RAW

SRC = RAW / "extra" / "proteomics_schniers2019.xlsx"
STUDY = "Schniers et al. 2019, Clin Proteomics 16:4"


def run() -> dict:
    from .expression import resolve_symbol
    d = pd.read_excel(SRC, "All proteins", header=7)
    cols = {c.split(" (")[0].split("  ")[0].strip(): c for c in d.columns}
    g = d[cols["Gene names"]].astype(str).str.split(";").str[0].str.strip()
    out = pd.DataFrame({"gene": g.map(lambda s: resolve_symbol(s) if s and s != "nan" else None),
                        "gene_name_published": g, "ratio_uc_h": pd.to_numeric(d[cols["Ratio UC/H"]], errors="coerce"),
                        "mlog10p": pd.to_numeric(d[[c for c in d.columns if c.startswith("-Log p-value")][0]], errors="coerce"),
                        "n_samples": pd.to_numeric(d[[c for c in d.columns if c.startswith("Number of samples")][0]], errors="coerce"),
                        "quantified_70": d[[c for c in d.columns if c.startswith("Quantified in")][0]].astype(str).str.startswith("+"),
                        "significant_authors": d[[c for c in d.columns if c.startswith("Significant")][0]].astype(str).str.startswith("+"),
                        "protein_ids": d[cols["Protein IDs"]].astype(str)})
    out = out.dropna(subset=["gene", "ratio_uc_h"])
    out = out[out.ratio_uc_h > 0]
    out["log2_ratio"] = np.log2(out.ratio_uc_h)
    out["p"] = 10 ** (-out.mlog10p)
    # one row per gene: the protein group with most samples, then the strongest
    out = out.sort_values(["gene", "n_samples", "mlog10p"], ascending=[True, False, False]).drop_duplicates("gene")
    q = out[out.quantified_70 & out.p.notna()].copy()
    q["padj"] = multipletests(q.p, method="fdr_bh")[1]
    out = out.merge(q[["gene", "padj"]], on="gene", how="left")
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    out["meta_log2fc"] = out.gene.map(meta.meta_log2fc)
    out["meta_padj"] = out.gene.map(meta.meta_padj)
    rna_sig = (out.meta_padj < 0.05) & (out.meta_log2fc.abs() >= 1)
    prot_sig = out.padj < 0.05
    out["replicates"] = np.where(rna_sig & out.padj.notna(),
                                 np.where(prot_sig & (np.sign(out.log2_ratio) == np.sign(out.meta_log2fc)), "yes", "no"), "")
    out.to_parquet(PROCESSED / "proteomics.parquet", index=False)
    both = out.dropna(subset=["meta_log2fc", "padj"])
    tested = both[(both.meta_padj < 0.05) & (both.meta_log2fc.abs() >= 1)]
    rep = (tested.padj < 0.05) & (np.sign(tested.log2_ratio) == np.sign(tested.meta_log2fc))
    large = tested[tested.meta_log2fc.abs() >= 2]
    rep_large = (large.padj < 0.05) & (np.sign(large.log2_ratio) == np.sign(large.meta_log2fc))
    s = {"study": STUDY, "proteins_mapped": int(len(out)), "quantified_70": int(out.quantified_70.sum()),
         "significant_padj05": int((out.padj < 0.05).sum()), "significant_authors": int(out.significant_authors.sum()),
         "rna_changes_tested": int(len(tested)), "replicated": int(rep.sum()),
         "replication_rate": round(float(rep.mean()), 3) if len(tested) else None,
         "replication_rate_log2fc_ge2": round(float(rep_large.mean()), 3) if len(large) else None,
         "spearman_log2fc": round(float(stats.spearmanr(both.meta_log2fc, both.log2_ratio).statistic), 3),
         "genes_compared": int(len(both))}
    (PROCESSED / "proteomics_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "proteomics.parquet")
    d = d[d.padj.notna()]
    n = 0
    with store.batch():
        store.add_node("DATASET:PROTEOMICS_SCHNIERS2019", "Dataset",
                       f"Colon mucosa proteomics, untreated UC vs healthy ({STUDY})")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "study": STUDY, "contrast": "inflamed UC at diagnosis vs healthy colon mucosa",
                   "log2_ratio": round(float(r.log2_ratio), 3), "p": float(f"{r.p:.3g}"), "padj": float(f"{r.padj:.3g}"),
                   "significant_by_authors": bool(r.significant_authors), "samples_quantified": int(r.n_samples),
                   "rna_meta_log2fc": None if pd.isna(r.meta_log2fc) else round(float(r.meta_log2fc), 3),
                   "replicates_rna_change": r.replicates or None}
            ev = store.add_evidence(evidence_id("proteomics_de", r.gene), "proteomics_de", "DATASET:PROTEOMICS_SCHNIERS2019",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "protein_changed_in", "DATASET:PROTEOMICS_SCHNIERS2019", ev)
            n += 1
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "proteomics.parquet")
    print(d[d.gene.isin(["S100A8", "S100A9", "LCN2", "DUOX2", "CHI3L1", "MMP1", "IDO1", "CXCL8", "AQP8", "HMGCS2",
                         "SLC26A2", "CXCR2", "REG4", "CD274", "OSMR", "CHRDL2", "XKR9"])][
        ["gene", "log2_ratio", "padj", "meta_log2fc", "replicates"]].to_string())
