"""Covariate-adjusted discovery contrast (added 2026-11-10, fourth round).

The MSCCR contrast (expression.py) compares patient means of inflamed UC and control biopsies. Controls are older
(median 56 vs 44 years) and come mostly from the rectum and caecum, inflamed biopsies more often from the sigmoid and
left colon, and colon segments differ in expression. Here, per gene, an ordinary least-squares model on biopsy-level
log2 normalised counts:

    expression ~ inflamed UC (vs control) + colon segment (6 levels) + age + sex

with standard errors clustered by patient (several biopsies per patient; CR0 sandwich), BH across genes. The
adjusted log2FC is compared with the unadjusted one; a candidate "holds" when the adjusted change keeps its sign, at
least half of its size and padj < 0.05.
Outputs: data/processed/de_adjusted.parquet, de_adjusted_summary.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED


def run() -> dict:
    x = pd.read_parquet(PROCESSED / "expr_GSE193677.parquet")
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    s = s[(s.cohort == "GSE193677") & (s.timepoint == "baseline") & s.group.isin(["UC_inflamed", "healthy"])
          & s["sample"].isin(x.columns)].dropna(subset=["region", "age", "sex"])
    y = x[s["sample"]].to_numpy(np.float64).T                       # biopsies x genes
    reg = pd.get_dummies(s["region"], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(s)), (s.group == "UC_inflamed").to_numpy(float), reg.to_numpy(),
                         s["age"].to_numpy(float) / 10, (s["sex"] == "Male").to_numpy(float)])
    xtx_inv = np.linalg.pinv(X.T @ X)
    B = xtx_inv @ X.T @ y
    E = y - X @ B
    h = (xtx_inv @ X.T)[1]                                          # row giving the inflamed coefficient
    clusters = pd.factorize(s["patient"])[0]
    C = np.zeros((clusters.max() + 1, len(s)))
    C[clusters, np.arange(len(s))] = 1
    S = C @ (h[:, None] * E)                                        # per-patient score contributions
    se = np.sqrt((S ** 2).sum(axis=0))
    beta = B[1]
    dof = clusters.max()
    p = 2 * stats.t.sf(np.abs(beta / np.maximum(se, 1e-12)), dof)
    out = pd.DataFrame({"gene": x.index, "adjusted_log2fc": beta, "adjusted_se": se, "adjusted_p": p})
    out["adjusted_padj"] = multipletests(out.adjusted_p.fillna(1), method="fdr_bh")[1]
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    un = de[(de.cohort == "GSE193677") & (de.contrast == "inflamed_vs_control")].set_index("gene")
    out["unadjusted_log2fc"] = out.gene.map(un.log2fc)
    out["unadjusted_padj"] = out.gene.map(un.padj)
    out["holds"] = ((np.sign(out.adjusted_log2fc) == np.sign(out.unadjusted_log2fc))
                    & (out.adjusted_log2fc.abs() >= 0.5 * out.unadjusted_log2fc.abs()) & (out.adjusted_padj < 0.05))
    out.to_parquet(PROCESSED / "de_adjusted.parquet", index=False)
    sig = out[(out.unadjusted_padj < 0.05) & (out.unadjusted_log2fc.abs() >= 1)]
    summ = {"biopsies": int(len(s)), "patients": int(clusters.max() + 1),
            "covariates": ["colon segment", "age (per 10 years)", "sex"],
            "spearman_adjusted_vs_unadjusted": round(float(stats.spearmanr(out.adjusted_log2fc, out.unadjusted_log2fc,
                                                                            nan_policy="omit").statistic), 3),
            "changes_tested": int(len(sig)), "changes_holding": int(sig.holds.sum()),
            "share_holding": round(float(sig.holds.mean()), 3)}
    (PROCESSED / "de_adjusted_summary.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    return summ


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "de_adjusted.parquet")
    keep = set(pd.read_parquet(PROCESSED / "meta.parquet").query("meta_padj < 0.05")["gene"])
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE193677", "Dataset", "MSCCR biopsy RNA-seq (GSE193677)")
        for r in d[d.gene.isin(keep)].itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "cohort": "GSE193677", "contrast": "inflamed UC vs control, biopsy level",
                   "covariates": "colon segment, age, sex; SE clustered by patient",
                   "adjusted_log2fc": round(float(r.adjusted_log2fc), 3), "adjusted_padj": float(f"{r.adjusted_padj:.3g}"),
                   "unadjusted_log2fc": None if pd.isna(r.unadjusted_log2fc) else round(float(r.unadjusted_log2fc), 3),
                   "holds_after_adjustment": bool(r.holds)}
            ev = store.add_evidence(evidence_id("expression_adjusted", r.gene), "expression_adjusted", "DATASET:GSE193677",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "changed_in_adjusted_model", "DATASET:GSE193677", ev)
            n += 1
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "de_adjusted.parquet").set_index("gene")
    print(d.loc[[g for g in ["CXCR2", "CXCL8", "S100A8", "LCN2", "DUOX2", "MMP1", "CHI3L1", "IDO1", "AQP8", "HMGCS2",
                             "SLC26A2", "CD274", "OSMR", "CHRDL2", "LAX1", "SLC26A3"] if g in d.index],
                ["unadjusted_log2fc", "adjusted_log2fc", "adjusted_padj", "holds"]].round(3).to_string())
