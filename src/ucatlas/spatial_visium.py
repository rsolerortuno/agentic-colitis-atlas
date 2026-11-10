"""Whole-transcriptome spatial context (10x Visium, GSE189184; added 2026-11-10).

The Xenium 5K panel lacks several candidates (for example CXCL8, S100A8, CHI3L1, LCN2, MMP3). Visium measures every
gene, at the price of resolution: one 55 um spot holds several cells. Sections: seven inflamed UC colon sections (five
biopsies, two resections) and two healthy colon biopsies. Checkpoint-inhibitor colitis sections are not used.

Per gene:
1. In situ change: per section, counts per million over all tissue spots; UC inflamed vs healthy as log2 ratio of the
   mean CPM and a Welch t-test over sections (7 vs 2, so a descriptive check, not a discovery test).
2. Compartment by co-variation: in each UC section, the Spearman correlation across spots between the gene
   (log CPM) and each lineage score (mean of the lineage markers, scaled by their mean; the same marker sets as the
   Xenium analysis, never including the gene itself). The lineage with the highest median correlation over the UC
   sections is the gene's spatial compartment; the record keeps every lineage's median correlation and the number of
   sections where that lineage ranks first.
Outputs: data/processed/spatial_visium.parquet (one row per gene) and spatial_visium_sections.parquet.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse, stats

from .paths import PROCESSED, RAW
from .spatial import MARKERS, read_10x_h5
from .spatial_download import VISIUM

D = RAW / "spatial" / "GSE189184"
MIN_SPOTS_EXPRESSING = 20


def load(gsm: str) -> tuple[sparse.csr_matrix, list[str]]:
    from .expression import resolve_symbol
    m, names, bcs, types = read_10x_h5(D / f"{gsm}_filtered_feature_bc_matrix.h5")
    pos = pd.read_csv(D / f"{gsm}_tissue_positions_list.csv", header=None, index_col=0)
    keep = np.asarray(pos.reindex(bcs)[1].fillna(0).astype(int) == 1)
    m = m[keep][:, types == "Gene Expression"]
    genes = [resolve_symbol(g) or g for g, t in zip(names, types) if t == "Gene Expression"]
    # one column per gene symbol
    s = pd.Series(range(len(genes)), index=genes)
    if s.index.has_duplicates:
        first = ~s.index.duplicated()
        m, genes = m[:, np.flatnonzero(first)], list(s.index[first])
    return m.tocsr(), genes


def _lineage_scores(logcpm: sparse.csr_matrix, idx: dict[str, int], exclude: str | None = None) -> dict[str, np.ndarray]:
    out = {}
    for lin, mk in MARKERS.items():
        present = [idx[g] for g in mk if g in idx and g != exclude]
        if len(present) < 2:
            continue
        sub = logcpm[:, present]
        mu = np.asarray(sub.mean(axis=0)).ravel()
        out[lin] = np.asarray((sub @ sparse.diags(1 / np.maximum(mu, 1e-9))).sum(axis=1)).ravel() / len(present)
    return out


def run(targets: list[str]) -> dict:
    sections, cpm_by_section, cor_rows = [], {}, []
    for gsm, (stem, state) in VISIUM.items():
        m, genes = load(gsm)
        lib = np.asarray(m.sum(axis=1)).ravel()
        cpm = sparse.diags(1e6 / np.maximum(lib, 1)) @ m
        logcpm = cpm.copy()
        logcpm.data = np.log1p(logcpm.data)
        logcpm = logcpm.tocsc()
        idx = {g: i for i, g in enumerate(genes)}
        sections.append({"gsm": gsm, "sample": stem, "state": state, "spots": int(m.shape[0]),
                         "median_counts_per_spot": float(np.median(lib)), "genes_detected": int((m.getnnz(axis=0) > 0).sum())})
        tot = np.asarray(m.sum(axis=0)).ravel()
        cpm_by_section[gsm] = pd.Series(1e6 * tot / tot.sum(), index=genes)
        if state != "UC_inflamed":
            continue
        base_scores = _lineage_scores(logcpm, idx)
        marker_genes = {g for mk in MARKERS.values() for g in mk}
        for g in targets:
            if g not in idx:
                continue
            col = logcpm[:, idx[g]].toarray().ravel()
            if (col > 0).sum() < MIN_SPOTS_EXPRESSING:
                continue
            scores = _lineage_scores(logcpm, idx, exclude=g) if g in marker_genes else base_scores
            for lin, sc in scores.items():
                rho = stats.spearmanr(col, sc).statistic
                cor_rows.append({"gene": g, "gsm": gsm, "lineage": lin, "rho": float(rho),
                                 "spots_expressing": int((col > 0).sum()), "spots": int(len(col))})
        print("visium", gsm, state, m.shape, flush=True)

    sec = pd.DataFrame(sections)
    cpm = pd.DataFrame(cpm_by_section).fillna(0)
    uc = [g for g, (_, st) in VISIUM.items() if st == "UC_inflamed"]
    hc = [g for g, (_, st) in VISIUM.items() if st == "healthy"]
    cor = pd.DataFrame(cor_rows)
    rows = []
    for g, d in cor.groupby("gene"):
        med = d.groupby("lineage")["rho"].median().sort_values(ascending=False)
        firsts = d.loc[d.groupby("gsm")["rho"].idxmax(), "lineage"].value_counts()
        a, b = np.log2(cpm.loc[g, uc] + 1), np.log2(cpm.loc[g, hc] + 1)
        t = stats.ttest_ind(a, b, equal_var=False)
        rows.append({"gene": g, "compartment": med.index[0], "compartment_rho": round(float(med.iloc[0]), 3),
                     "second_compartment": med.index[1], "second_rho": round(float(med.iloc[1]), 3),
                     "sections_compartment_first": int(firsts.get(med.index[0], 0)), "uc_sections": int(d["gsm"].nunique()),
                     "lineage_rho": {k: round(float(v), 3) for k, v in med.items()},
                     "mean_cpm_uc": round(float(cpm.loc[g, uc].mean()), 2), "mean_cpm_healthy": round(float(cpm.loc[g, hc].mean()), 2),
                     "log2_ratio_uc_vs_healthy": round(float(np.log2((cpm.loc[g, uc].mean() + 1) / (cpm.loc[g, hc].mean() + 1))), 3),
                     "p_sections": float(t.pvalue) if np.isfinite(t.pvalue) else None,
                     "median_spots_expressing_pct": round(float((d.groupby("gsm").first()["spots_expressing"]
                                                                 / d.groupby("gsm").first()["spots"]).median() * 100), 1)})
    out = pd.DataFrame(rows)
    out["lineage_rho"] = out["lineage_rho"].map(lambda x: __import__("json").dumps(x))
    out.to_parquet(PROCESSED / "spatial_visium.parquet", index=False)
    sec.to_parquet(PROCESSED / "spatial_visium_sections.parquet", index=False)
    return {"sections": len(sec), "genes": len(out)}


def load_graph(store=None) -> int:
    """One spatial_visium evidence record per gene."""
    import json
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    v = pd.read_parquet(PROCESSED / "spatial_visium.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE189184", "Dataset", "10x Visium whole transcriptome, UC inflamed and healthy colon (GSE189184)")
        for r in v.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "platform": "10x Visium (55 um spots, several cells per spot)",
                   "compartment": r.compartment, "compartment_rho": r.compartment_rho,
                   "second_compartment": r.second_compartment, "second_rho": r.second_rho,
                   "sections_compartment_first": r.sections_compartment_first, "uc_sections": r.uc_sections,
                   "lineage_rho": json.loads(r.lineage_rho), "mean_cpm_uc": r.mean_cpm_uc,
                   "mean_cpm_healthy": r.mean_cpm_healthy, "log2_ratio_uc_vs_healthy": r.log2_ratio_uc_vs_healthy,
                   "p_sections": None if r.p_sections is None or pd.isna(r.p_sections) else float(f"{r.p_sections:.3g}"),
                   "median_spots_expressing_pct": r.median_spots_expressing_pct,
                   "sections": "7 UC inflamed, 2 healthy"}
            ev = store.add_evidence(evidence_id("spatial_visium", r.gene), "spatial_visium", "DATASET:GSE189184",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "observed_in_situ_in", "DATASET:GSE189184", ev)
            n += 1
    return n


if __name__ == "__main__":
    import json
    from .cli import expression_top
    extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
    print(json.dumps(run(list(dict.fromkeys(expression_top(300) + extra))), indent=1))
    print(pd.read_parquet(PROCESSED / "spatial_visium.parquet").drop(columns="lineage_rho").head(30).to_string())
