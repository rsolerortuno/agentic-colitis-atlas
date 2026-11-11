"""M11: does the cell-type attribution of Smillie 2019 hold in a second dataset (Garrido-Trigo et al. 2023)?

The second dataset has five coarse lineages, so the comparison is made at lineage level: Smillie clusters are mapped
to epithelial / stromal / T and NK / plasma / myeloid (B cells, absent from the second dataset, are left out). For
each gene the dominant lineage (largest share of the gene's transcripts in UC mucosa) is compared between datasets.
Output: data/processed/sc_replication.parquet.

Four atlases (added 2026-11-11, run4): the same comparison across Smillie 2019, Garrido-Trigo 2023, Gubatan 2026
(SCP3755, with neutrophils) and GSE282122 (Thomas 2024, full count matrices typed in sc_gse282122.py). All are folded
into the five lineages the second atlas has: epithelial, stromal (fibroblasts, endothelium, pericytes, glia, smooth
muscle), T/NK, plasma and myeloid (macrophages, dendritic cells, monocytes, mast cells, and neutrophils and
eosinophils where an atlas has them); B cells are left out, as before. Per gene: the dominant lineage in inflamed UC
mucosa in each atlas, and in how many of the other three the Smillie lineage is also dominant (an atlas that does not
measure the gene counts as not agreeing). sc_replication.parquet is kept unchanged (the sc_summary records cite it).
Outputs: data/processed/sc_replication4.parquet, sc_replication4_summary.json.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

LINEAGE_2023 = {"colon epithelial cell": "epithelial", "stromal cell of lamina propria of colon": "stromal",
                "T cell of anorectum": "T/NK", "plasma cell": "plasma", "myeloid cell": "myeloid"}
MYELOID = ("Macrophages", "Inflammatory Monocytes", "DC1", "DC2", "CD69+ Mast", "CD69- Mast", "Cycling Monocytes")
B_CELLS = ("Follicular", "GC", "Cycling B")


def lineage_2019(cluster: str, compartment: str) -> str | None:
    if compartment == "epithelial":
        return "epithelial"
    if compartment == "stromal":
        return "stromal"
    if cluster == "Plasma":
        return "plasma"
    if cluster in MYELOID:
        return "myeloid"
    if cluster in B_CELLS:
        return None
    return "T/NK"


def shares_2019(health: str = "inflamed") -> pd.DataFrame:
    s = pd.read_parquet(PROCESSED / "sc_celltype_share.parquet")
    s = s[s["health"] == health].copy()
    s["lineage"] = [lineage_2019(c, k) for c, k in zip(s["cluster"], s["compartment"])]
    s = s.dropna(subset=["lineage"]).groupby(["gene", "lineage"], as_index=False)["cp10k_sum"].sum()
    s["share"] = s["cp10k_sum"] / s.groupby("gene")["cp10k_sum"].transform("sum")
    return s


def shares_2023(genes: list[str]) -> pd.DataFrame:
    import anndata as ad
    from .expression import resolve_symbol
    a = ad.read_h5ad(RAW / "sc_garrido2023" / "ibd_hc_10x.h5ad")
    a = a[a.obs["disease"] == "ulcerative colitis"]
    sym = a.var["feature_name"].astype(str).map(lambda g: resolve_symbol(g) or g)
    keep = np.where(sym.isin(genes))[0]
    x = a.X[:, keep]
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    if x.max() < 50:  # log-normalised (CELLxGENE convention): back to the normalised scale before summing
        x = np.expm1(x)
    lin = a.obs["cell_type"].astype(str).map(LINEAGE_2023).to_numpy()
    rows = []
    for L in sorted(set(lin)):
        tot = x[lin == L].sum(axis=0)
        rows.append(pd.DataFrame({"gene": sym.to_numpy()[keep], "lineage": L, "sum": tot}))
    s = pd.concat(rows).groupby(["gene", "lineage"], as_index=False)["sum"].sum()
    s["share"] = s["sum"] / s.groupby("gene")["sum"].transform("sum")
    return s


LINEAGE_GUBATAN = {"Enterocytes": "epithelial", "Goblet": "epithelial", "Stem Cycling": "epithelial",
                   "Enteroendocrine": "epithelial", "Tuft": "epithelial", "Stromal": "stromal", "Endothelial": "stromal",
                   "Pericyte": "stromal", "Glial": "stromal", "T Cell and NK": "T/NK", "Plasma": "plasma",
                   "Macrophage and DC": "myeloid", "Mast Cell": "myeloid", "Neutrophil": "myeloid",
                   "Eosinophil": "myeloid"}   # B Cell and Erythrocyte left out
LINEAGE_GSE282122 = {"epithelial": "epithelial", "stromal": "stromal", "endothelial": "stromal", "T/NK": "T/NK",
                     "plasma": "plasma", "myeloid": "myeloid", "mast": "myeloid", "neutrophil": "myeloid"}   # B left out
ATLASES = {"smillie": "Smillie 2019", "garrido": "Garrido-Trigo 2023", "gubatan": "Gubatan 2026",
           "gse282122": "GSE282122 (Thomas 2024)"}


def _fold(s: pd.DataFrame, mapping: dict, col: str) -> pd.DataFrame:
    s = s.assign(lineage=s[col].map(mapping)).dropna(subset=["lineage"])
    s = s.groupby(["gene", "lineage"], as_index=False)["cp10k_sum"].sum()
    s["share"] = s["cp10k_sum"] / s.groupby("gene")["cp10k_sum"].transform("sum")
    return s


def shares_gubatan(genes: list[str]) -> pd.DataFrame:
    s = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    return _fold(s[(s["health"] == "inflamed") & s["gene"].isin(genes)], LINEAGE_GUBATAN, "compartment")


def shares_gse282122(genes: list[str]) -> pd.DataFrame:
    s = pd.read_parquet(PROCESSED / "sc4_celltype_share.parquet")
    return _fold(s[(s["health"] == "inflamed") & s["gene"].isin(genes)], LINEAGE_GSE282122, "cluster")


def run4(genes: list[str]) -> dict:
    import json
    top = lambda s: s.sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    a = shares_2019()
    full = {"smillie": a[a["gene"].isin(genes)], "garrido": shares_2023(genes), "gubatan": shares_gubatan(genes),
            "gse282122": shares_gse282122(genes)}
    tabs = {k: top(v) for k, v in full.items()}
    out = pd.DataFrame({"gene": tabs["smillie"].index})
    for k, t in tabs.items():
        out[f"lineage_{k}"] = out.gene.map(t["lineage"])
        out[f"share_{k}"] = out.gene.map(t["share"]).astype(float).round(3)
    for k in ("garrido", "gubatan", "gse282122"):   # share of the Smillie lineage in each other atlas (for the figure)
        s = full[k].set_index(["gene", "lineage"])["share"]
        out[f"share_same_{k}"] = [round(float(s.get((g, L), 0.0)), 3) if pd.notna(out.at[i, f"lineage_{k}"]) else None
                                  for i, (g, L) in enumerate(zip(out.gene, out.lineage_smillie))]
    others = [k for k in ATLASES if k != "smillie"]
    for k in others:
        out[f"agree_{k}"] = out[f"lineage_{k}"] == out["lineage_smillie"]
    out["n_measured"] = out[[f"lineage_{k}" for k in others]].notna().sum(axis=1).astype(int)
    out["n_agree"] = out[[f"agree_{k}" for k in others]].sum(axis=1).astype(int)
    out["replicated"] = out["n_agree"] >= 2   # the Smillie lineage is dominant in at least 2 of the 3 other atlases
    out.to_parquet(PROCESSED / "sc_replication4.parquet", index=False)
    s = {"genes": int(len(out)), "atlases": ATLASES,
         "agree_count": {str(i): int((out.n_agree == i).sum()) for i in range(4)},
         "replicated_2_of_3": int(out.replicated.sum()),
         "agree_with_smillie": {k: {"agree": int(out[f"agree_{k}"].sum()), "measured": int(out[f"lineage_{k}"].notna().sum())}
                                for k in others},
         "all_four_measured": int((out.n_measured == 3).sum())}
    (PROCESSED / "sc_replication4_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def genes4() -> list[str]:
    """Genes compared across the four atlases: the 500 largest meta-analysis changes and every gene of the discovery
    table (the shortlist genes are among them)."""
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    top = m.reindex(m["meta_log2fc"].abs().sort_values(ascending=False).index).head(500)["gene"].tolist()
    # discovery.parquet (discovery step), not discovery_scored.parquet: the latter is written later in run-all by the
    # report layers, and reading it here gave 390 genes instead of 1,003 in the first clean rebuild (2026-11-11)
    p = PROCESSED / "discovery.parquet"
    disc = pd.read_parquet(p)["gene"].tolist() if p.exists() else []
    return list(dict.fromkeys(top + disc))


def load4(store=None) -> dict:
    """One sc_replication_4atlas record per compared gene (a new evidence type: existing records do not change)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "sc_replication4.parquet")
    n = 0
    with store.batch():
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "comparison": "dominant lineage in inflamed UC mucosa, five shared lineages",
                   **{f"lineage_{k}": getattr(r, f"lineage_{k}") for k in ATLASES},
                   **{f"share_{k}": (None if pd.isna(getattr(r, f"share_{k}")) else float(getattr(r, f"share_{k}")))
                      for k in ATLASES},
                   "other_atlases_agreeing": int(r.n_agree), "other_atlases_measuring": int(r.n_measured),
                   "replicated_2_of_3": bool(r.replicated)}
            ev = store.add_evidence(evidence_id("sc_replication_4atlas", r.gene), "sc_replication_4atlas",
                                    "DATASET:SCP259", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "cell_attribution_checked_in", "DATASET:GSE282122", ev)
            n += 1
    return {"sc_replication_4atlas_records": n}


def run(genes: list[str]) -> dict:
    a = shares_2019()
    a = a[a["gene"].isin(genes)]
    b = shares_2023(genes)
    top = lambda s: s.sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    ta, tb = top(a), top(b)
    common = ta.index.intersection(tb.index)
    out = pd.DataFrame({"gene": common, "lineage_2019": ta.loc[common, "lineage"].values,
                        "share_2019": ta.loc[common, "share"].values.round(3),
                        "lineage_2023": tb.loc[common, "lineage"].values, "share_2023": tb.loc[common, "share"].values.round(3)})
    out["same_lineage"] = out["lineage_2019"] == out["lineage_2023"]
    out.to_parquet(PROCESSED / "sc_replication.parquet", index=False)
    return {"genes_compared": len(out), "same_dominant_lineage": int(out["same_lineage"].sum()),
            "share_same": round(float(out["same_lineage"].mean()), 3)}


if __name__ == "__main__":
    import json
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    m = m.reindex(m["meta_log2fc"].abs().sort_values(ascending=False).index).head(500)
    print(json.dumps(run(m["gene"].tolist()), indent=1))
