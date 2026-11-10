"""CosMx single-cell spatial data (GSE312415; added 2026-11-10, third round): a second imaging platform at single-cell
resolution, with a 6,000-gene panel that covers candidates missing from the Xenium 5K panel.

Sections: each slide holds two rectum biopsies and the cell files do not say which cell belongs to which biopsy, so only
the slides whose two biopsies have the same state are used: the 6K-panel slide of two inflamed UC biopsies (I0380,
I0387) and a 1K-panel slide of two inflamed UC biopsies (I0286, I0277). Slides that mix inflamed and non-inflamed or
healthy tissue are downloaded but not used (assigning cells to a biopsy from their expression would be circular).
Result (2026-11-10): about a third of the cells can be typed, and most candidates are expressed by 5-15 % of the cells of
every lineage alike (background); a gene counts as detected above background only when its best lineage expresses it
in at least twice the cells of the median other lineage (cosmx_detection.parquet). LCN2 and other epithelial genes pass;
the neutrophil and stromal candidates do not, so CosMx does not resolve them in these biopsies.

Per slide: cells typed against the single-cell reference profiles (cell_types below; candidates are not markers),
then with the Xenium functions (spatial.py): share of each candidate's transcripts per lineage, and lineages
enriched among the 10 nearest neighbours of expressing cells (permutation test). Negative-probe and system-control
columns are dropped; cells with fewer than 20 counts are left out.
Outputs: data/processed/cosmx_attribution.parquet, cosmx_neighbourhood.parquet, cosmx_sections.parquet.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from .paths import PROCESSED, RAW
from .spatial import attribution, neighbourhood

D = RAW / "spatial" / "GSE312415"
USED = {"Cos1380R387RRNA6K4u2716": ("6k", "inflamed", "I0380+I0387"),
        "I0286RI0277R1119": ("1k", "inflamed", "I0286+I0277")}
MIN_COUNTS = 20


R_MIN, MARGIN = 0.1, 0.02
SPECIFIC = 2.0   # best lineage expresses the gene in >= 2x the cells of the median other lineage
# reference cell groups folded into the Xenium lineage names, so attribution tables are comparable across platforms
TO_LINEAGE = {"Epithelial": "epithelial", "Fibroblast": "fibroblast", "Inflammatory_fibroblast": "fibroblast",
              "Smooth_muscle": "smooth_muscle", "Endothelial": "endothelial", "Myeloid": "myeloid",
              "Neutrophil": "neutrophil", "T_NK": "T_NK", "B": "B", "Plasma": "plasma", "Mast": "mast",
              "Eosinophil": "myeloid"}


def cell_types(m: sparse.csr_matrix, genes: list[str], exclude: set[str]) -> np.ndarray:
    """Reference-based typing. Marker-score rules (the Xenium rule, or z-scored markers) gave implausible CosMx
    lineage sizes because of background counts. Each cell is instead compared with the 12 single-cell reference
    profiles of the Visium deconvolution (spatial_deconv.reference; up to 40 markers per group present on the panel,
    candidate genes excluded): Pearson correlation of log1p counts per 1,000 with log1p reference CPM over the
    markers; the cell takes the best group if r >= R_MIN and it beats the second by >= MARGIN, else unassigned."""
    from .spatial_deconv import reference, select_markers
    ref = reference()
    mk = select_markers(ref, exclude, set(genes))
    allm = list(dict.fromkeys(g for v in mk.values() for g in v))
    idx = {g: i for i, g in enumerate(genes)}
    lib = np.asarray(m.sum(axis=1)).ravel()
    x = np.log1p((sparse.diags(1e3 / np.maximum(lib, 1)) @ m[:, [idx[g] for g in allm]]).toarray())
    r = np.log1p(ref.loc[allm].to_numpy())
    xc = x - x.mean(axis=1, keepdims=True)
    rc = r - r.mean(axis=0, keepdims=True)
    cor = (xc @ rc) / np.maximum(np.linalg.norm(xc, axis=1, keepdims=True) * np.linalg.norm(rc, axis=0, keepdims=True), 1e-9)
    order = np.sort(cor, axis=1)
    groups = np.array([TO_LINEAGE[c] for c in ref.columns], dtype=object)
    lab = groups[cor.argmax(axis=1)]
    lab[(order[:, -1] < R_MIN) | (order[:, -1] - order[:, -2] < MARGIN)] = "unassigned"
    return lab


def load_slide(stem: str) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame]:
    parts, genes, ids = [], None, []
    for chunk in pd.read_csv(D / f"{stem}_exprMat_file.csv", chunksize=20000):
        if genes is None:
            genes = [c for c in chunk.columns if c not in ("fov", "cell_ID", "cell")
                     and not c.lower().startswith(("negative", "systemcontrol", "negprb", "falsecode"))]
        ids.append(chunk[["fov", "cell_ID"]])
        parts.append(sparse.csr_matrix(chunk[genes].to_numpy(dtype=np.float32)))
    m = sparse.vstack(parts).tocsr()
    idx = pd.concat(ids, ignore_index=True)
    meta = pd.read_csv(D / f"{stem}_metadata_file.csv", usecols=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"])
    cells = idx.merge(meta, on=["fov", "cell_ID"], how="left")
    keep = (np.asarray(m.sum(axis=1)).ravel() >= MIN_COUNTS) & (idx["cell_ID"].to_numpy() != 0) & cells["CenterX_global_px"].notna().to_numpy()
    return m[keep], genes, cells[keep].reset_index(drop=True)


def run(targets: list[str]) -> dict:
    att, nbh, sec = [], [], []
    for stem, (panel, state, biopsies) in USED.items():
        m, genes, cells = load_slide(stem)
        lab = cell_types(m, genes, set(targets))
        xy = cells[["CenterX_global_px", "CenterY_global_px"]].to_numpy(float)
        sec.append({"slide": stem, "panel": panel, "state": state, "biopsies": biopsies, "cells": int(m.shape[0]),
                    "genes": len(genes), "median_counts_per_cell": float(np.median(np.asarray(m.sum(axis=1)).ravel())),
                    **{f"n_{k}": int(v) for k, v in zip(*np.unique(lab, return_counts=True))}})
        att.append(attribution(m, genes, lab, targets).assign(slide=stem, panel=panel, state=state))
        nbh.append(neighbourhood(m, genes, lab, xy, targets).assign(slide=stem, panel=panel, state=state))
        print(stem, m.shape, flush=True)
    a, n, s = pd.concat(att), pd.concat(nbh), pd.DataFrame(sec)
    s["assigned_fraction"] = 1 - s.get("n_unassigned", 0) / s["cells"]
    # above background: among typed cells, % expressing in the best lineage over the median % of the other lineages
    t = a[a.lineage != "unassigned"]
    det = []
    for (gene, panel), d in t.groupby(["gene", "panel"]):
        d = d.sort_values("pct_expressing", ascending=False)
        spec = float(d.pct_expressing.iloc[0] / max(d.pct_expressing.iloc[1:].median(), 0.1))
        det.append({"gene": gene, "panel": panel, "best_lineage": d.lineage.iloc[0],
                    "pct_best": round(float(d.pct_expressing.iloc[0]), 1),
                    "pct_others_median": round(float(d.pct_expressing.iloc[1:].median()), 1),
                    "specificity": round(spec, 2), "above_background": spec >= SPECIFIC})
    pd.DataFrame(det).to_parquet(PROCESSED / "cosmx_detection.parquet", index=False)
    a.to_parquet(PROCESSED / "cosmx_attribution.parquet", index=False)
    n.to_parquet(PROCESSED / "cosmx_neighbourhood.parquet", index=False)
    s.to_parquet(PROCESSED / "cosmx_sections.parquet", index=False)
    return {"slides": len(s), "cells": int(s.cells.sum()), "genes_attributed": int(a["gene"].nunique())}


def load(store=None) -> int:
    """One spatial_cosmx record per gene: lineage shares of its transcripts and enriched neighbours, per panel."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    a = pd.read_parquet(PROCESSED / "cosmx_attribution.parquet")
    nb = pd.read_parquet(PROCESSED / "cosmx_neighbourhood.parquet")
    det = pd.read_parquet(PROCESSED / "cosmx_detection.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE312415", "Dataset", "CosMx 6K and 1K panels, inflamed UC rectum biopsies (GSE312415)")
        for gene, d in a.groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            shares = {p: dd.groupby("lineage")["share"].mean().round(3).nlargest(3).to_dict() for p, dd in d.groupby("panel")}
            nn = nb[(nb.gene == gene) & (nb.p_enriched < 0.05) & (nb.enrichment > 1.2)]
            dd = det[det.gene == gene]
            rec = {"gene": gene, "platform": "CosMx (single cell)", "state": "inflamed UC rectum",
                   "detection": [{"panel": r.panel, "best_lineage": r.best_lineage, "pct_cells_best": r.pct_best,
                                  "pct_cells_other_median": r.pct_others_median, "specificity": r.specificity,
                                  "above_background": bool(r.above_background)} for r in dd.itertuples()],
                   "lineage_share_by_panel": shares,
                   "enriched_neighbours": [{"lineage": r.neighbour_lineage, "enrichment": r.enrichment, "p": r.p_enriched,
                                            "panel": r.panel, "n_cells": r.n_expressing_cells} for r in nn.itertuples()]}
            ev = store.add_evidence(evidence_id("spatial_cosmx", gene), "spatial_cosmx", "DATASET:GSE312415",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "observed_in_situ_in", "DATASET:GSE312415", ev)
            n += 1
    return n


if __name__ == "__main__":
    import json
    from .cli import expression_top
    from .spatial_deconv import EXTRA
    extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
    print(json.dumps(run(list(dict.fromkeys(expression_top(300) + extra + EXTRA))), indent=1))
    a = pd.read_parquet(PROCESSED / "cosmx_attribution.parquet")
    top = a[a.panel == "6k"].sort_values("share").groupby("gene").tail(1).set_index("gene")
    for g in ["S100A8", "CXCL8", "CXCR2", "LCN2", "DUOX2", "MMP1", "MMP3", "CHI3L1", "IDO1", "IL11", "OSMR", "OSM", "CD274"]:
        print(g, top.loc[g, ["lineage", "share"]].to_dict() if g in top.index else "not on panel")
