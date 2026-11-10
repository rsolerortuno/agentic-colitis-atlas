"""M11 spatial: what surrounds each candidate in UC colon (Xenium 5K panel).

Sections: GSE312420 (healthy, non-inflamed UC and two inflamed UC rectum sections) and, added 2026-11-10 because two
GSE312420 sections had too few transcripts per cell to type most cells, GSE339323 (Xenium Prime 5K, one UC resection,
most and least inflamed region; median about 300 transcripts per cell). Sections where fewer than half the cells
can be typed are still analysed but left out of the knowledge graph (named in each record). Genes that are not on
the 5K panel are covered by the whole-transcriptome Visium sections (spatial_visium.py). Images are not used.

1. Cell typing without the candidate genes: each cell gets the lineage whose marker set (from the single-cell
   datasets, restricted to the panel) has the largest normalised sum; cells with no marker counts stay 'unassigned'.
   Markers never include the gene being tested (the marker sets are fixed before any candidate is examined, and a
   candidate that is a marker is reported with that flag).
2. Attribution in situ: share of a gene's transcripts in each lineage, per section.
3. Neighbourhood: for the cells expressing a gene (>= 1 count), the lineage mix of their 10 nearest neighbours,
   compared with the same number of random cells of the section (1,000 permutations, seed fixed per gene). When
   more than 5,000 cells express a gene, a random subset of 5,000 of them is tested (same seed).
Outputs: data/processed/spatial_attribution.parquet, spatial_neighbourhood.parquet, spatial_sections.parquet.
"""
from __future__ import annotations

import gzip
import io
import zlib

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.spatial import cKDTree

from .paths import PROCESSED, RAW

D = RAW / "spatial" / "GSE312420"
MARKERS = {  # lineage markers commonly present on Xenium panels (checked against each panel at run time)
    "epithelial": ["EPCAM", "KRT8", "KRT20", "CDH1", "PIGR", "CEACAM5", "AGR2", "MUC2", "LGALS4", "FABP1"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "ADAMDEC1", "CXCL14"],
    "endothelial": ["PECAM1", "VWF", "PLVAP", "CDH5", "RAMP2", "ACKR1"],
    "smooth_muscle": ["ACTA2", "DES", "MYH11", "CNN1"],
    "T_NK": ["CD3E", "CD3D", "CD2", "TRAC", "IL7R", "NKG7", "GZMA", "CD8A"],
    "B": ["MS4A1", "CD79A", "CD19", "BANK1", "PAX5"],
    "plasma": ["JCHAIN", "MZB1", "IGHA1", "XBP1", "TNFRSF17", "SDC1"],
    "myeloid": ["CD68", "CD163", "C1QA", "C1QB", "LYZ", "CD14", "FCGR3A", "MS4A7"],
    # CXCR1/CXCR2/S100A8/S100A9 are candidates, so they are not used to define neutrophils
    "neutrophil": ["FCGR3B", "CSF3R", "MNDA", "FPR1", "CMTM2", "PROK2", "MMP25"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
}
SECTIONS_5K = {"GSM9345505": "healthy", "GSM9345506": "noninflamed", "GSM9345507": "inflamed", "GSM9345508": "inflamed",
               "GSM9893496": "inflamed", "GSM9893497": "noninflamed"}
D2 = RAW / "spatial" / "GSE339323"
STUDY = {"GSM9345505": "GSE312420", "GSM9345506": "GSE312420", "GSM9345507": "GSE312420", "GSM9345508": "GSE312420",
         "GSM9893496": "GSE339323", "GSM9893497": "GSE339323"}


def read_10x_h5(path) -> tuple[sparse.csr_matrix, list[str], list[str], np.ndarray]:
    """10x Genomics h5 (Xenium or Visium): cells x features CSR, feature names, barcodes, feature types."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(g["shape"][:])
        m = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).T.tocsr()
        names = [x.decode() for x in g["features/name"][:]]
        types = np.array([x.decode() for x in g["features/feature_type"][:]])
        bcs = [x.decode() for x in g["barcodes"][:]]
    return m, names, bcs, types


def _sample_files(gsm: str) -> dict:
    files = {p.name.split("_")[-1].split(".")[0]: p for p in D.glob(f"{gsm}_*")}
    return {k: files[k] for k in ("matrix", "barcodes", "features", "cells")}


def _read_gz(path) -> bytes:
    raw = path.read_bytes()
    try:
        return gzip.decompress(raw)
    except (OSError, zlib.error):
        return raw


def load_section(gsm: str) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame]:
    """cells x genes counts (CSR), gene names and the cell table (x, y), in the same cell order."""
    from .expression import resolve_symbol
    if STUDY[gsm] == "GSE339323":
        h5 = next(D2.glob(f"{gsm}_*_cell_feature_matrix.h5"))
        m, names, bcs, types = read_10x_h5(h5)
        cells = pd.read_parquet(io.BytesIO(_read_gz(next(D2.glob(f"{gsm}_*_cells.parquet.gz")))))
        cells["cell_id"] = cells["cell_id"].astype(str)
        cells = cells.set_index("cell_id").reindex(bcs)
        is_gene = types == "Gene Expression"
        genes = [resolve_symbol(g) or g for g, k in zip(names, is_gene) if k]
        return m[:, is_gene], genes, cells[["x_centroid", "y_centroid"]].reset_index(drop=True)
    f = _sample_files(gsm)
    m = mmread(io.BytesIO(_read_gz(f["matrix"]))).tocsr().T.tocsr()          # cells x features
    feats = pd.read_csv(io.BytesIO(_read_gz(f["features"])), sep="\t", header=None)
    bc = pd.read_csv(io.BytesIO(_read_gz(f["barcodes"])), header=None)[0].astype(str)
    cells = pd.read_parquet(io.BytesIO(_read_gz(f["cells"])))
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id").reindex(bc)
    is_gene = feats[2].eq("Gene Expression").to_numpy()
    genes = [resolve_symbol(g) or g for g in feats.loc[is_gene, 1]]
    return m[:, is_gene], genes, cells[["x_centroid", "y_centroid"]].reset_index(drop=True)


def cell_types(m: sparse.csr_matrix, genes: list[str]) -> np.ndarray:
    lib = np.asarray(m.sum(axis=1)).ravel()
    norm = sparse.diags(1e2 / np.maximum(lib, 1)) @ m
    idx = {g: i for i, g in enumerate(genes)}
    names, cols = [], []
    for lin, mk in MARKERS.items():
        present = [idx[g] for g in mk if g in idx]
        if len(present) >= 2:
            names.append(lin)
            # scale each marker by its mean so that one very abundant marker does not dominate
            sub = norm[:, present]
            mu = np.asarray(sub.mean(axis=0)).ravel()
            cols.append(np.asarray((sub @ sparse.diags(1 / np.maximum(mu, 1e-9))).sum(axis=1)).ravel() / len(present))
    s = np.column_stack(cols)
    lab = np.array(names, dtype=object)[s.argmax(axis=1)]
    lab[s.max(axis=1) == 0] = "unassigned"
    return lab


def attribution(m, genes, lab, targets) -> pd.DataFrame:
    idx = {g: i for i, g in enumerate(genes)}
    rows = []
    for g in targets:
        if g not in idx:
            continue
        v = np.asarray(m[:, idx[g]].todense()).ravel()
        tot = v.sum()
        if tot == 0:
            continue
        for L in np.unique(lab):
            rows.append({"gene": g, "lineage": L, "share": float(v[lab == L].sum() / tot),
                         "pct_expressing": float(100 * (v[lab == L] > 0).mean())})
    return pd.DataFrame(rows)


MAX_EXPRESSING = 5000  # added 2026-11-10 for the 190,000-280,000-cell GSE339323 sections


def neighbourhood(m, genes, lab, xy, targets, k=10, n_perm=1000) -> pd.DataFrame:
    idx = {g: i for i, g in enumerate(genes)}
    ok = ~np.isnan(xy).any(axis=1)
    tree = cKDTree(xy[ok])
    _, nn = tree.query(xy[ok], k=k + 1)
    nn = nn[:, 1:]
    labs = lab[ok]
    lineages = sorted(set(labs) - {"unassigned"})
    onehot = np.stack([(labs == L).astype(float) for L in lineages], axis=1)
    mix = onehot[nn].mean(axis=1)                       # per cell: lineage mix of its k neighbours
    rows = []
    for g in targets:
        if g not in idx:
            continue
        expr = np.asarray(m[:, idx[g]].todense()).ravel()[ok] > 0
        n_all = int(expr.sum())
        if n_all < 20:
            continue
        rng = np.random.default_rng(zlib.crc32(g.encode()))  # seed fixed per gene (reproducible)
        where = np.flatnonzero(expr)
        if n_all > MAX_EXPRESSING:  # large sections: a random subset of the expressing cells, same test
            where = np.sort(rng.choice(where, MAX_EXPRESSING, replace=False))
        n = len(where)
        obs = mix[where].mean(axis=0)
        null = np.stack([mix[rng.choice(len(mix), n, replace=False)].mean(axis=0) for _ in range(n_perm)])
        for j, L in enumerate(lineages):
            p = (np.sum(null[:, j] >= obs[j]) + 1) / (n_perm + 1)
            rows.append({"gene": g, "neighbour_lineage": L, "observed": round(float(obs[j]), 4),
                         "expected": round(float(null[:, j].mean()), 4),
                         "enrichment": round(float(obs[j] / max(null[:, j].mean(), 1e-9)), 3),
                         "p_enriched": round(float(p), 4), "n_expressing_cells": n_all, "n_tested_cells": n})
    return pd.DataFrame(rows)


def run(targets: list[str]) -> dict:
    att, nbh, sec = [], [], []
    for gsm, state in SECTIONS_5K.items():
        m, genes, cells = load_section(gsm)
        lab = cell_types(m, genes)
        xy = cells[["x_centroid", "y_centroid"]].to_numpy(float)
        sec.append({"gsm": gsm, "study": STUDY[gsm], "state": state, "cells": int(m.shape[0]), "genes": len(genes),
                    "median_counts_per_cell": float(np.median(np.asarray(m.sum(axis=1)).ravel())),
                    **{f"n_{k}": int(v) for k, v in zip(*np.unique(lab, return_counts=True))}})
        att.append(attribution(m, genes, lab, targets).assign(gsm=gsm, state=state))
        if state != "healthy":
            nbh.append(neighbourhood(m, genes, lab, xy, targets).assign(gsm=gsm, state=state))
        print(gsm, state, m.shape, flush=True)
    a, n, s = pd.concat(att), pd.concat(nbh), pd.DataFrame(sec)
    s["assigned_fraction"] = 1 - s["n_unassigned"] / s["cells"]
    a.to_parquet(PROCESSED / "spatial_attribution.parquet", index=False)
    n.to_parquet(PROCESSED / "spatial_neighbourhood.parquet", index=False)
    s.to_parquet(PROCESSED / "spatial_sections.parquet", index=False)
    return {"sections": len(s), "genes_attributed": a["gene"].nunique(), "genes_neighbourhood": n["gene"].nunique()}


if __name__ == "__main__":
    import json
    from .cli import expression_top
    extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
    print(json.dumps(run(list(dict.fromkeys(expression_top(300) + extra))), indent=1))


def load(store=None, min_assigned: float = 0.5) -> int:
    """One spatial_summary record per gene from the sections with >= 50 % typed cells (low-count sections are left
    out and named in the record)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    sec = pd.read_parquet(PROCESSED / "spatial_sections.parquet")
    good = set(sec.loc[sec["assigned_fraction"] >= min_assigned, "gsm"])
    dropped = sorted(set(sec["gsm"]) - good)
    a = pd.read_parquet(PROCESSED / "spatial_attribution.parquet")
    a = a[a["gsm"].isin(good)]
    nb = pd.read_parquet(PROCESSED / "spatial_neighbourhood.parquet")
    nb = nb[nb["gsm"].isin(good)]
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE312420", "Dataset", "Xenium 5K panel, UC and healthy colon (GSE312420 and GSE339323)")
        for gene, d in a.groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            att = {st: dd.groupby("lineage")["share"].mean().round(3).nlargest(3).to_dict() for st, dd in d.groupby("state")}
            nn = nb[(nb["gene"] == gene) & (nb["p_enriched"] < 0.05) & (nb["enrichment"] > 1.2)]
            rec = {"gene": gene, "attribution_by_state": att, "panel": "Xenium 5K", "studies": ["GSE312420", "GSE339323"],
                   "enriched_neighbours": [{"lineage": r.neighbour_lineage, "enrichment": r.enrichment,
                                            "p": r.p_enriched, "section_state": r.state, "n_cells": r.n_expressing_cells}
                                           for r in nn.itertuples()],
                   "sections_used": sorted(good), "sections_left_out_low_counts": dropped}
            ev = store.add_evidence(evidence_id("spatial_summary", gene), "spatial_summary", "DATASET:GSE312420",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "observed_in_situ_in", "DATASET:GSE312420", ev)
            n += 1
    return n
