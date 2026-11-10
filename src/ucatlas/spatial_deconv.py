"""Visium spot deconvolution and in-situ cell attribution (added 2026-11-10; closes the "spots mix cells" limit).

spatial_visium.py tells which compartment a gene co-varies with across spots; one 55 um spot holds several cells, so
that is a compartment, not a cell. This module estimates each spot's cell make-up and then which cell type explains a
gene's expression in situ, with the design of the kidney atlas (Li et al. 2022 Visium):

1. Reference: the colon single-cell atlas with neutrophils (Gubatan et al. 2026, SCP3755, inflamed and non-inflamed
   biopsies), clusters folded into 12 groups a spot can separate (CELL_GROUPS); profile = cell-weighted mean CP10K,
   rescaled to CPM. Markers: per group up to MARKERS_PER_GROUP genes with log2((group + 1) / (highest other + 1)) > 1
   and at least MIN_MARKER_CPM. Candidate genes (the list passed to run and the discovery pool) are never markers, so
   their in-situ attribution is an independent read-out of the deconvolution.
2. Deconvolution: non-negative least squares of each spot's CPM over the markers on the group profiles (each marker
   scaled by its mean); weights normalised to 1 = estimated share of the spot's transcripts per group.
3. In-situ expresser: NNLS of the gene's spot CPM on the spot shares over all UC inflamed spots; beta_k is expression
   per unit share of group k and contribution_k = beta_k x mean share_k (normalised).
4. Composition-corrected source: in-situ mean group share x single-cell profile, normalised (removes dissociation bias:
   the single-cell atlas over-represents immune cells and under-represents epithelium and stroma).
5. Neighbourhood: top decile spots of the gene (>= 1 UMI) per UC section; mean share of each group in their ring-1
   neighbours (6 adjacent spots, about 100 um) against random spots of the same sections (observed / expected),
   PERMUTATIONS section-stratified permutations (per-gene seed), Benjamini-Hochberg across genes x groups.
Outputs: data/processed/spatial_visium_deconvolution.parquet, spatial_visium_expressers.parquet,
spatial_visium_neighbors.parquet, spatial_visium_summary.json, spatial_viewer.json (tissue viewer of the report).
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.optimize import nnls

from .paths import PROCESSED
from .spatial import read_10x_h5
from .spatial_download import VISIUM
from .spatial_visium import D, MIN_SPOTS_EXPRESSING, load

CELL_GROUPS = {
    "Epithelial": ["Enterocytes", "Enterocytes BEST4", "Goblet", "Tuft", "Stem LGR5 SMOC2", "Precursor Cells",
                   "Cycling PCNA", "Cycling TOP2A", "Enterochromaffin Cells", "L Cells"],
    "Fibroblast": ["Activated crypt-top fibroblast", "Adventitial fibroblast", "Crypt-bottom fibroblast", "LP fibroblast",
                   "Resting crypt-top fibroblast"],
    "Inflammatory_fibroblast": ["Inflammatory fibroblast"],
    "Smooth_muscle": ["SMC-1", "SMC-2", "Pericyte"],
    "Endothelial": ["ACKR1 endothelial", "CD36 endothelial", "CXCL12 endothelial", "Lymphatics"],
    "Myeloid": ["Macrophage", "Monocyte-derived macrophage", "Resident macrophage", "MDC", "DC", "CDC1", "CDC2", "PDC",
                "TlsDC"],
    "Neutrophil": ["Neutrophil", "Neutrophil CXCR4", "Neutrophil MX1", "Neutrophil OSM", "Neutrophil PADI4"],
    "T_NK": ["CD4 FOS T Cell", "CD4 IL17 T Cell", "CD4 Naive T Cell", "CD4 PD1 T Cell", "CD4 T Cell",
             "CD8 Cytotoxic T Cell", "CD8 Gamma Delta T Cell", "CD8 IEL", "CD8 KLRG1 T Cell", "Central Memory T Cell",
             "Cycling T Cell", "T Regs", "ILC", "NK"],
    "B": ["Follicular B", "Memory B", "GC-DZ", "GC-LZ", "Cycling B"],
    "Plasma": ["Plasma IgA", "Plasma IgG", "Plasma IgM"],
    "Mast": ["Mast cell"],
    "Eosinophil": ["Eosinophil"],
}
MARKERS_PER_GROUP, MIN_MARKER_CPM = 40, 20.0
HIGH_QUANTILE, MIN_HIGH_SPOTS, PERMUTATIONS, SEED = 0.9, 10, 1000, 20261110
# lineage genes whose correlation with the matching share is reported as a sanity check (they may be markers)
VALIDATION = {"EPCAM": "Epithelial", "COL1A1": "Fibroblast", "PECAM1": "Endothelial", "CD68": "Myeloid",
              "CD3E": "T_NK", "JCHAIN": "Plasma", "FCGR3B": "Neutrophil", "MS4A1": "B", "TPSAB1": "Mast",
              "DES": "Smooth_muscle"}
# candidates whose cell is known from single-cell data; never markers, so this is the independent check
HELD_OUT = {"S100A8": "Neutrophil", "CXCR2": "Neutrophil", "CXCR1": "Neutrophil", "LCN2": "Epithelial",
            "DUOX2": "Epithelial", "MMP3": "Inflammatory_fibroblast", "IL11": "Inflammatory_fibroblast",
            "CHI3L1": "Inflammatory_fibroblast"}
EXTRA = ["S100A8", "LCN2", "DUOX2", "MMP3", "CHI3L1", "MMP1", "IDO1", "CXCL8", "IL11", "OSMR", "CXCR1", "CXCR2",
         "OSM", "TREM1", "IL1B", "TNF", "IL23A", "ITGB7", "CD274", "TNFRSF6B", "IL12B", "JAK2"]
VIEWER_GENES = ["S100A8", "CXCL8", "LCN2", "DUOX2", "MMP3", "MMP1", "CHI3L1", "IDO1", "IL11", "OSMR", "AQP8",
                "CXCR2", "IL1B", "REG1A", "SAA1", "TNFRSF6B", "CD274"]


def reference() -> pd.DataFrame:
    """Genes x cell groups CPM from the SCP3755 atlas (all biopsies)."""
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet",
                         columns=["gene", "health", "cluster", "cp10k_sum", "n_cells"])
    cl2g = {c: g for g, cs in CELL_GROUPS.items() for c in cs}
    sh["group"] = sh["cluster"].map(cl2g)
    sh = sh.dropna(subset=["group"])
    tot = sh.groupby(["gene", "group"])["cp10k_sum"].sum().unstack(fill_value=0)
    cells = sh[["health", "cluster", "group", "n_cells"]].drop_duplicates().groupby("group")["n_cells"].sum()
    prof = (tot / cells.reindex(tot.columns).to_numpy())[list(CELL_GROUPS)]
    return prof / prof.sum(axis=0) * 1e6


def select_markers(ref: pd.DataFrame, exclude: set[str], present: set[str]) -> dict[str, list[str]]:
    x = ref[ref.index.isin(present) & ~ref.index.isin(exclude) & ~ref.index.str.match(r"^(MT-|RPL|RPS|IG[HKL])")]
    out = {}
    for g in x.columns:
        score = np.log2((x[g] + 1) / (x.drop(columns=g).max(axis=1) + 1))
        keep = score[(x[g] >= MIN_MARKER_CPM) & (score > 1)].sort_values(ascending=False)
        out[g] = keep.head(MARKERS_PER_GROUP).index.tolist()
    return out


def all_spots() -> tuple[pd.DataFrame, sparse.csr_matrix, list[str]]:
    """In-tissue spots of the nine sections on a common gene axis, with array and pixel positions."""
    loaded = {gsm: load(gsm) for gsm in VISIUM}
    common = sorted(set.intersection(*[set(g) for _, g in loaded.values()]))
    obs, mats = [], []
    for gsm, (stem, state) in VISIUM.items():
        m, genes = loaded[gsm]
        _, _, bcs, _ = read_10x_h5(D / f"{gsm}_filtered_feature_bc_matrix.h5")
        pos = pd.read_csv(D / f"{gsm}_tissue_positions_list.csv", header=None, index_col=0).reindex(bcs)
        pos = pos[pos[1].fillna(0).astype(int) == 1]
        assert len(pos) == m.shape[0]
        mats.append(m[:, pd.Index(genes).get_indexer(common)])
        obs.append(pd.DataFrame({"spot": [f"{gsm}_{b}" for b in pos.index], "section": gsm, "sample": stem,
                                 "state": state, "x": pos[5].to_numpy(float), "y": pos[4].to_numpy(float)}))
    return pd.concat(obs, ignore_index=True), sparse.vstack(mats).tocsr(), common


def deconvolve(counts: sparse.csr_matrix, genes: list[str], ref: pd.DataFrame, markers: dict) -> np.ndarray:
    mk = list(dict.fromkeys(g for m in markers.values() for g in m))
    col = pd.Index(genes).get_indexer(mk)
    lib = np.asarray(counts.sum(axis=1)).ravel()
    y = (sparse.diags(1e6 / np.maximum(lib, 1)) @ counts[:, col]).toarray()
    R = ref.loc[mk].to_numpy()
    scale = R.mean(axis=1, keepdims=True)
    W = np.vstack([nnls(R / scale, y[i] / scale.ravel())[0] for i in range(len(y))])
    s = W.sum(axis=1, keepdims=True)
    return np.divide(W, s, out=np.zeros_like(W), where=s > 0)


def neighbors(spots: pd.DataFrame) -> sparse.csr_matrix:
    """Row-normalised ring-1 adjacency within each section (centre distance <= 1.25 x median nearest neighbour)."""
    from scipy.spatial import cKDTree
    rows, cols = [], []
    for _, idx in spots.groupby("section").indices.items():
        xy = spots[["x", "y"]].to_numpy()[idx]
        tree = cKDTree(xy)
        d, _ = tree.query(xy, k=2)
        for i, nb in enumerate(tree.query_ball_point(xy, 1.25 * np.median(d[:, 1]))):
            for j in nb:
                if j != i:
                    rows.append(idx[i])
                    cols.append(idx[j])
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(spots), len(spots)))
    return sparse.diags(1 / np.maximum(np.asarray(A.sum(axis=1)).ravel(), 1)) @ A


def _rng(gene: str) -> np.random.Generator:
    return np.random.default_rng([SEED, int(hashlib.sha256(gene.encode()).hexdigest()[:8], 16)])


def _rho(a, b) -> float:
    return round(float(stats.spearmanr(a, b).statistic), 3)


def run(targets: list[str]) -> dict:
    from statsmodels.stats.multitest import multipletests
    ref = reference()
    spots, counts, genes = all_spots()
    disc = set(pd.read_parquet(PROCESSED / "discovery.parquet", columns=["gene"])["gene"])
    exclude = (set(targets) | disc | set(HELD_OUT)) - set(VALIDATION)
    markers = select_markers(ref, exclude, set(genes))
    groups = list(CELL_GROUPS)
    P = deconvolve(counts, genes, ref, markers)
    dec = spots.copy()
    for k, g in enumerate(groups):
        dec[g] = P[:, k].astype("float32")
    dec.to_parquet(PROCESSED / "spatial_visium_deconvolution.parquet", index=False)
    A = neighbors(spots)
    NP = A @ P
    has_nb = np.asarray((A > 0).sum(axis=1)).ravel() > 0
    sections = spots["section"].to_numpy()
    uc = (spots["state"] == "UC_inflamed").to_numpy()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    gi = {g: i for i, g in enumerate(genes)}
    mean_share = P[uc].mean(axis=0)
    ex_rows, nb_rows, summ = [], [], []
    for gene in dict.fromkeys(targets):
        if gene not in gi:
            continue
        umi = counts[:, gi[gene]].toarray().ravel()
        cpm = umi / np.maximum(lib, 1) * 1e6
        if (umi[uc] > 0).sum() < MIN_SPOTS_EXPRESSING:
            continue
        beta = nnls(P[uc], cpm[uc])[0]
        contrib = beta * mean_share
        contrib = contrib / contrib.sum() if contrib.sum() > 0 else contrib
        r = ref.loc[gene].to_numpy() if gene in ref.index else np.zeros(len(groups))
        src = mean_share * r
        src = src / src.sum() if src.sum() > 0 else src
        for k, g in enumerate(groups):
            ex_rows.append({"gene": gene, "cell_group": g, "in_situ_beta_cpm": float(beta[k]),
                            "in_situ_contribution": float(contrib[k]), "composition_corrected_share": float(src[k]),
                            "spearman_with_share": _rho(cpm[uc], P[uc, k])})
        high = np.zeros(len(spots), bool)
        used = []
        for s in np.unique(sections[uc]):
            idx = np.flatnonzero(sections == s)
            if (umi[idx] > 0).sum() < MIN_HIGH_SPOTS:
                continue
            h = idx[(cpm[idx] >= np.quantile(cpm[idx], HIGH_QUANTILE)) & (umi[idx] > 0)]
            if len(h) >= MIN_HIGH_SPOTS:
                high[h] = True
                used.append(s)
        entry = {"gene": gene, "spots_with_umi": int((umi[uc] > 0).sum()), "sections_used": len(used),
                 "n_high_spots": int(high.sum()),
                 "top_in_situ_expresser": groups[int(np.argmax(contrib))] if contrib.sum() > 0 else None,
                 "in_situ_contribution": round(float(contrib.max()), 3) if contrib.sum() > 0 else None,
                 "top_corrected_source": groups[int(np.argmax(src))] if src.sum() > 0 else None}
        if used:
            rng = _rng(gene)
            obs_s, exp_s = np.zeros(len(groups)), np.zeros(len(groups))
            null, n_high = np.zeros((PERMUTATIONS, len(groups))), 0
            for s in used:
                pool = np.flatnonzero((sections == s) & has_nb)
                h = np.flatnonzero(high & (sections == s) & has_nb)
                if len(h) == 0:
                    continue
                n_high += len(h)
                obs_s += NP[h].sum(axis=0)
                exp_s += len(h) * NP[pool].mean(axis=0)
                draws = np.vstack([rng.choice(pool, size=len(h), replace=False) for _ in range(PERMUTATIONS)])
                null += NP[draws].sum(axis=1)
            obs, exp, null = obs_s / n_high, exp_s / n_high, null / n_high
            lr = np.log2((obs + 1e-4) / (exp + 1e-4))
            lrn = np.log2((null + 1e-4) / (exp + 1e-4))
            p = (1 + (np.abs(lrn) >= np.abs(lr)).sum(axis=0)) / (1 + PERMUTATIONS)
            for k, g in enumerate(groups):
                nb_rows.append({"gene": gene, "cell_group": g, "observed_neighbor_share": float(obs[k]),
                                "expected_neighbor_share": float(exp[k]), "log2_enrichment": float(lr[k]),
                                "pvalue": float(p[k]), "n_high_spots": int(n_high), "n_sections": len(used)})
        summ.append(entry)
    nb = pd.DataFrame(nb_rows)
    nb["padj"] = multipletests(nb["pvalue"], method="fdr_bh")[1]
    nb.to_parquet(PROCESSED / "spatial_visium_neighbors.parquet", index=False)
    pd.DataFrame(ex_rows).to_parquet(PROCESSED / "spatial_visium_expressers.parquet", index=False)
    for t in summ:
        s = nb[(nb.gene == t["gene"]) & (nb.padj < 0.05) & (nb.log2_enrichment > 0)].sort_values("log2_enrichment",
                                                                                                  ascending=False)
        t["enriched_neighbors"] = [{"cell_group": r.cell_group, "log2_enrichment": round(r.log2_enrichment, 3),
                                    "padj": float(f"{r.padj:.3g}")} for r in s.itertuples()]

    def check(table):
        out = {}
        for g, grp in table.items():
            if g in gi:
                cpm = counts[:, gi[g]].toarray().ravel() / np.maximum(lib, 1) * 1e6
                out[g] = {"group": grp, "spearman": _rho(cpm[uc], P[uc, groups.index(grp)])}
        return out
    mk_all = {m for ms in markers.values() for m in ms}
    val = check(VALIDATION)
    for g in val:
        val[g]["marker"] = g in mk_all
    summary = {"spots": int(len(spots)), "uc_spots": int(uc.sum()), "sections": int(spots.section.nunique()),
               "groups": groups, "markers": {g: len(m) for g, m in markers.items()}, "marker_genes": markers,
               "mean_share_uc": dict(zip(groups, np.round(mean_share, 3).tolist())),
               "mean_share_healthy": dict(zip(groups, np.round(P[~uc].mean(axis=0), 3).tolist())),
               "validation_markers": val, "validation_held_out": check(HELD_OUT), "targets": summ}
    (PROCESSED / "spatial_visium_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return {k: v for k, v in summary.items() if k not in ("marker_genes", "targets")}


def viewer(genes: list[str], max_px: int = 900, expr_scale: float = 20.0) -> dict:
    """Tissue-viewer data: cropped H&E (JPEG, base64), spot positions, shares, gene log1p CPM and neighbours."""
    import base64
    import io
    from PIL import Image
    dec = pd.read_parquet(PROCESSED / "spatial_visium_deconvolution.parquet")
    groups = list(CELL_GROUPS)
    nb = pd.read_parquet(PROCESSED / "spatial_visium_neighbors.parquet")
    genes = [g for g in dict.fromkeys(genes) if g in set(nb.gene)]
    secs = []
    for gsm, (stem, state) in VISIUM.items():
        m, gnames = load(gsm)
        gi = {g: i for i, g in enumerate(gnames)}
        sf = json.loads((D / f"{gsm}_scalefactors_json.json").read_text())
        img = Image.open(D / f"{gsm}_tissue_lowres_image.png").convert("RGB")
        s = dec[dec.section == gsm].reset_index(drop=True)
        f0 = float(sf["tissue_lowres_scalef"])
        xy = s[["x", "y"]].to_numpy() * f0
        diam = float(sf["spot_diameter_fullres"]) * f0
        pad = diam * 2
        x0, y0 = max(int(xy[:, 0].min() - pad), 0), max(int(xy[:, 1].min() - pad), 0)
        x1, y1 = min(int(xy[:, 0].max() + pad), img.size[0]), min(int(xy[:, 1].max() + pad), img.size[1])
        crop = img.crop((x0, y0, x1, y1))
        f = max_px / max(crop.size)
        crop = crop.resize((max(int(crop.size[0] * f), 1), max(int(crop.size[1] * f), 1)), Image.LANCZOS)
        buf = io.BytesIO()
        crop.save(buf, format="JPEG", quality=72, optimize=True)
        lib = np.asarray(m.sum(axis=1)).ravel()
        expr = {g: np.round(np.log1p(m[:, gi[g]].toarray().ravel() / np.maximum(lib, 1) * 1e6) * expr_scale)
                .astype(int).tolist() for g in genes if g in gi}
        A = neighbors(s).tocsr()
        secs.append({"id": gsm, "region": "UC inflamed" if state == "UC_inflamed" else "healthy", "donor": stem,
                     "n_spots": len(s), "width": crop.size[0], "height": crop.size[1], "spot_px": round(diam * f, 2),
                     "image": "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode(),
                     "x": np.round((xy[:, 0] - x0) * f, 1).tolist(), "y": np.round((xy[:, 1] - y0) * f, 1).tolist(),
                     "shares": {g: np.round(s[g].to_numpy() * 100).astype(int).tolist() for g in groups},
                     "expr": expr, "neighbors": [A.indices[A.indptr[i]:A.indptr[i + 1]].tolist() for i in range(len(s))]})
    st = {g: [{"cell_group": r.cell_group, "log2": round(float(r.log2_enrichment), 3), "padj": float(f"{r.padj:.3g}"),
               "observed": round(float(r.observed_neighbor_share), 4), "expected": round(float(r.expected_neighbor_share), 4)}
              for r in d.itertuples()] for g, d in nb[nb.gene.isin(genes)].groupby("gene")}
    n_high = nb.groupby("gene")[["n_high_spots", "n_sections"]].first().to_dict("index")
    summ = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())
    out = {"groups": groups, "genes": genes, "expr_scale": expr_scale, "sections": secs, "stats": st,
           "stats_n": {g: n_high[g] for g in genes}, "validation": summ["validation_held_out"]}
    p = PROCESSED / "spatial_viewer.json"
    p.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    return {"sections": len(secs), "genes": len(genes), "megabytes": round(p.stat().st_size / 1e6, 1)}


def load_graph(store=None) -> int:
    """One spatial_visium_cell record per gene: deconvolved in-situ expresser and enriched neighbours."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    summ = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())
    ex = pd.read_parquet(PROCESSED / "spatial_visium_expressers.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE189184", "Dataset", "10x Visium whole transcriptome, UC inflamed and healthy colon (GSE189184)")
        for t in summ["targets"]:
            node = store.resolve(t["gene"])
            if not node or node["type"] != "Gene":
                continue
            e = ex[ex.gene == t["gene"]].sort_values("in_situ_contribution", ascending=False)
            rec = {"gene": t["gene"], "platform": "10x Visium, spots deconvolved with SCP3755 single-cell profiles (NNLS)",
                   "top_in_situ_expresser": t["top_in_situ_expresser"],
                   "in_situ_contribution": t["in_situ_contribution"],
                   "in_situ_expressers": [{"cell_group": r.cell_group, "contribution": round(r.in_situ_contribution, 3),
                                           "spearman_with_share": r.spearman_with_share}
                                          for r in e.head(3).itertuples() if r.in_situ_contribution > 0.05],
                   "composition_corrected_source": t["top_corrected_source"],
                   "enriched_neighbors": t["enriched_neighbors"][:4], "n_high_spots": t["n_high_spots"],
                   "sections_used": t["sections_used"], "uc_spots": summ["uc_spots"]}
            ev = store.add_evidence(evidence_id("spatial_visium_cell", t["gene"]), "spatial_visium_cell",
                                    "DATASET:GSE189184", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "observed_in_situ_in", "DATASET:GSE189184", ev)
            n += 1
    return n


if __name__ == "__main__":
    from .cli import expression_top
    extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
    print(json.dumps(run(list(dict.fromkeys(expression_top(300) + extra + EXTRA))), indent=1))
    print(viewer(VIEWER_GENES))
    s = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())
    for t in s["targets"][:60]:
        print(t["gene"], t["top_in_situ_expresser"], t["in_situ_contribution"], t["top_corrected_source"],
              [(e["cell_group"], e["log2_enrichment"]) for e in t["enriched_neighbors"]][:3])
