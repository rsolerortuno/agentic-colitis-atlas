"""M11: fourth single-cell atlas, from the full count matrices of GSE282122 (Thomas et al. 2024, Nat Immunol 25:2152;
added 2026-11-11).

The anti-TNF atlas (adalimumab, Crohn's disease, UC and healthy gut, 216 biopsies) is also used through the authors'
cell-state tables (antitnf_atlas.py). Here its raw data are used as a fourth atlas for the cell attribution. GEO holds
one filtered 10x matrix (.h5) per biopsy in a 2.8 GB tar and the sample metadata in the series matrix (disease,
patient, site, inflammation, before or after treatment); it has no cell annotations, so the cells are typed here.
Remission after adalimumab comes from the authors' Supplementary Table 1A iii.

Samples: the 108 UC biopsies and the 12 healthy ones (colon sites only; ileum excluded). Cells with < 200 detected
genes or > 50 % mitochondrial reads are dropped (gut epithelium runs high on mitochondrial reads).

Cell typing (no authors' annotation): every cell is scored on a panel of canonical markers of nine lineages
(epithelial, stromal, endothelial, myeloid, T/NK, B, plasma, mast, neutrophil; MARKERS), log1p of counts per 10,000.
The panel is standardised, reduced to 15 principal components (SVD of a random 100,000-cell subsample) and clustered
with k-means (k = 60, k-means++ seeding, fixed seed); each cluster takes the lineage whose markers have the highest
mean z-score in it. Clustering on a marker panel is coarse: it types lineages, not cell states, which is the level
the replication works at. A cluster whose best and second lineage scores differ by < 0.25 z is called "mixed"
(likely doublets) and left out. Sanity check: for every canonical marker, the lineage with its highest mean
expression is its own (typing_check(), written to the summary).
Result (2026-11-11): 694,801 cells typed in 60 clusters; every canonical marker of the eight lineages found is highest in
its own lineage. No neutrophil cluster emerges (neutrophils carry little RNA and are mostly removed by droplet
filtering); 51,036 cells in four clusters without a clear lineage are left out. The series mixes two Cell Ranger
references (33,538 and 33,808 features), so genes are matched across samples by Ensembl id.

Per gene and lineage: transcripts summed after scaling each cell to 10,000 counts, in the format of
sc_celltype_share (gene, health, cluster, compartment, cp10k_sum, n_expressing, n_cells, mean_cp10k, pct_expressing,
share; cluster = compartment = lineage; health inflamed / noninflamed / healthy), and the same by treatment time
(pre / post) and remission for UC.
Outputs (data/processed/): sc4_pseudobulk.parquet (patient x health x time x lineage), sc4_celltype_share.parquet,
sc4_response_share.parquet, sc4_clusters.parquet, sc4_summary.json.
"""
from __future__ import annotations

import gzip
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

D = RAW / "gse282122"
TAR = D / "GSE282122_filtered_processed_data.tar.gz"
H5 = D / "h5"   # UC and healthy colon samples extracted from the tar (derived, not in the manifest)
STUDY = "Thomas et al. 2024, Nat Immunol (GSE282122), full count matrices"
HEALTH = {"Inflamed": "inflamed", "Non_Inflamed": "noninflamed", "Healthy": "healthy"}
LINEAGES = ["epithelial", "stromal", "endothelial", "myeloid", "T/NK", "B", "plasma", "mast", "neutrophil"]
MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT20", "CDH1", "PIGR", "CLDN7"],
    "stromal": ["COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A2", "SFRP2"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5", "EGFL7"],
    "myeloid": ["CD68", "CD14", "C1QA", "C1QB", "AIF1", "CSF1R", "LYZ"],
    "T/NK": ["CD3E", "CD3D", "CD2", "TRAC", "NKG7", "GNLY", "CD7"],
    "B": ["MS4A1", "CD79B", "CD19", "BANK1", "CD22", "PAX5"],
    "plasma": ["MZB1", "JCHAIN", "IGHA1", "DERL3", "TNFRSF17", "IGKC"],
    "mast": ["TPSAB1", "CPA3", "KIT", "MS4A2", "HDC", "TPSB2"],
    "neutrophil": ["FCGR3B", "CXCR2", "CSF3R", "CMTM2", "FFAR2", "PROK2"],
}
K, NPC, SEED, MIXED_Z = 60, 15, 282122, 0.25


def samples() -> pd.DataFrame:
    """GEO sample table with remission (UC patients) from the authors' Supplementary Table 1A iii."""
    rows: dict[str, list[str]] = {}
    with gzip.open(D / "GSE282122_series_matrix.txt.gz", "rt", encoding="utf-8") as f:
        for line in f:
            p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if p[0] in ("!Sample_title", "!Sample_geo_accession"):
                rows[p[0][8:]] = p[1:]
            elif p[0] == "!Sample_characteristics_ch1":
                k = p[1].split(": ", 1)[0].lower()
                rows[k] = [x.split(": ", 1)[1] if ": " in x else "" for x in p[1:]]
    s = pd.DataFrame(rows).rename(columns={"geo_accession": "gsm"})
    from .antitnf_atlas import SRC, _sheet
    if SRC.exists():
        u = _sheet("41590_2024_1994_MOESM3_ESM.xlsx", "Supp_Table_1A_iii_UCspecific")
        u = u[u.Patient.astype(str).str.match(r"UC\d+$")].dropna(subset=["Remission"]).drop_duplicates("Patient")
        rem = dict(zip(u.Patient, u.Remission.str.lower()))
        s["response"] = s.patient.map(rem).fillna("")
    else:
        s["response"] = ""
    s["health"] = s.inflammation.map(HEALTH)
    s["time"] = s.treatment.str.lower().replace({"none": ""})
    return s


def kept_samples() -> pd.DataFrame:
    s = samples()
    return s[s.disease.isin(["UC", "Healthy"]) & (s.site != "Terminal_Ileum")].reset_index(drop=True)


def extract() -> dict:
    """Stream the tar once and write the h5 files of the kept samples to data/raw/gse282122/h5/."""
    keep = kept_samples()
    titles = set(keep.title)
    H5.mkdir(parents=True, exist_ok=True)
    done, seen = 0, []
    with tarfile.open(TAR, mode="r|gz") as t:
        for m in t:
            if not m.isfile() or not m.name.endswith(".h5"):
                continue
            seen.append(m.name)
            title = _title_of(m.name, titles)
            if title is None:
                continue
            dest = H5 / f"{title}.h5"
            if dest.exists() and dest.stat().st_size == m.size:
                done += 1
                continue
            f = t.extractfile(m)
            tmp = dest.with_suffix(".part")
            with open(tmp, "wb") as out:
                while chunk := f.read(1 << 22):
                    out.write(chunk)
            tmp.replace(dest)
            done += 1
    (D / "h5_members.txt").write_text("\n".join(seen), encoding="utf-8")
    return {"h5_in_tar": len(seen), "extracted": done, "kept_samples": len(keep)}


def _title_of(member: str, titles: set[str]) -> str | None:
    """Sample title of a tar member (filtered_processed_data/<sample>/filtered_feature_bc_matrix.h5); some GEO titles
    carry a "-reup" suffix (re-uploaded) that the folder name may lack."""
    alias = {t: t for t in titles} | {re.sub(r"-reup$", "", t): t for t in titles}
    for part in Path(member).parts[:-1]:
        if part in alias:
            return alias[part]
    return None


def read_h5(path: Path):
    """10x HDF5 (Cell Ranger v3+): CSC matrix genes x cells -> (data, indices, indptr, n_genes, symbols, ids)."""
    import h5py
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = m["shape"][:]
        names = m["features"]["name"][:].astype(str)
        ids = m["features"]["id"][:].astype(str)
        ftype = m["features"]["feature_type"][:].astype(str)
        return (m["data"][:].astype(np.float64), m["indices"][:].astype(np.int64), m["indptr"][:].astype(np.int64),
                int(shape[0]), names, ids, ftype)


def _qc(data, indices, indptr, names):
    col = np.repeat(np.arange(len(indptr) - 1), np.diff(indptr))
    total = np.bincount(col, weights=data, minlength=len(indptr) - 1)
    mito = np.char.startswith(names.astype(str), "MT-")
    mt = np.bincount(col, weights=data * mito[indices], minlength=len(indptr) - 1)
    ngenes = np.diff(indptr)
    ok = (ngenes >= 200) & (mt <= 0.5 * np.maximum(total, 1))
    return col, total, ok


def _panel(data, indices, col, total, ok, names, panel):
    """log1p(cp10k) of the panel genes for the QC-passing cells (cells x panel, float32)."""
    pos = {g: i for i, g in enumerate(names)}
    gi = np.array([pos.get(g, -1) for g in panel])
    lut = np.full(len(names), -1)
    lut[gi[gi >= 0]] = np.where(gi >= 0)[0]
    sel = lut[indices] >= 0
    X = np.zeros((len(total), len(panel)), np.float32)
    np.add.at(X, (col[sel], lut[indices[sel]]), (data[sel] * 1e4 / total[col[sel]]).astype(np.float32))
    return np.log1p(X[ok])


def _kmeans(Z, k, seed):
    from scipy.cluster.vq import kmeans2
    return kmeans2(Z, k, minit="++", seed=seed, iter=30)


def type_cells(files: list[Path]) -> tuple[dict, pd.DataFrame]:
    """Pass 1: marker panel per cell, PCA + k-means, cluster -> lineage. Returns {sample: lineage per cell (or '')}."""
    panel = [g for v in MARKERS.values() for g in v]
    Xs, ns, oks = [], [], {}
    for p in files:
        data, indices, indptr, ng, names, ids, ftype = read_h5(p)
        col, total, ok = _qc(data, indices, indptr, names)
        Xs.append(_panel(data, indices, col, total, ok, names, panel))
        ns.append(int(ok.sum()))
        oks[p.stem] = ok
    X = np.vstack(Xs)
    mu, sd = X.mean(0), X.std(0) + 1e-6
    Z = np.clip((X - mu) / sd, -10, 10)
    rng = np.random.default_rng(SEED)
    sub = rng.choice(len(Z), size=min(100_000, len(Z)), replace=False)
    _, _, vt = np.linalg.svd(Z[sub] - Z[sub].mean(0), full_matrices=False)
    P = (Z - Z[sub].mean(0)) @ vt[:NPC].T
    cent, _ = _kmeans(P[sub].astype(np.float64), K, SEED)
    lab = np.empty(len(P), np.int32)
    for i in range(0, len(P), 50_000):   # assign every cell to its nearest centroid
        dd = ((P[i:i + 50_000, None, :] - cent[None]) ** 2).sum(-1)
        lab[i:i + 50_000] = dd.argmin(1)
    # cluster -> lineage by the mean z-score of each lineage's markers
    gidx = {L: [panel.index(g) for g in MARKERS[L]] for L in LINEAGES}
    rows = []
    for c in range(K):
        m = lab == c
        if not m.any():
            continue
        zc = Z[m].mean(0)
        sc = {L: float(zc[gidx[L]].mean()) for L in LINEAGES}
        best = sorted(sc, key=sc.get, reverse=True)
        lin = best[0] if sc[best[0]] - sc[best[1]] >= MIXED_Z else "mixed"
        rows.append({"cluster": c, "n_cells": int(m.sum()), "lineage": lin, "best": best[0], "second": best[1],
                     "margin": round(sc[best[0]] - sc[best[1]], 3), **{f"z_{L}": round(v, 3) for L, v in sc.items()}})
    cl = pd.DataFrame(rows)
    lin_of = dict(zip(cl.cluster, cl.lineage))
    cell_lin = np.array([lin_of.get(c, "mixed") for c in lab], dtype=object)
    out, i = {}, 0
    for p, n in zip(files, ns):
        full = np.full(len(oks[p.stem]), "", dtype=object)
        full[oks[p.stem]] = cell_lin[i:i + n]
        out[p.stem] = full
        i += n
    # sanity check data: mean log1p(cp10k) of each marker by lineage
    chk = pd.DataFrame(X, columns=panel).assign(lineage=cell_lin).groupby("lineage").mean()
    chk.to_parquet(PROCESSED / "sc4_marker_check.parquet")
    return out, cl


def typing_check() -> dict:
    """For every canonical marker: is the lineage with its highest mean expression its own lineage?"""
    chk = pd.read_parquet(PROCESSED / "sc4_marker_check.parquet").drop(index="mixed", errors="ignore")
    res = {}
    for L, genes in MARKERS.items():
        if L not in chk.index:
            res[L] = {"present": False}
            continue
        top = {g: chk[g].idxmax() for g in genes}
        res[L] = {"present": True, "own_top": sum(v == L for v in top.values()), "markers": len(genes),
                  "misses": {g: v for g, v in top.items() if v != L}}
    return res


def build() -> dict:
    from .expression import resolve_symbol
    keep = kept_samples()
    files = [H5 / f"{t}.h5" for t in keep.title if (H5 / f"{t}.h5").exists()]
    if len(files) < len(keep):
        print(f"warning: {len(keep) - len(files)} kept samples have no h5", flush=True)
    lin_cells, cl = type_cells(files)
    cl.to_parquet(PROCESSED / "sc4_clusters.parquet", index=False)
    meta = keep.set_index("title")
    LI = {L: i for i, L in enumerate(LINEAGES)}
    nL = len(LINEAGES)
    acc: dict[tuple, dict] = {}
    # two Cell Ranger references are mixed in the series (33,538 and 33,808 features): genes are matched by Ensembl id
    feats: dict[str, tuple[str, bool]] = {}
    for p in files:
        import h5py
        with h5py.File(p, "r") as f:
            fe = f["matrix"]["features"]
            for i, n, t in zip(fe["id"][:].astype(str), fe["name"][:].astype(str), fe["feature_type"][:].astype(str)):
                feats.setdefault(i, (n, t == "Gene Expression"))
    gid = {i: k for k, i in enumerate(feats)}
    ng = len(gid)
    sym = np.array([resolve_symbol(n) or n for n, _ in feats.values()], dtype=object)
    gene_ok = np.array([ok for _, ok in feats.values()])
    for p in files:
        data, indices, indptr, _, names, ids, ftype = read_h5(p)
        indices = np.array([gid[i] for i in ids])[indices]
        lin = lin_cells[p.stem]
        code = np.array([LI.get(x, -1) for x in lin])
        col = np.repeat(np.arange(len(indptr) - 1), np.diff(indptr))
        total = np.bincount(col, weights=data, minlength=len(indptr) - 1)
        c = code[col]
        sel = c >= 0
        scaled = data[sel] * 1e4 / total[col[sel]]
        idx = indices[sel] * nL + c[sel]
        r = meta.loc[p.stem]
        key = (r.patient, r.health, r.time, r.response)
        a = acc.setdefault(key, {"sum": np.zeros(ng * nL), "nexp": np.zeros(ng * nL), "ncell": np.zeros(nL)})
        a["sum"] += np.bincount(idx, weights=scaled, minlength=ng * nL)
        a["nexp"] += np.bincount(idx, minlength=ng * nL)
        a["ncell"] += np.bincount(code[code >= 0], minlength=nL)
    parts = []
    for (pat, health, time, resp), a in acc.items():
        s = a["sum"].reshape(-1, nL)
        gi, li = np.nonzero(s)
        keepg = gene_ok[gi]
        gi, li = gi[keepg], li[keepg]
        parts.append(pd.DataFrame({"gene": sym[gi], "subject": pat, "health": health, "time": time, "response": resp,
                                   "cluster": np.array(LINEAGES, dtype=object)[li],
                                   "cp10k_sum": s[gi, li].astype(np.float32),
                                   "n_expressing": a["nexp"].reshape(-1, nL)[gi, li].astype(np.int32),
                                   "n_cells": a["ncell"][li].astype(np.int32)}))
    pb = pd.concat(parts, ignore_index=True)
    pb["compartment"] = pb["cluster"]
    # a few symbols map to the same approved name: sum them
    pb = pb.groupby(["gene", "subject", "health", "time", "response", "cluster", "compartment", "n_cells"],
                    as_index=False, observed=True)[["cp10k_sum", "n_expressing"]].sum()
    pb.to_parquet(PROCESSED / "sc4_pseudobulk.parquet", index=False)
    share(pb, ["health"]).to_parquet(PROCESSED / "sc4_celltype_share.parquet", index=False)
    share(pb[pb.health != "healthy"], ["health", "time", "response"]).to_parquet(
        PROCESSED / "sc4_response_share.parquet", index=False)
    cells = pb.drop_duplicates(["subject", "health", "time", "cluster"])
    s = {"study": STUDY, "samples": len(files), "uc_samples": int((keep.disease == "UC").sum()),
         "healthy_samples": int((keep.disease == "Healthy").sum()),
         "patients": int(keep.patient.nunique()), "cells_typed": int(cells.n_cells.sum()),
         "cells_by_lineage": cells.groupby("cluster").n_cells.sum().astype(int).to_dict(),
         "cells_mixed": int(sum((v == "mixed").sum() for v in lin_cells.values())),
         "cells_failed_qc": int(sum((v == "").sum() for v in lin_cells.values())),
         "clusters": int(len(cl)), "clusters_by_lineage": cl.lineage.value_counts().to_dict(),
         "typing_check": typing_check(),
         "uc_patients_by_response": keep[keep.disease == "UC"].drop_duplicates("patient").response.replace("", "unknown")
         .value_counts().to_dict()}
    (PROCESSED / "sc4_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def share(pb: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """The sc_celltype_share format, grouped by `by` (+ gene, cluster, compartment)."""
    g = pb.groupby(["gene", *by, "cluster", "compartment"], as_index=False, observed=True)[["cp10k_sum", "n_expressing"]].sum()
    cells = pb.drop_duplicates(["subject", "health", "time", "cluster"]).groupby([*by, "cluster"], observed=True).n_cells.sum()
    g["n_cells"] = [cells.get(tuple(k) if len(k) > 1 else k[0], np.nan)
                    for k in g[[*by, "cluster"]].itertuples(index=False, name=None)]
    g["mean_cp10k"] = g.cp10k_sum / g.n_cells
    g["pct_expressing"] = 100 * g.n_expressing / g.n_cells
    g["share"] = g.cp10k_sum / g.groupby(["gene", *by], observed=True).cp10k_sum.transform("sum")
    return g


def gene_view(gene: str, sh: pd.DataFrame | None = None, rs: pd.DataFrame | None = None) -> dict:
    sh = pd.read_parquet(PROCESSED / "sc4_celltype_share.parquet", filters=[("gene", "==", gene)]) if sh is None else sh
    rs = pd.read_parquet(PROCESSED / "sc4_response_share.parquet", filters=[("gene", "==", gene)]) if rs is None else rs
    out = {}
    for h in ("inflamed", "noninflamed", "healthy"):
        d = sh[sh.health == h].sort_values("share", ascending=False)
        out[f"lineage_share_{h}"] = {r.cluster: round(float(r.share), 3) for r in d.itertuples() if r.share >= 0.01}
        out[f"mean_cp10k_{h}"] = {r.cluster: round(float(r.mean_cp10k), 3) for r in d.head(3).itertuples()}
    resp = {}
    for (h, t, r_), d in rs[rs.response != ""].groupby(["health", "time", "response"]):
        top = d.sort_values("share", ascending=False).iloc[0]
        resp[f"{h} {t} {r_}"] = {"dominant_lineage": top.cluster, "share": round(float(top.share), 3),
                                 "mean_cp10k": round(float(top.mean_cp10k), 3)}
    out["by_time_and_remission"] = resp
    return out


def load(store=None) -> dict:
    """One sc_lineage_gse282122 record per gene of the meta-analysis (new evidence type: no existing record changes)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    genes = set(pd.read_parquet(PROCESSED / "meta.parquet")["gene"])
    sh = pd.read_parquet(PROCESSED / "sc4_celltype_share.parquet")
    sh = sh[sh.gene.isin(genes)]
    rs = pd.read_parquet(PROCESSED / "sc4_response_share.parquet")
    rs = rs[rs.gene.isin(genes)]
    rs_g = dict(tuple(rs.groupby("gene")))
    s = json.loads((PROCESSED / "sc4_summary.json").read_text(encoding="utf-8"))
    n = 0
    with store.batch():
        if store.get_node("DATASET:GSE282122") is None:
            store.add_node("DATASET:GSE282122", "Dataset", f"Anti-TNF (adalimumab) longitudinal single-cell atlas of IBD ({STUDY})")
        for gene, d in sh.groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": gene, "study": STUDY, "disease": "ulcerative colitis (and healthy colon)",
                   "cell_typing": "lineages typed in this project from canonical markers (no authors' annotation on GEO)",
                   "cells_typed": s["cells_typed"],
                   **gene_view(gene, d, rs_g.get(gene, rs.iloc[:0]))}
            ev = store.add_evidence(evidence_id("sc_lineage_gse282122", gene), "sc_lineage_gse282122",
                                    "DATASET:GSE282122", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "expressed_by_lineage_in", "DATASET:GSE282122", ev)
            n += 1
    return {"sc_lineage_gse282122_records": n}


def run() -> dict:
    if not H5.exists() or len(list(H5.glob("*.h5"))) < len(kept_samples()):
        print(extract(), flush=True)
    return build()


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
