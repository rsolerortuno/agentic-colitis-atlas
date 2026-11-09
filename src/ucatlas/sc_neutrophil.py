"""M11: third single-cell dataset, with neutrophils (Gubatan et al. 2026, Front Immunol; Broad SCP3755).

408,930 colon cells from 20 UC patients, inflamed and uninflamed biopsies (BD Rhapsody), 47,516 neutrophils in five
states. The portal provides a Seurat object (.RDS) and the cell metadata (.csv). No R is installed: the RDS is
parsed with the pure-Python `rdata` package, the RNA count matrix (dgCMatrix, genes x cells) is located in the
object and summed into pseudobulk per patient x inflammation x cell state, after scaling each cell to 10,000 counts.

Outputs (data/processed/): sc3_pseudobulk.parquet (same columns as sc_pseudobulk) and sc3_celltype_share.parquet.
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

D = RAW / "scp3755"
STEM = "UC.tissue.merged.annotated.all.objects.harmonized.finalversion.updated.metadata.09.18.2023_upload_ready"
HEALTH = {"Inflamed": "inflamed", "Uninflamed": "noninflamed"}


def _walk(obj, path="", depth=0, found=None):
    """Find every dgCMatrix in the parsed R object tree: (path, dims, nnz, the object)."""
    from rdata.parser import RObjectType
    found = [] if found is None else found
    if depth > 40 or obj is None:
        return found
    attrs = {}
    a = getattr(obj, "attributes", None)
    while a is not None and getattr(a, "info", None) is not None and a.info.type == RObjectType.LIST:
        tag = a.tag.value.value.decode() if a.tag is not None and a.tag.value is not None else None
        attrs[tag] = a.value[0]
        a = a.value[1] if len(a.value) > 1 else None
    cls = attrs.get("class")
    if cls is not None and getattr(cls, "value", None) is not None:
        try:
            names = [c.value.decode() if hasattr(c, "value") else str(c) for c in cls.value]
        except Exception:
            names = []
        if "dgCMatrix" in names:
            found.append((path, attrs))
            return found
    for k, v in attrs.items():
        if k not in ("class",):
            _walk(v, f"{path}@{k}", depth + 1, found)
    val = getattr(obj, "value", None)
    if isinstance(val, list):
        for i, v in enumerate(val):
            if hasattr(v, "info"):
                _walk(v, f"{path}[{i}]", depth + 1, found)
    return found


def _arr(robj):
    return np.asarray(robj.value)


def load_counts():
    import rdata
    print("parsing RDS (this takes a while)...", flush=True)
    parsed = rdata.parser.parse_file(D / f"{STEM}.RDS", extension="rds")
    mats = _walk(parsed.object)
    print("dgCMatrix objects:", [(p, _arr(a["Dim"]).tolist(), len(_arr(a["x"]))) for p, a in mats], flush=True)
    meta = pd.read_csv(D / f"{STEM}.csv", index_col=0)
    cand = [(p, a) for p, a in mats if _arr(a["Dim"])[1] == len(meta)]
    # raw counts: integer-valued x; prefer the path that names counts
    def is_int(a):
        x = _arr(a["x"])[:100000]
        return np.allclose(x, np.round(x))
    cand.sort(key=lambda pa: (not is_int(pa[1]), "counts" not in pa[0]))
    path, a = cand[0]
    dimnames = a["Dimnames"].value
    genes = np.array([s.value.decode() for s in dimnames[0].value])
    cells = np.array([s.value.decode() for s in dimnames[1].value])
    print("using", path, "integer:", is_int(a), flush=True)
    return genes, cells, _arr(a["i"]), _arr(a["p"]), _arr(a["x"]), meta


def build() -> dict:
    from .expression import resolve_symbol
    genes, cells, i, p, x, meta = load_counts()
    meta = meta.reindex(cells)
    key = (meta["PatientID"].astype(str) + "|" + meta["Biopsy_INFLAMMATION"].map(HEALTH).astype(str) + "|"
           + meta["Annotations_Level_3"].astype(str))
    codes, groups = pd.factorize(key)
    G = len(groups)
    n_per_cell = np.diff(p)
    col = np.repeat(np.arange(len(cells)), n_per_cell)
    totals = np.bincount(col, weights=x, minlength=len(cells))
    scaled = x * (1e4 / totals[col])
    idx = i.astype(np.int64) * G + codes[col]
    sums = np.bincount(idx, weights=scaled, minlength=len(genes) * G)
    nexp = np.bincount(idx, minlength=len(genes) * G)
    ncell = np.bincount(codes, minlength=G)
    gi, gr = np.nonzero(sums.reshape(len(genes), G))
    flat = gi * G + gr
    parts = pd.Series(groups).str.split("|", expand=True)
    sym = pd.Series(genes).map(lambda g: resolve_symbol(g) or g).to_numpy()
    lvl1 = meta.groupby("Annotations_Level_3")["Annotations_Level_1"].first()
    pb = pd.DataFrame({"gene": sym[gi], "subject": parts[0].to_numpy()[gr], "health": parts[1].to_numpy()[gr],
                       "cluster": parts[2].to_numpy()[gr], "compartment": parts[2].map(lvl1).to_numpy()[gr],
                       "cp10k_sum": sums[flat].astype(np.float32), "n_expressing": nexp[flat].astype(np.int32),
                       "n_cells": ncell[gr].astype(np.int32)})
    pb.to_parquet(PROCESSED / "sc3_pseudobulk.parquet", index=False)
    g = pb.groupby(["gene", "health", "cluster", "compartment"], as_index=False)[["cp10k_sum", "n_expressing"]].sum()
    cells_n = pb.drop_duplicates(["subject", "health", "cluster"]).groupby(["health", "cluster"])["n_cells"].sum()
    g["n_cells"] = [cells_n.get((h, c), np.nan) for h, c in zip(g["health"], g["cluster"])]
    g["mean_cp10k"] = g["cp10k_sum"] / g["n_cells"]
    g["pct_expressing"] = 100 * g["n_expressing"] / g["n_cells"]
    g["share"] = g["cp10k_sum"] / g.groupby(["gene", "health"])["cp10k_sum"].transform("sum")
    g.to_parquet(PROCESSED / "sc3_celltype_share.parquet", index=False)
    return {"cells": len(cells), "genes": len(genes), "groups": G, "pseudobulk_rows": len(pb)}


if __name__ == "__main__":
    print(build())
    sys.exit(0)
