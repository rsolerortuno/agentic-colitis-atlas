"""M11: does the cell-type attribution of Smillie 2019 hold in a second dataset (Garrido-Trigo et al. 2023)?

The second dataset has five coarse lineages, so the comparison is made at lineage level: Smillie clusters are mapped
to epithelial / stromal / T and NK / plasma / myeloid (B cells, absent from the second dataset, are left out). For
each gene the dominant lineage (largest share of the gene's transcripts in UC mucosa) is compared between datasets.
Output: data/processed/sc_replication.parquet.
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
