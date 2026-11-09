"""M11: single-cell layer from Smillie et al. 2019 (SCP259, colon mucosa of UC patients and healthy donors).

The three Matrix Market files (epithelial, stromal, immune; genes x cells, raw UMI counts) are streamed in chunks
and summed into pseudobulk per subject x health state x cell cluster, after scaling each cell to 10,000 UMI
(CP10K, the per-cell nUMI comes from the portal metadata). Nothing is held as a full matrix.

Outputs (data/processed/):
- sc_pseudobulk.parquet: gene, subject, health, cluster, compartment, cp10k_sum, n_expressing, n_cells
- sc_celltype_share.parquet: per gene and health state, the share of the gene's transcripts made by each cell
  cluster (pooled over donors) and the cluster's mean CP10K
- sc_celltype_de.parquet: per gene and cluster, inflamed vs healthy donors (pseudobulk mean CP10K per donor,
  log2(x+1), Welch t-test, BH)
Caveat (as in the paper): epithelial and lamina propria fractions were dissociated separately, so cell proportions
reflect the protocol, not the tissue; shares are comparable between health states, not absolute.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .paths import PROCESSED, RAW

SC = RAW / "scp259"
COMPARTMENTS = {"Epi": "epithelial", "Fib": "stromal", "Imm": "immune"}
HEALTH = {"Healthy": "healthy", "Non-inflamed": "noninflamed", "Inflamed": "inflamed"}


def meta() -> pd.DataFrame:
    m = pd.read_csv(SC / "all.meta2.txt", sep="\t", skiprows=[1])
    m["health"] = m["Health"].map(HEALTH)
    return m.set_index("NAME")


def _pseudobulk(comp: str, m: pd.DataFrame, chunksize: int = 20_000_000) -> pd.DataFrame:
    cells = pd.read_csv(SC / f"{comp}.barcodes2.tsv", header=None)[0]
    genes = pd.read_csv(SC / f"{comp}.genes.tsv", header=None)[0]
    # the 2019 annotation uses older symbols (IL8 for CXCL8): map to approved HGNC symbols where unambiguous
    from .expression import resolve_symbol
    genes = genes.map(lambda g: resolve_symbol(g) or g)
    cm = m.reindex(cells)
    key = cm["Subject"].astype(str) + "|" + cm["health"].astype(str) + "|" + cm["Cluster"].astype(str)
    codes, groups = pd.factorize(key)
    n_groups = len(groups)
    scale = (1e4 / cm["nUMI"].to_numpy(float))
    sums = np.zeros(len(genes) * n_groups)
    nexp = np.zeros(len(genes) * n_groups)
    reader = pd.read_csv(SC / f"gene_sorted-{comp}.matrix.mtx", sep=" ", comment="%", header=None, skiprows=2,
                         names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.float32},
                         chunksize=chunksize)
    for chunk in reader:
        g = chunk["g"].to_numpy() - 1
        c = chunk["c"].to_numpy() - 1
        idx = g.astype(np.int64) * n_groups + codes[c]
        sums += np.bincount(idx, weights=chunk["v"].to_numpy() * scale[c], minlength=sums.size)
        nexp += np.bincount(idx, minlength=nexp.size)
    ncell = np.bincount(codes, minlength=n_groups)
    gi, gr = np.nonzero(sums.reshape(len(genes), n_groups))
    flat = gi * n_groups + gr
    parts = pd.Series(groups).str.split("|", expand=True)
    out = pd.DataFrame({"gene": genes.to_numpy()[gi], "subject": parts[0].to_numpy()[gr], "health": parts[1].to_numpy()[gr],
                        "cluster": parts[2].to_numpy()[gr], "compartment": COMPARTMENTS[comp],
                        "cp10k_sum": sums[flat].astype(np.float32), "n_expressing": nexp[flat].astype(np.int32),
                        "n_cells": ncell[gr].astype(np.int32)})
    return out


def build_pseudobulk() -> pd.DataFrame:
    m = meta()
    pb = pd.concat([_pseudobulk(c, m) for c in COMPARTMENTS], ignore_index=True)
    pb.to_parquet(PROCESSED / "sc_pseudobulk.parquet", index=False)
    cells = m.groupby(["Subject", "health", "Cluster"]).size().rename("n_cells").reset_index()
    cells.to_parquet(PROCESSED / "sc_cells.parquet", index=False)
    return pb


def celltype_share(pb: pd.DataFrame) -> pd.DataFrame:
    """Share of each gene's transcripts per cluster within a health state, and the cluster mean CP10K."""
    g = pb.groupby(["gene", "health", "cluster", "compartment"], as_index=False)[["cp10k_sum", "n_expressing"]].sum()
    cells = pb.drop_duplicates(["subject", "health", "cluster"]).groupby(["health", "cluster"])["n_cells"].sum()
    g["n_cells"] = [cells.get((h, c), np.nan) for h, c in zip(g["health"], g["cluster"])]
    g["mean_cp10k"] = g["cp10k_sum"] / g["n_cells"]
    g["pct_expressing"] = 100 * g["n_expressing"] / g["n_cells"]
    g["share"] = g["cp10k_sum"] / g.groupby(["gene", "health"])["cp10k_sum"].transform("sum")
    g.to_parquet(PROCESSED / "sc_celltype_share.parquet", index=False)
    return g


def celltype_de(pb: pd.DataFrame, min_donors: int = 3) -> pd.DataFrame:
    """Per gene and cluster: inflamed vs healthy, one value per donor (mean CP10K of the donor's cells, log2(x+1))."""
    from statsmodels.stats.multitest import multipletests
    pb = pb[pb["health"].isin(["inflamed", "healthy"])].copy()
    pb["mean"] = np.log2(pb["cp10k_sum"] / pb["n_cells"] + 1)
    # donors that have the cluster but no reads for the gene count as 0
    donors = pb.drop_duplicates(["subject", "health", "cluster"])[["subject", "health", "cluster"]]
    rows = []
    for (gene, cl), d in pb.groupby(["gene", "cluster"]):
        dn = donors[donors["cluster"] == cl]
        v = dn.merge(d[["subject", "health", "mean"]], on=["subject", "health"], how="left").fillna({"mean": 0.0})
        a, b = v.loc[v.health == "inflamed", "mean"], v.loc[v.health == "healthy", "mean"]
        if len(a) < min_donors or len(b) < min_donors:
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append((gene, cl, a.mean() - b.mean(), p, len(a), len(b)))
    de = pd.DataFrame(rows, columns=["gene", "cluster", "log2fc", "p", "n_inflamed", "n_healthy"])
    ok = de["p"].notna()
    de["padj"] = np.nan
    de.loc[ok, "padj"] = multipletests(de.loc[ok, "p"], method="fdr_bh")[1]
    de.to_parquet(PROCESSED / "sc_celltype_de.parquet", index=False)
    return de


def run(genes: list[str] | None = None) -> dict:
    pb = build_pseudobulk() if not (PROCESSED / "sc_pseudobulk.parquet").exists() else pd.read_parquet(PROCESSED / "sc_pseudobulk.parquet")
    sh = celltype_share(pb)
    if genes is not None:
        pb = pb[pb["gene"].isin(genes)]
    de = celltype_de(pb)
    return {"pseudobulk_rows": len(pb), "genes": pb["gene"].nunique(), "clusters": pb["cluster"].nunique(),
            "share_rows": len(sh), "de_tests": len(de)}


if __name__ == "__main__":
    import json
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    print(json.dumps(run(sorted(set(m["gene"]))), indent=1))
