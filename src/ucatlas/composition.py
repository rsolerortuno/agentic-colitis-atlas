"""M12: cell composition, the equivalent of tumour purity. Is a bulk change of inflamed UC mucosa a change inside
the tissue, or does it come from cells that flood in (neutrophils, monocytes, plasma cells)?

1. Lineage markers from Gubatan 2026 (SCP3755, includes neutrophils): per lineage, the 25 genes with the largest
   share of their transcripts in that lineage (>= 80 %) among genes expressed in >= 20 % of its cells.
2. Per MSCCR patient-level sample (inflamed UC and control colon): lineage score = mean z-score of log2 expression of
   its markers (the gene under test is left out of its own lineage's score).
3. Per gene: ordinary least squares, expression ~ inflamed + lineage scores; the inflamed coefficient is the change
   adjusted for composition. Verdict: 'intrinsic' when the adjusted change stays significant (BH padj < 0.05) in the
   same direction and keeps at least half of its size, 'partly_composition' when it stays significant but loses more
   than half, 'composition' when it is no longer significant. (The 50 % rule was added after the first run showed
   CXCR2 keeping significance with 17 % of its change; logged in evals/history.md.)
Output: data/processed/composition_effect.parquet and composition_markers.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED

LINEAGES = {"Neutrophil": "neutrophil", "Macrophage and DC": "myeloid", "Mast Cell": "mast", "Plasma": "plasma",
            "B Cell": "B", "T Cell and NK": "T/NK", "Stromal": "stromal", "Endothelial": "endothelial",
            "Enterocytes": "epithelial", "Goblet": "epithelial", "Stem Cycling": "epithelial"}
N_MARKERS = 25


def markers() -> dict[str, list[str]]:
    s = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    s = s[s["health"] == "inflamed"].copy()
    s["lineage"] = s["compartment"].map(LINEAGES)
    s = s.dropna(subset=["lineage"])
    g = s.groupby(["gene", "lineage"]).agg(cp=("cp10k_sum", "sum"), nexp=("n_expressing", "sum")).reset_index()
    cells = s.drop_duplicates(["cluster"]).groupby("lineage")["n_cells"].sum()
    g["share"] = g["cp"] / g.groupby("gene")["cp"].transform("sum")
    g["pct"] = g["nexp"] / g["lineage"].map(cells)
    out = {}
    for lin, d in g[(g["share"] >= 0.8) & (g["pct"] >= 0.2)].groupby("lineage"):
        out[lin] = d.sort_values("cp", ascending=False)["gene"].head(N_MARKERS).tolist()
    return out


def _patient_matrix() -> tuple[pd.DataFrame, pd.Series]:
    from .expression import patient_means
    x = pd.read_parquet(PROCESSED / "expr_GSE193677.parquet")
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    s = s[(s["cohort"] == "GSE193677") & (s["timepoint"] == "baseline") & s["group"].isin(["UC_inflamed", "healthy"])]
    pm = patient_means(x, s)
    infl = pd.Series([c.endswith("|UC_inflamed") for c in pm.columns], index=pm.columns).astype(float)
    return pm, infl


def _zscores(pm: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    sub = pm.reindex([g for g in genes if g in pm.index])
    return sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)


def run(genes: list[str]) -> pd.DataFrame:
    mk = markers()
    pm, infl = _patient_matrix()
    z = {lin: _zscores(pm, g) for lin, g in mk.items()}
    base_scores = pd.DataFrame({lin: zz.mean(axis=0) for lin, zz in z.items()})
    rows = []
    for gene in genes:
        if gene not in pm.index:
            continue
        sc = base_scores.copy()
        for lin, zz in z.items():
            if gene in zz.index:  # leave the gene out of its own lineage score
                sc[lin] = zz.drop(index=gene).mean(axis=0)
        y = pm.loc[gene].to_numpy(float)
        X0 = np.column_stack([np.ones(len(y)), infl.to_numpy()])
        X = np.column_stack([X0, sc.to_numpy()])
        out = {"gene": gene}
        for name, M in (("unadjusted", X0), ("adjusted", X)):
            beta, *_ = np.linalg.lstsq(M, y, rcond=None)
            resid = y - M @ beta
            dof = len(y) - M.shape[1]
            sigma2 = resid @ resid / dof
            cov = sigma2 * np.linalg.pinv(M.T @ M)
            se = float(np.sqrt(cov[1, 1]))
            from scipy import stats
            t = beta[1] / se
            out[f"{name}_log2fc"] = float(beta[1])
            out[f"{name}_se"] = se
            out[f"{name}_p"] = float(2 * stats.t.sf(abs(t), dof))
        rows.append(out)
    d = pd.DataFrame(rows)
    for name in ("unadjusted", "adjusted"):
        d[f"{name}_padj"] = multipletests(d[f"{name}_p"], method="fdr_bh")[1]
    same = np.sign(d["adjusted_log2fc"]) == np.sign(d["unadjusted_log2fc"])
    # share of the bulk change that remains after adjustment (1 = untouched, 0 = fully explained by composition)
    d["retained"] = (d["adjusted_log2fc"] / d["unadjusted_log2fc"]).clip(lower=0, upper=1.5)
    sig = (d["adjusted_padj"] < 0.05) & same
    d["verdict"] = np.select([sig & (d["retained"] >= 0.5), sig], ["intrinsic", "partly_composition"], "composition")
    d["n_inflamed"], d["n_control"] = int(infl.sum()), int((1 - infl).sum())
    d.to_parquet(PROCESSED / "composition_effect.parquet", index=False)
    (PROCESSED / "composition_markers.json").write_text(json.dumps(mk, indent=1))
    return d


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "composition_effect.parquet")
    n = 0
    with store.batch():
        store.add_node("ANALYSIS:composition", "Dataset", "Composition-adjusted change, MSCCR inflamed UC vs control "
                                                          "(lineage scores from Gubatan 2026 markers)")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            ev = store.add_evidence(evidence_id("composition_effect", r.gene), "composition_effect", "ANALYSIS:composition",
                                    date.today().isoformat(), gene=r.gene, layer="rna",
                                    unadjusted_log2fc=round(r.unadjusted_log2fc, 4), unadjusted_padj=float(f"{r.unadjusted_padj:.3g}"),
                                    adjusted_log2fc=round(r.adjusted_log2fc, 4), adjusted_padj=float(f"{r.adjusted_padj:.3g}"),
                                    retained_fraction=round(float(r.retained), 3), verdict=r.verdict, n_inflamed=int(r.n_inflamed), n_control=int(r.n_control))
            store.add_claim(node["id"], "composition_assessed_in", "ANALYSIS:composition", ev)
            n += 1
    return n


if __name__ == "__main__":
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    d = run(sorted(m["gene"]))
    print(d["verdict"].value_counts().to_dict(), "records:", load())
