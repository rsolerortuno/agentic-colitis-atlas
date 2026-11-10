"""M2: differential expression per cohort and random-effects meta-analysis (deterministic, no LLM).

Every cohort is analysed on its own scale (RNA-seq counts or log2 microarray intensities); only effect sizes and
their standard errors are combined. Patients with several biopsies in one group are averaged first, so a patient
counts once (biopsies of one person are not independent).

Contrasts
- GSE193677 (RNA-seq): UC inflamed vs control, UC non-inflamed vs control, and inflamed vs non-inflamed within the
  same patient (paired; 135 patients).
- Microarray cohorts: UC at baseline vs healthy controls.
- RNA-seq replication cohorts (added 2026-11-10): PROTECT (GSE109142) and RISK (GSE117993), paediatric rectum at
  diagnosis, and GSE165512, adult colon; UC vs healthy controls, NCBI-generated raw counts. They let genes that are not
  on the arrays (mostly lncRNAs and newer genes) meet more than one cohort.
Outputs: data/processed/de_results.parquet, meta.parquet, symbol_map.parquet, expr_<cohort>.parquet (gene x sample).
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .literature_entities import symbol_lookup
from .paths import PROCESSED, RAW

ARRAYS = {"GSE206285": "GPL13158", "GSE92415": "GPL13158", "GSE16879": "GPL570", "GSE73661": "GPL6244"}
RNASEQ = ["GSE109142", "GSE117993", "GSE165512"]
# treatment-response cohorts without controls (added 2026-11-10, third round): expression matrices only, not in the
# UC vs control meta-analysis
TREATMENT_ARRAYS = {"GSE23597": "GPL570"}
TREATMENT_RNASEQ = ["GSE72819"]
COHORTS = ["GSE193677", *ARRAYS, *RNASEQ]


# ------------------------------------------------------------------ gene identifiers
def _hgnc() -> pd.DataFrame:
    return pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id", "locus_group"],
                       dtype=str)


def ensembl_to_symbol() -> dict[str, str]:
    h = _hgnc().dropna(subset=["ensembl_gene_id"])
    return dict(zip(h["ensembl_gene_id"], h["symbol"]))


def resolve_symbol(sym: str) -> str | None:
    """Approved HGNC symbol for a platform symbol (approved, previous or unique alias, e.g. IL8 -> CXCL8)."""
    if not isinstance(sym, str) or not sym or "///" in sym:
        return None
    return symbol_lookup().get(sym.strip().upper())


# ------------------------------------------------------------------ loading
def series_table(gse: str) -> pd.DataFrame:
    with gzip.open(RAW / gse / f"{gse}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
        df = pd.read_csv(f, sep="\t", index_col=0, comment="!", low_memory=False)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def platform_symbols(gpl: str) -> pd.Series:
    with gzip.open(RAW / "platforms" / f"{gpl}.annot.gz", "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        a = pd.read_csv(f, sep="\t", usecols=["ID", "Gene symbol"], dtype=str, comment="!")
    a = a.dropna()
    return a.set_index("ID")["Gene symbol"]


def array_genes(gse: str) -> tuple[pd.DataFrame, dict[str, str]]:
    """Probe matrix -> gene matrix (log2). Per gene the probe with the highest mean expression is kept."""
    x = series_table(gse)
    if np.nanmax(x.values) > 100:  # not yet on log2 scale
        x = np.log2(x.clip(lower=1))
    sym = platform_symbols({**ARRAYS, **TREATMENT_ARRAYS}[gse]).reindex(x.index)
    approved = sym.map(resolve_symbol)
    smap = {s: a for s, a in zip(sym, approved) if isinstance(s, str) and isinstance(a, str)}
    x = x[approved.notna()]
    x["gene"] = approved[approved.notna()]
    x["_mean"] = x.drop(columns="gene").mean(axis=1)
    x = x.sort_values("_mean", ascending=False).drop_duplicates("gene").drop(columns="_mean").set_index("gene")
    return x, smap


def rnaseq_genes() -> pd.DataFrame:
    """MSCCR counts -> log2 normalised counts (median-of-ratios size factors), protein-coding and lncRNA HGNC genes."""
    c = pd.read_csv(RAW / "GSE193677" / "GSE193677_MSCCR_Biopsy_counts.txt.gz", sep=r"\s+", index_col=0)
    c.index = c.index.str.replace('"', "").str.split(".").str[0]
    c.columns = c.columns.str.replace('"', "")
    sym = c.index.map(ensembl_to_symbol())
    c = c[sym.notna()]
    c.index = sym[sym.notna()]
    return _normalise(c.groupby(level=0).sum())


def ncbi_counts(gse: str) -> pd.DataFrame:
    """NCBI-generated raw counts (Entrez GeneID rows, GSM columns) -> log2 normalised counts on HGNC symbols."""
    c = pd.read_csv(RAW / gse / f"{gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="	", index_col=0)
    # every HGNC locus group, as for the MSCCR Ensembl ids
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="	", usecols=["symbol", "entrez_id"], dtype=str).dropna()
    sym = c.index.astype(str).map(dict(zip(hg["entrez_id"], hg["symbol"])))
    c = c[sym.notna()]
    c.index = sym[sym.notna()]
    return _normalise(c.groupby(level=0).sum())


def _normalise(c: pd.DataFrame) -> pd.DataFrame:
    c = c[(c >= 10).sum(axis=1) >= 0.1 * c.shape[1]]  # expressed in at least 10 % of biopsies
    logc = np.log(c.replace(0, np.nan))
    geo_mean = logc.mean(axis=1)
    ok = geo_mean.notna() & np.isfinite(geo_mean) & c.gt(0).all(axis=1)
    sf = np.exp((logc[ok].sub(geo_mean[ok], axis=0)).median(axis=0))
    return np.log2(c.div(sf, axis=1) + 1)


# ------------------------------------------------------------------ statistics
def _bh(p: pd.Series) -> np.ndarray:
    out = np.full(len(p), np.nan)
    ok = p.notna().values
    if ok.any():
        out[ok] = multipletests(p.values[ok], method="fdr_bh")[1]
    return out


def welch(x: pd.DataFrame, a: list[str], b: list[str]) -> pd.DataFrame:
    """Group a vs group b per gene: difference of means (log2FC), Welch SE, t-test p."""
    xa, xb = x[a], x[b]
    ma, mb = xa.mean(axis=1), xb.mean(axis=1)
    va, vb = xa.var(axis=1, ddof=1), xb.var(axis=1, ddof=1)
    se = np.sqrt(va / len(a) + vb / len(b))
    t, p = stats.ttest_ind(xa, xb, axis=1, equal_var=False, nan_policy="omit")
    return pd.DataFrame({"log2fc": ma - mb, "se": se, "t": t, "p": p, "n1": len(a), "n2": len(b),
                         "mean_disease": ma, "mean_control": mb}, index=x.index)


def paired(x: pd.DataFrame, a: list[str], b: list[str]) -> pd.DataFrame:
    d = x[a].values - x[b].values
    n = d.shape[1]
    md = d.mean(axis=1)
    se = d.std(axis=1, ddof=1) / np.sqrt(n)
    t, p = stats.ttest_1samp(d, 0, axis=1)
    return pd.DataFrame({"log2fc": md, "se": se, "t": t, "p": p, "n1": n, "n2": n,
                         "mean_disease": x[a].mean(axis=1).values, "mean_control": x[b].mean(axis=1).values},
                        index=x.index)


def patient_means(x: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Average the biopsies of one patient within one group (columns named patient|group)."""
    s = samples.set_index("sample").loc[[c for c in x.columns if c in set(samples["sample"])]]
    key = s["patient"].astype(str) + "|" + s["group"]
    return x[s.index].T.groupby(key.values).mean().T


def random_effects(fc: np.ndarray, se: np.ndarray) -> tuple[float, float, float, float, int]:
    """DerSimonian-Laird random-effects meta-analysis: (estimate, se, p, I2, k)."""
    ok = np.isfinite(fc) & np.isfinite(se) & (se > 0)
    fc, se = fc[ok], se[ok]
    k = len(fc)
    if k == 0:
        return np.nan, np.nan, np.nan, np.nan, 0
    w = 1 / se ** 2
    fixed = np.sum(w * fc) / np.sum(w)
    q = np.sum(w * (fc - fixed) ** 2)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    est = np.sum(wr * fc) / np.sum(wr)
    se_est = np.sqrt(1 / np.sum(wr))
    p = 2 * stats.norm.sf(abs(est / se_est))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    return est, se_est, p, i2, k


# ------------------------------------------------------------------ run
def run() -> dict:
    samples = pd.read_parquet(PROCESSED / "samples.parquet")
    base = samples[samples["timepoint"] == "baseline"]
    results, smap_all = [], {}

    # discovery: MSCCR RNA-seq
    x = rnaseq_genes()
    msc = base[base["cohort"] == "GSE193677"]
    pm = patient_means(x, msc)
    cols = lambda g: [c for c in pm.columns if c.endswith("|" + g)]
    infl, noninfl, ctrl = cols("UC_inflamed"), cols("UC_noninflamed"), cols("healthy")
    for name, a in (("inflamed_vs_control", infl), ("noninflamed_vs_control", noninfl)):
        r = welch(pm, a, ctrl).assign(cohort="GSE193677", contrast=name)
        results.append(r)
    both = sorted({c.split("|")[0] for c in infl} & {c.split("|")[0] for c in noninfl})
    r = paired(pm, [f"{p}|UC_inflamed" for p in both], [f"{p}|UC_noninflamed" for p in both])
    results.append(r.assign(cohort="GSE193677", contrast="inflamed_vs_noninflamed_paired"))
    x.astype("float32").to_parquet(PROCESSED / "expr_GSE193677.parquet")

    # replication: microarrays, UC baseline vs healthy (patient = one baseline biopsy)
    for gse in ARRAYS:
        xa, smap = array_genes(gse)
        smap_all.update(smap)
        s = base[base["cohort"] == gse].set_index("geo")
        xa = xa[[c for c in xa.columns if c in s.index]]
        uc = [c for c in xa.columns if s.loc[c, "group"] == "UC"]
        hc = [c for c in xa.columns if s.loc[c, "group"] == "healthy"]
        results.append(welch(xa, uc, hc).assign(cohort=gse, contrast="uc_vs_control"))
        xa_all, _ = array_genes(gse)  # every sample (incl. post-treatment) on the gene level, for later milestones
        xa_all.astype("float32").to_parquet(PROCESSED / f"expr_{gse}.parquet")

    # replication: RNA-seq cohorts, UC vs healthy (one biopsy per patient)
    for gse in RNASEQ:
        xr = ncbi_counts(gse)
        s = base[(base["cohort"] == gse)].set_index("geo")
        uc = [c for c in xr.columns if c in s.index and s.loc[c, "group"] == "UC"]
        hc = [c for c in xr.columns if c in s.index and s.loc[c, "group"] == "healthy"]
        results.append(welch(xr, uc, hc).assign(cohort=gse, contrast="uc_vs_control"))
        xr.astype("float32").to_parquet(PROCESSED / f"expr_{gse}.parquet")

    for gse in TREATMENT_ARRAYS:
        xa, _ = array_genes(gse)
        xa.astype("float32").to_parquet(PROCESSED / f"expr_{gse}.parquet")
    for gse in TREATMENT_RNASEQ:
        ncbi_counts(gse).astype("float32").to_parquet(PROCESSED / f"expr_{gse}.parquet")

    de = pd.concat(results).rename_axis("gene").reset_index()
    de["padj"] = de.groupby(["cohort", "contrast"])["p"].transform(lambda p: _bh(p))
    de.to_parquet(PROCESSED / "de_results.parquet", index=False)
    pd.DataFrame({"dataset_symbol": list(smap_all), "symbol": list(smap_all.values())}).to_parquet(
        PROCESSED / "symbol_map.parquet", index=False)

    # meta-analysis of disease vs control: MSCCR inflamed vs control + the four array and three RNA-seq cohorts
    use = de[de["contrast"].isin(["inflamed_vs_control", "uc_vs_control"])]
    rows = []
    for gene, g in use.groupby("gene"):
        est, se, p, i2, k = random_effects(g["log2fc"].values, g["se"].values)
        sig = g[(g["padj"] < 0.05)]
        rows.append({"gene": gene, "meta_log2fc": est, "meta_se": se, "meta_p": p, "i2": i2, "k": k,
                     "n_sig_up": int((sig["log2fc"] > 0).sum()), "n_sig_down": int((sig["log2fc"] < 0).sum())})
    meta = pd.DataFrame(rows)
    # A pooled estimate needs at least two cohorts (review 2026-11-10): a gene measured in one cohort keeps that
    # cohort's estimate but gets no pooled significance, so it cannot enter the expression stream on one cohort alone.
    meta["meta_padj"] = np.nan
    multi = meta["k"] >= 2
    meta.loc[multi, "meta_padj"] = _bh(meta.loc[multi, "meta_p"])
    meta.to_parquet(PROCESSED / "meta.parquet", index=False)
    replication(de).to_parquet(PROCESSED / "replication.parquet", index=False)
    return {"de": de, "meta": meta}


def replication(de: pd.DataFrame, padj: float = 0.05, min_fc: float = 1.0) -> pd.DataFrame:
    """Share of discovery changes (MSCCR inflamed vs control) that are significant in the same direction elsewhere."""
    d = de[(de.cohort == "GSE193677") & (de.contrast == "inflamed_vs_control")].set_index("gene")
    disc = d[(d.padj < padj) & (d.log2fc.abs() >= min_fc)]
    out = []
    for gse in COHORTS[1:]:
        r = de[(de.cohort == gse)].set_index("gene").reindex(disc.index)
        tested = r["padj"].notna()
        same = tested & (r["padj"] < padj) & (np.sign(r["log2fc"]) == np.sign(disc["log2fc"]))
        rho = stats.spearmanr(disc.loc[tested, "log2fc"], r.loc[tested, "log2fc"]).statistic
        out.append({"cohort": gse, "discovery_changes": len(disc), "tested": int(tested.sum()),
                    "replicated_same_direction": int(same.sum()), "share": same.sum() / tested.sum(), "rho": rho})
    return pd.DataFrame(out)


if __name__ == "__main__":
    out = run()
    print(replication(out["de"]).to_string(index=False))
