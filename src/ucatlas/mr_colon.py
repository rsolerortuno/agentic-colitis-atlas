"""Mendelian randomisation with gut-tissue eQTL (added 2026-11-10; closes the "plasma, not gut" limit of mr.py).

The plasma pQTL analysis (mr.py) asks whether genetically higher blood protein changes UC risk. Most candidates act in
the colon wall, so the same design is repeated with instruments for gene expression measured in the gut itself:

- colon tissue: GTEx v10 gene-level eQTL of transverse and sigmoid colon (497 and 430 donors, no IBD) and CEDAR
  microarray eQTL of transverse colon and rectum biopsies (about 300 healthy volunteers);
- gut cell types: IBDverse single-cell eQTL (OTAR2057, gut biopsies of people with and without IBD; about 60 epithelial,
  immune and stromal cell types).

Instruments, outcome and estimate are those of mr.py (eQTL credible sets have no z-score, so the instrument standard
error is |beta| / z with z from the p-value): Open Targets 26.09 fine-mapped cis credible sets (lead variant,
p < 5e-8, MHC and palindromic variants out, a tag of the same 95 % credible set when the lead is not in the GWAS), the
de Lange 2017 UC GWAS, Wald ratio or IVW. Effects are per SD of normalised expression. One estimate per gene and
study; per gene, the colon-tissue layer keeps the study with the most instruments (then the strongest), the cell layer
keeps every cell type and reports the strongest. Benjamini-Hochberg within each layer across genes. Significant genes
are checked for colocalisation with an IBD GWAS credible set in Open Targets, as in mr.py.

Outputs: data/raw/extra/ot_gut_eqtl_credible_sets.parquet, data/processed/mr_colon.parquet (one row per gene and
layer), mr_colon_celltypes.parquet (every gene x cell type estimate), mr_colon_instruments.parquet.
"""
from __future__ import annotations

import io
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from . import mr
from .paths import PROCESSED, RAW

EXTRACT = RAW / "extra" / "ot_gut_eqtl_credible_sets.parquet"
CACHE = EXTRACT.parent / "_gut_eqtl_parts"
COLON = {"UBERON_0001157": "transverse colon", "UBERON_0001159": "sigmoid colon", "UBERON_0001052": "rectum"}


def studies() -> pd.DataFrame:
    """The gut eQTL studies used, with gene symbol, layer and a readable tissue or cell label."""
    st = pd.read_parquet(RAW / "extra" / "ot_study.parquet",
                         columns=["studyId", "projectId", "studyType", "biosampleId", "geneId", "nSamples"])
    gtex = (st.projectId == "GTEx-v10") & (st.studyType == "eqtl") & st.biosampleId.isin(COLON) \
        & st.studyId.str.startswith("gtex-v10_ge_")
    cedar = (st.projectId == "CEDAR") & (st.studyType == "eqtl") & st.biosampleId.isin(COLON)
    ibd = (st.projectId == "OTAR2057_IBDverse") & (st.studyType == "sceqtl")
    s = st[gtex | cedar | ibd].copy()
    s["layer"] = np.where(s.projectId == "OTAR2057_IBDverse", "gut_cell", "colon_tissue")
    s["label"] = np.where(s.layer == "gut_cell", s.studyId.str.extract(r"IBDverse_ge_(.*)_ENSG", expand=False),
                          s.projectId.str.replace("-v10", "") + " " + s.biosampleId.map(COLON).fillna(""))
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    s["gene"] = s["geneId"].map(dict(zip(hg["ensembl_gene_id"], hg["symbol"])))
    return s.dropna(subset=["gene"])


def _filter_part(name: str, ids) -> pd.DataFrame:
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    cache = CACHE / name
    if cache.exists():
        return pd.read_parquet(cache)
    r = requests.get(f"{mr.OT}/credible_set/{name}", timeout=900)
    r.raise_for_status()
    tab = pq.read_table(io.BytesIO(r.content), columns=mr.COLS)
    tab = tab.filter(pc.is_in(tab["studyId"], value_set=ids))
    t = tab.to_pandas()
    del tab
    t = t[t["pValueMantissa"].astype(float) * 10.0 ** t["pValueExponent"].astype(float) < mr.P_INSTRUMENT]
    t["locus"] = t["locus"].map(lambda xs: json.dumps(
        [[x["variantId"], round(float(x["posteriorProbability"] or 0), 4), x["beta"]]
         for x in (xs if xs is not None else []) if x["is95CredibleSet"] and x["beta"] is not None]))
    CACHE.mkdir(parents=True, exist_ok=True)
    t.to_parquet(cache, index=False)
    return t


def download(workers: int = 6) -> dict:
    """Stream the credible-set parts and keep the significant gut eQTL credible sets (parts cached, resumable)."""
    import pyarrow as pa
    from concurrent.futures import ThreadPoolExecutor
    from .download import sha256
    if EXTRACT.exists():
        return {"rows": len(pd.read_parquet(EXTRACT, columns=["studyLocusId"]))}
    ids = pa.array(sorted(studies()["studyId"]))
    parts = mr._parts("credible_set")
    keep = []
    with ThreadPoolExecutor(workers) as ex:
        for i, t in enumerate(ex.map(lambda n: _filter_part(n, ids), parts)):
            keep.append(t)
            if i % 20 == 0:
                print(f"part {i + 1}/{len(parts)}: {sum(len(k) for k in keep)} gut eQTL credible sets", flush=True)
    out = pd.concat(keep, ignore_index=True).sort_values("studyLocusId").reset_index(drop=True)
    out.to_parquet(EXTRACT, index=False)
    import shutil
    shutil.rmtree(CACHE, ignore_errors=True)
    mpath = RAW / "manifest.json"
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    manifest["extra/ot_gut_eqtl_credible_sets.parquet"] = {
        "url": f"{mr.OT}/credible_set/ ({len(parts)} parts, filtered to GTEx-v10 colon gene eQTL, CEDAR colon/rectum "
               f"eQTL and OTAR2057 IBDverse sc-eQTL studies, p < 5e-8)",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256(EXTRACT), "bytes": EXTRACT.stat().st_size}
    mpath.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    return {"rows": len(out)}


def instruments() -> pd.DataFrame:
    cs = pd.read_parquet(EXTRACT)
    cs["p"] = cs["pValueMantissa"].astype(float) * 10.0 ** cs["pValueExponent"].astype(float)
    cs = cs[(~cs["isTransQtl"].fillna(False)) & (cs["p"] < mr.P_INSTRUMENT)].copy()
    st = studies().set_index("studyId")
    cs = cs[cs["studyId"].isin(st.index)]
    cs["gene"] = cs["studyId"].map(st["gene"])
    # eQTL credible sets carry no z-score: |z| from the p-value (mantissa and exponent, so no underflow)
    from scipy.special import ndtri_exp
    logp = np.log(cs["pValueMantissa"].astype(float)) + cs["pValueExponent"].astype(float) * np.log(10)
    cs["zScore"] = cs["zScore"].fillna(pd.Series(-ndtri_exp(logp - np.log(2)), index=cs.index))
    cs = cs.dropna(subset=["gene", "beta", "zScore"])
    cs = cs[cs["zScore"] != 0]
    cs["se"] = (cs["beta"] / cs["zScore"]).abs()
    mhc = (cs["chromosome"] == mr.MHC[0]) & cs["position"].between(mr.MHC[1], mr.MHC[2])
    return cs[~mhc]


def _per_study(h: pd.DataFrame) -> pd.DataFrame:
    """mr.estimate without its per-gene collapse: one row per gene and study."""
    from scipy import stats
    out = []
    for (gene, study), d in h.groupby(["gene", "studyId"]):
        ratio = d["beta_outcome"] / d["beta_exposure"]
        se = d["se_outcome"] / d["beta_exposure"].abs()
        w = 1 / se ** 2
        b = float((w * ratio).sum() / w.sum())
        s = float(np.sqrt(1 / w.sum()))
        out.append({"gene": gene, "studyId": study, "n_instruments": len(d),
                    "method": "IVW" if len(d) > 1 else "Wald ratio", "beta": b, "se": s,
                    "p": float(2 * stats.norm.sf(abs(b / s))), "or_per_sd": float(np.exp(b)),
                    "min_p_exposure": float(d["p_exposure"].min()),
                    "lead_instrument": d.sort_values("p_exposure").iloc[0]["lead_variant"],
                    "uses_tag_variant": bool(d["substitute_pp"].notna().any())})
    return pd.DataFrame(out)


def run(q_coloc: float = 0.05) -> dict:
    from statsmodels.stats.multitest import multipletests
    cs = instruments()
    g = mr.outcome(mr._variant_keys(cs))
    h = mr.harmonise(cs, g)
    st = studies().set_index("studyId")
    est = _per_study(h)
    est["layer"] = est["studyId"].map(st["layer"])
    est["label"] = est["studyId"].map(st["label"])
    cells = est[est.layer == "gut_cell"].copy()
    if len(cells):
        cells["q_within_gene"] = cells.groupby("gene")["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    cells.to_parquet(PROCESSED / "mr_colon_celltypes.parquet", index=False)
    rows = []
    for layer, d in est.groupby("layer"):
        if layer == "colon_tissue":
            best = d.sort_values(["gene", "n_instruments", "min_p_exposure"], ascending=[True, False, True])
        else:  # the strongest cell type, with a Bonferroni factor for the number of cell types tested for that gene
            d = d.assign(n_celltypes=d.groupby("gene")["studyId"].transform("size"))
            d["p_raw"] = d["p"]
            d["p"] = (d["p"] * d["n_celltypes"]).clip(upper=1)
            best = d.sort_values(["gene", "p"])
        best = best.drop_duplicates("gene").copy()
        best["q"] = multipletests(best["p"], method="fdr_bh")[1]
        if layer == "colon_tissue":
            agree = d.groupby("gene")["beta"].agg(lambda b: int((np.sign(b) == np.sign(b.iloc[0])).sum()))
            best["studies"] = best["gene"].map(d.groupby("gene").size())
        rows.append(best)
    r = pd.concat(rows, ignore_index=True)
    parts = []
    for layer, d in r.groupby("layer"):   # FinnGen replication and combined estimate, BH within each layer
        f = d["n_celltypes"] if layer == "gut_cell" and "n_celltypes" in d else None
        parts.append(mr.attach_liu2023(mr.attach_finngen(d.reset_index(drop=True), cs,
                                                         f.reset_index(drop=True) if f is not None else None), cs))
    r = pd.concat(parts, ignore_index=True)
    sig = r[r["q"] < q_coloc]
    loci = h.merge(sig[["gene", "studyId"]], on=["gene", "studyId"])
    cl = mr.coloc(sorted(loci["studyLocusId"].unique()))
    best = {}
    for row in loci.itertuples():
        c = cl.get(row.studyLocusId) or {}
        b = best.setdefault((row.gene, row.studyId), {"coloc_h4": None, "coloc_clpp": None, "coloc_trait": None})
        for k_src, k in (("h4", "coloc_h4"), ("clpp", "coloc_clpp")):
            if c.get(k_src) is not None and (b[k] is None or c[k_src] > b[k]):
                b[k] = round(float(c[k_src]), 4)
                b["coloc_trait"] = c.get("trait")
    for k in ("coloc_h4", "coloc_clpp", "coloc_trait"):
        r[k] = [best.get((g_, s_), {}).get(k) for g_, s_ in zip(r["gene"], r["studyId"])]
    yes = (r["coloc_h4"].astype(float).fillna(0) >= 0.8) | (r["coloc_clpp"].astype(float).fillna(0) >= 0.01)
    r["colocalises"] = pd.Series(np.where(r["q"] < q_coloc, yes, None), index=r.index, dtype=object)
    # agreement with the plasma pQTL estimate where both exist
    pq_ = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").set_index("gene") if (PROCESSED / "mr_pqtl.parquet").exists() else None
    if pq_ is not None:
        r["plasma_beta"] = r["gene"].map(pq_["beta"])
        r["plasma_q"] = r["gene"].map(pq_["q"])
    r = r.sort_values(["layer", "p"]).reset_index(drop=True)
    r.to_parquet(PROCESSED / "mr_colon.parquet", index=False)
    h.to_parquet(PROCESSED / "mr_colon_instruments.parquet", index=False)
    out = {}
    for layer, d in r.groupby("layer"):
        out[layer] = {"genes_tested": len(d), "significant_q05": int((d.q < 0.05).sum()),
                      "colocalising": int(d.colocalises.fillna(False).astype(bool).sum())}
    out["instruments"] = len(h)
    return out


def load(store=None) -> int:
    """One mendelian_randomization_gut record per gene and layer (significant or not)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    r = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    n = 0
    names = {"colon_tissue": "colon tissue eQTL (GTEx v10 transverse/sigmoid colon, CEDAR colon and rectum)",
             "gut_cell": "gut cell-type eQTL (IBDverse single-cell, IBD and control biopsies)"}
    with store.batch():
        store.add_node("DATASET:MR_GUT_EQTL_UC", "Dataset",
                       "Mendelian randomisation: gut eQTL (GTEx v10 colon, CEDAR, IBDverse; Open Targets 26.09) on UC GWAS")
        for x in r.itertuples():
            node = store.resolve(x.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": x.gene, "exposure": f"expression in {names[x.layer]}, per SD", "layer": x.layer,
                   "tissue_or_cell": x.label, "outcome": "ulcerative colitis", "method": x.method,
                   "n_instruments": int(x.n_instruments), "beta": round(float(x.beta), 4), "se": round(float(x.se), 4),
                   "odds_ratio_per_sd": round(float(x.or_per_sd), 3), "p": float(f"{x.p:.3g}"),
                   "q": float(f"{x.q:.3g}"), "significant": bool(x.q < 0.05),
                   "direction": "higher expression raises UC risk" if x.beta > 0 else "higher expression lowers UC risk",
                   "lead_instrument": x.lead_instrument,
                   "colocalises_with_ibd_gwas": None if x.colocalises is None or pd.isna(x.colocalises) else bool(x.colocalises),
                   "colocalisation_h4": None if pd.isna(x.coloc_h4) else float(x.coloc_h4),
                   "colocalisation_clpp": None if pd.isna(x.coloc_clpp) else float(x.coloc_clpp),
                   "study": x.studyId}
            if pd.notna(getattr(x, "beta_finngen", np.nan)):
                rec.update({"finngen_odds_ratio_per_sd": round(float(np.exp(x.beta_finngen)), 3),
                            "finngen_p": float(f"{x.p_finngen:.3g}"),
                            "replicates_in_finngen": bool(x.replicates_finngen),
                            "combined_odds_ratio_per_sd": round(float(np.exp(x.beta_combined)), 3),
                            "combined_q": float(f"{x.q_combined:.3g}")})
            if pd.notna(getattr(x, "beta_liu2023", np.nan)):
                rec.update({"liu2023_odds_ratio_per_sd": round(float(np.exp(x.beta_liu2023)), 3),
                            "liu2023_p": float(f"{x.p_liu2023:.3g}"), "replicates_in_liu2023": bool(x.replicates_liu2023)})
            if x.layer == "gut_cell":
                rec["cell_types_tested"] = int(x.n_celltypes)
                rec["p_strongest_cell_type_uncorrected"] = float(f"{x.p_raw:.3g}")
            ev = store.add_evidence(evidence_id("mendelian_randomization_gut", f"{x.gene}|{x.layer}"),
                                    "mendelian_randomization_gut", "DATASET:MR_GUT_EQTL_UC", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "causal_estimate_for", "DATASET:MR_GUT_EQTL_UC", ev)
            n += 1
    return n


if __name__ == "__main__":
    import sys
    if "download" in sys.argv:
        print(download())
    else:
        print(json.dumps(run(), indent=1))
        r = pd.read_parquet(PROCESSED / "mr_colon.parquet")
        print(r[r.q < 0.05].head(40).to_string())
