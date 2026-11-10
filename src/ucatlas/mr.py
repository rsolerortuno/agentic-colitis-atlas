"""Mendelian randomisation with plasma cis-pQTL (added 2026-11-10; design fixed before measuring, evals/history.md).

Question: does genetically higher plasma protein level change the risk of ulcerative colitis? A protein whose genetic
instruments raise UC risk is a candidate for inhibition; one that lowers risk is a candidate for activation (or a
safety flag for inhibition).

Exposure: UK Biobank Pharma Proteomics Project (UKB-PPP, Olink, European ancestry) fine-mapped pQTL credible sets from
Open Targets 26.09 (SuSiE). Instruments: the lead variant of each cis credible set (isTransQtl false) with
p < 5e-8, outside the MHC (chr6:25-34 Mb), not palindromic (A/T, C/G). Effects are per SD of protein level; the
standard error is |beta / z|.
Outcome: UC GWAS of de Lange et al. 2017 (GWAS Catalog GCST004133, harmonised to GRCh38; log odds). When the lead
variant is missing from the GWAS, the variant of the same 95 % credible set with the highest posterior probability
that the GWAS has is used, with its own pQTL effect (so the alleles stay in phase).
Estimate: Wald ratio (one instrument) or inverse-variance weighted fixed effect (several), first-order standard error.
Benjamini-Hochberg across proteins. Check: colocalisation of the pQTL credible set with a UC or IBD GWAS credible set
in Open Targets (coloc H4 or eCAVIAR CLPP), fetched for the significant proteins only. A significant MR without
colocalisation can come from a nearby causal variant in linkage disequilibrium; the record says which.

Outputs: data/raw/extra/ot_pqtl_credible_sets.parquet (filtered extract, with its source parts in the manifest),
data/processed/mr_pqtl.parquet, mr_instruments.parquet.
"""
from __future__ import annotations

import gzip
import io
import json
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests
from scipy import stats

from .paths import PROCESSED, RAW

OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"
GWAS = RAW / "gwas" / "GCST004133_uc_delange2017.h.tsv.gz"
# replication outcome (added 2026-11-10, third round): FinnGen R12 UC (strict definition), GRCh38, beta for ALT,
# no sample overlap with the de Lange 2017 meta-analysis
FINNGEN = RAW / "gwas" / "finngen_R12_K11_UC_STRICT2.gz"
# largest UC GWAS (added 2026-11-10, fourth round): Liu et al. 2023, East Asian + European meta-analysis (23,252
# cases). Its European part includes the IIBDGC samples of de Lange 2017 and FinnGen, so it is reported as a separate
# outcome and not combined with the other two.
LIU2023 = RAW / "gwas" / "GCST90446794_uc_liu2023.h.tsv.gz"
EXTRACT = RAW / "extra" / "ot_pqtl_credible_sets.parquet"
PROJECT = "UKB_PPP_EUR"
P_INSTRUMENT = 5e-8
MHC = ("6", 25_000_000, 34_000_000)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
COLS = ["studyLocusId", "studyId", "variantId", "chromosome", "position", "beta", "zScore", "pValueMantissa",
        "pValueExponent", "credibleSetIndex", "finemappingMethod", "studyType", "isTransQtl", "locus"]
API = "https://api.platform.opentargets.org/api/v4/graphql"


# ------------------------------------------------------------------ download
def _parts(dataset: str) -> list[str]:
    html = requests.get(f"{OT}/{dataset}/", timeout=60).text
    return sorted(set(re.findall(r'href="(part-[^"]+\.parquet)"', html)))


def _filter_part(name: str) -> pd.DataFrame:
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    cache = EXTRACT.parent / "_pqtl_parts" / name
    if cache.exists():
        return pd.read_parquet(cache)
    r = requests.get(f"{OT}/credible_set/{name}", timeout=600)
    r.raise_for_status()
    tab = pq.read_table(io.BytesIO(r.content), columns=COLS)
    mask = pc.and_(pc.equal(tab["studyType"], "pqtl"), pc.starts_with(tab["studyId"], PROJECT))
    t = tab.filter(mask).to_pandas()   # filter before converting: the locus lists are large
    del tab
    t["locus"] = t["locus"].map(lambda xs: json.dumps(
        [[x["variantId"], round(float(x["posteriorProbability"] or 0), 4), x["beta"]]
         for x in (xs if xs is not None else []) if x["is95CredibleSet"] and x["beta"] is not None]))
    cache.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(cache, index=False)
    return t


def download(workers: int = 8) -> dict:
    """Stream the 200 credible-set parts (about 1.8 GB) and keep only the UKB-PPP pQTL rows (parts cached, resumable)."""
    from concurrent.futures import ThreadPoolExecutor
    from .download import sha256
    if EXTRACT.exists():
        return {"rows": len(pd.read_parquet(EXTRACT, columns=["studyLocusId"]))}
    parts = _parts("credible_set")
    with ThreadPoolExecutor(workers) as ex:
        keep = []
        for i, t in enumerate(ex.map(_filter_part, parts)):
            keep.append(t)
            if i % 20 == 0:
                print(f"part {i + 1}/{len(parts)}: {sum(len(k) for k in keep)} pQTL credible sets", flush=True)
    out = pd.concat(keep, ignore_index=True).sort_values("studyLocusId").reset_index(drop=True)
    EXTRACT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(EXTRACT, index=False)
    import shutil
    shutil.rmtree(EXTRACT.parent / "_pqtl_parts", ignore_errors=True)
    mpath = RAW / "manifest.json"
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    manifest["extra/ot_pqtl_credible_sets.parquet"] = {
        "url": f"{OT}/credible_set/ ({len(parts)} parts, filtered to studyType pqtl and studyId {PROJECT}*)",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256(EXTRACT), "bytes": EXTRACT.stat().st_size}
    mpath.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    return {"rows": len(out)}


# ------------------------------------------------------------------ instruments and outcome
def _study_genes() -> pd.Series:
    st = pd.read_parquet(RAW / "extra" / "ot_study.parquet", columns=["studyId", "geneId", "studyType"])
    st = st[st["studyId"].str.startswith(PROJECT)]
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    sym = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    return st.set_index("studyId")["geneId"].map(sym)


def _palindromic(a: str, b: str) -> bool:
    return COMP.get(a) == b


def instruments() -> pd.DataFrame:
    cs = pd.read_parquet(EXTRACT)
    cs["p"] = cs["pValueMantissa"].astype(float) * 10.0 ** cs["pValueExponent"].astype(float)
    cs = cs[(~cs["isTransQtl"].fillna(False)) & (cs["p"] < P_INSTRUMENT)].copy()
    cs["gene"] = cs["studyId"].map(_study_genes())
    cs = cs.dropna(subset=["gene", "beta", "zScore"])
    cs = cs[cs["zScore"] != 0]
    cs["se"] = (cs["beta"] / cs["zScore"]).abs()
    mhc = (cs["chromosome"] == MHC[0]) & cs["position"].between(MHC[1], MHC[2])
    cs = cs[~mhc]
    return cs


def _variant_keys(cs: pd.DataFrame) -> set[str]:
    keys = set(cs["variantId"])
    for s in cs["locus"]:
        keys.update(v for v, _, _ in json.loads(s))
    return keys


def outcome(keys: set[str]) -> pd.DataFrame:
    """UC GWAS rows for the instrument variants and their tags (chr_pos_ref_alt keys in either allele order)."""
    want = {}
    for k in keys:
        c, p, a, b = k.split("_")
        want[(c, p)] = True
    rows = []
    with gzip.open(GWAS, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            c, p = x[ix["hm_chrom"]], x[ix["hm_pos"]]
            if (c, p) in want:
                rows.append({"chrom": c, "pos": p, "effect": x[ix["hm_effect_allele"]], "other": x[ix["hm_other_allele"]],
                             "beta": x[ix["hm_beta"]], "se": x[ix["standard_error"]], "p": x[ix["p_value"]],
                             "rsid": x[ix["hm_rsid"]]})
    g = pd.DataFrame(rows)
    for c in ("beta", "se", "p"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    return g.dropna(subset=["beta", "se"])


def outcome_finngen(keys: set[str]) -> pd.DataFrame:
    """FinnGen rows for the instrument variants (same columns as outcome(); the effect allele is ALT)."""
    want = {tuple(k.split("_")[:2]) for k in keys}
    rows = []
    with gzip.open(FINNGEN, "rt") as f:
        f.readline()
        for line in f:
            x = line.split("\t", 10)
            if (x[0], x[1]) in want:
                rows.append({"chrom": x[0], "pos": x[1], "effect": x[3], "other": x[2], "beta": x[8], "se": x[9],
                             "p": x[6], "rsid": x[4]})
    g = pd.DataFrame(rows)
    for c in ("beta", "se", "p"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    return g.dropna(subset=["beta", "se"])


def outcome_ssf(path, keys: set[str]) -> pd.DataFrame:
    """GWAS-SSF harmonised file (GWAS Catalog): rows for the instrument positions, same columns as outcome()."""
    want = {tuple(k.split("_")[:2]) for k in keys}
    rows = []
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            if (x[ix["chromosome"]], x[ix["base_pair_location"]]) in want:
                rows.append({"chrom": x[ix["chromosome"]], "pos": x[ix["base_pair_location"]],
                             "effect": x[ix["effect_allele"]], "other": x[ix["other_allele"]], "beta": x[ix["beta"]],
                             "se": x[ix["standard_error"]], "p": x[ix["p_value"]], "rsid": x[ix.get("rsid", 0)]})
    g = pd.DataFrame(rows)
    for c in ("beta", "se", "p"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    return g.dropna(subset=["beta", "se"])


def attach_liu2023(r: pd.DataFrame, cs: pd.DataFrame) -> pd.DataFrame:
    """Repeat each selected estimate with Liu 2023 as outcome (replication column only; overlapping samples)."""
    if not LIU2023.exists():
        return r
    e = per_study(harmonise(cs, outcome_ssf(LIU2023, _variant_keys(cs)))).rename(
        columns={"beta": "beta_liu2023", "se": "se_liu2023", "p": "p_liu2023", "n_instruments": "n_instruments_liu2023"})
    r = r.drop(columns=[c for c in e.columns if c not in ("gene", "studyId") and c in r], errors="ignore")
    r = r.merge(e, on=["gene", "studyId"], how="left")
    has = r["beta_liu2023"].notna()
    r["replicates_liu2023"] = pd.Series(np.where(has, (np.sign(r["beta_liu2023"]) == np.sign(r["beta"]))
                                                 & (r["p_liu2023"] < 0.05), None), index=r.index, dtype=object)
    return r


def per_study(h: pd.DataFrame) -> pd.DataFrame:
    """Wald ratio or IVW per gene and study (no collapse to one study per gene)."""
    out = []
    for (gene, study), d in h.groupby(["gene", "studyId"]):
        ratio = d["beta_outcome"] / d["beta_exposure"]
        se = d["se_outcome"] / d["beta_exposure"].abs()
        w = 1 / se ** 2
        b = float((w * ratio).sum() / w.sum())
        s = float(np.sqrt(1 / w.sum()))
        out.append({"gene": gene, "studyId": study, "beta": b, "se": s, "p": float(2 * stats.norm.sf(abs(b / s))),
                    "n_instruments": len(d)})
    return pd.DataFrame(out)


def attach_finngen(r: pd.DataFrame, cs: pd.DataFrame, p_factor: pd.Series | None = None) -> pd.DataFrame:
    """Repeat each selected estimate (gene, study) with FinnGen as outcome and combine the two outcomes by inverse
    variance (fixed effect; the samples do not overlap). BH over genes for the combined p. p_factor (the cell layer's
    Bonferroni factor) multiplies the combined p before BH, as for the primary estimate."""
    from statsmodels.stats.multitest import multipletests
    if not FINNGEN.exists():
        return r
    g2 = outcome_finngen(_variant_keys(cs))
    e2 = per_study(harmonise(cs, g2)).rename(columns={"beta": "beta_finngen", "se": "se_finngen", "p": "p_finngen",
                                                       "n_instruments": "n_instruments_finngen"})
    r = r.drop(columns=[c for c in e2.columns if c not in ("gene", "studyId") and c in r], errors="ignore")
    r = r.merge(e2, on=["gene", "studyId"], how="left")
    w1, w2 = 1 / r["se"] ** 2, 1 / r["se_finngen"] ** 2
    has = r["beta_finngen"].notna()
    r["beta_combined"] = np.where(has, (w1 * r["beta"] + w2 * r["beta_finngen"]) / (w1 + w2), r["beta"])
    r["se_combined"] = np.where(has, np.sqrt(1 / (w1 + w2)), r["se"])
    pc = 2 * stats.norm.sf(np.abs(r["beta_combined"] / r["se_combined"]))
    if p_factor is not None:
        pc = np.clip(pc * p_factor.reindex(r.index).fillna(1).to_numpy(), 0, 1)
    r["p_combined"] = pc
    r["q_combined"] = multipletests(r["p_combined"], method="fdr_bh")[1]
    r["replicates_finngen"] = np.where(has, (np.sign(r["beta_finngen"]) == np.sign(r["beta"])) & (r["p_finngen"] < 0.05), None)
    r["replicates_finngen"] = r["replicates_finngen"].astype(object)
    return r


def _lookup(gw: dict, vid: str):
    """GWAS effect of the variant's ALT allele (Open Targets variant ids are chr_pos_ref_alt, beta for alt)."""
    c, p, ref, alt = vid.split("_")
    for r in gw.get((c, p), []):
        if r["effect"] == alt and r["other"] == ref:
            return r["beta"], r["se"], r["p"], r["rsid"]
        if r["effect"] == ref and r["other"] == alt:
            return -r["beta"], r["se"], r["p"], r["rsid"]
    return None


def harmonise(cs: pd.DataFrame, g: pd.DataFrame) -> pd.DataFrame:
    gw: dict = {}
    for r in g.to_dict("records"):
        gw.setdefault((r["chrom"], r["pos"]), []).append(r)
    rows = []
    for r in cs.itertuples():
        c, p, ref, alt = r.variantId.split("_")
        if len(ref) != 1 or len(alt) != 1 or _palindromic(ref, alt):
            continue
        used, hit, b_exp, pp = r.variantId, _lookup(gw, r.variantId), r.beta, None
        if hit is None:  # another variant of the same 95 % credible set, with its own pQTL effect (same allele order)
            for tag, t_pp, t_beta in sorted(json.loads(r.locus), key=lambda x: -x[1]):
                tc, tp, tref, talt = tag.split("_")
                if tag == r.variantId or len(tref) != 1 or len(talt) != 1 or _palindromic(tref, talt) or not t_beta:
                    continue
                hit = _lookup(gw, tag)
                if hit is not None:
                    used, b_exp, pp = tag, t_beta, t_pp
                    break
            if hit is None:
                continue
        b_out, se_out, p_out, rsid = hit
        rows.append({"gene": r.gene, "studyId": r.studyId, "studyLocusId": r.studyLocusId, "lead_variant": r.variantId,
                     "gwas_variant": used, "rsid": rsid, "substitute_pp": pp, "beta_exposure": b_exp, "se_exposure": r.se,
                     "p_exposure": r.p, "beta_outcome": b_out, "se_outcome": se_out, "p_outcome": p_out})
    h = pd.DataFrame(rows)
    # one instrument per credible set; independent credible sets of one protein are combined
    return h.drop_duplicates(["studyLocusId"]).drop_duplicates(["studyId", "gwas_variant"])


def estimate(h: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (gene, study), d in h.groupby(["gene", "studyId"]):
        ratio = d["beta_outcome"] / d["beta_exposure"]
        se = d["se_outcome"] / d["beta_exposure"].abs()
        w = 1 / se ** 2
        b = float((w * ratio).sum() / w.sum())
        s = float(np.sqrt(1 / w.sum()))
        het_p = float(stats.chi2.sf(float((w * (ratio - b) ** 2).sum()), len(d) - 1)) if len(d) > 1 else None
        out.append({"gene": gene, "studyId": study, "n_instruments": len(d), "method": "IVW" if len(d) > 1 else "Wald ratio",
                    "beta": b, "se": s, "p": float(2 * stats.norm.sf(abs(b / s))), "or_per_sd": float(np.exp(b)),
                    "heterogeneity_p": het_p, "lead_instrument": d.sort_values("p_exposure").iloc[0]["lead_variant"],
                    "min_p_exposure": float(d["p_exposure"].min()), "uses_tag_variant": bool(d["substitute_pp"].notna().any())})
    r = pd.DataFrame(out)
    # a protein measured by several Olink assays keeps the assay with the most instruments, then the strongest
    r = r.sort_values(["gene", "n_instruments", "min_p_exposure"], ascending=[True, False, True]).drop_duplicates("gene")
    from statsmodels.stats.multitest import multipletests
    r["q"] = multipletests(r["p"], method="fdr_bh")[1]
    return r.sort_values("p").reset_index(drop=True)


# ------------------------------------------------------------------ colocalisation check (Open Targets API)
COLOC_Q = """query($id:String!){ credibleSet(studyLocusId:$id){ colocalisation(studyTypes:[gwas], page:{index:0,size:500}){
  rows{ h4 clpp colocalisationMethod otherStudyLocus{ studyLocusId study{ id traitFromSource diseases{ id name } } } } } } }"""
IBD = re.compile(r"ulcerative colitis|inflammatory bowel|crohn", re.I)


COLOC_CACHE = RAW / "extra" / "ot_coloc_cache.json"   # API answers kept, so a rebuild does not depend on the network


def _coloc_query(sl: str) -> dict:
    import time
    cache = json.loads(COLOC_CACHE.read_text(encoding="utf-8")) if COLOC_CACHE.exists() else {}
    if sl in cache:
        return cache[sl]
    for attempt in range(5):
        try:
            r = requests.post(API, json={"query": COLOC_Q, "variables": {"id": sl}}, timeout=120)
            if r.status_code == 200 and r.text.strip():
                data = r.json()
                cache[sl] = data
                COLOC_CACHE.write_text(json.dumps(cache), encoding="utf-8")
                return data
        except (requests.RequestException, ValueError):
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Open Targets colocalisation query failed for {sl}")


def coloc(study_locus_ids: list[str]) -> dict[str, dict]:
    out = {}
    for sl in study_locus_ids:
        r = _coloc_query(sl)
        rows = (((r.get("data") or {}).get("credibleSet") or {}).get("colocalisation") or {}).get("rows") or []
        best = {"h4": None, "clpp": None, "trait": None, "study": None}
        for x in rows:
            st = (x.get("otherStudyLocus") or {}).get("study") or {}
            names = " ".join([st.get("traitFromSource") or ""] + [d["name"] for d in st.get("diseases") or []])
            if not IBD.search(names):
                continue
            h4, clpp = x.get("h4"), x.get("clpp")
            if (h4 or 0) > (best["h4"] or 0) or ((clpp or 0) > (best["clpp"] or 0)):
                best = {"h4": max(h4 or 0, best["h4"] or 0) if h4 is not None or best["h4"] is not None else None,
                        "clpp": max(clpp or 0, best["clpp"] or 0) if clpp is not None or best["clpp"] is not None else None,
                        "trait": st.get("traitFromSource"), "study": st.get("id")}
        out[sl] = best
    return out


def run(q_coloc: float = 0.05) -> dict:
    cs = instruments()
    g = outcome(_variant_keys(cs))
    h = harmonise(cs, g)
    r = estimate(h)
    sig = r[r["q"] < q_coloc]
    loci = h[h["gene"].isin(sig["gene"]) & h["studyId"].isin(sig["studyId"])]
    cl = coloc(sorted(loci["studyLocusId"].unique()))
    best = {}
    for row in loci.itertuples():
        c = cl.get(row.studyLocusId) or {}
        b = best.setdefault(row.gene, {"coloc_h4": None, "coloc_clpp": None, "coloc_trait": None})
        for k_src, k in (("h4", "coloc_h4"), ("clpp", "coloc_clpp")):
            if c.get(k_src) is not None and (b[k] is None or c[k_src] > b[k]):
                b[k] = round(float(c[k_src]), 4)
                b["coloc_trait"] = c.get("trait")
    r["coloc_h4"] = r["gene"].map(lambda x: best.get(x, {}).get("coloc_h4"))
    r["coloc_clpp"] = r["gene"].map(lambda x: best.get(x, {}).get("coloc_clpp"))
    r["coloc_trait"] = r["gene"].map(lambda x: best.get(x, {}).get("coloc_trait"))
    coloc_yes = (r["coloc_h4"].astype(float).fillna(0) >= 0.8) | (r["coloc_clpp"].astype(float).fillna(0) >= 0.01)
    r["colocalises"] = pd.Series(np.where(r["q"] < q_coloc, coloc_yes, None), index=r.index, dtype=object)
    r = attach_finngen(r, cs)
    r = attach_liu2023(r, cs)
    r.to_parquet(PROCESSED / "mr_pqtl.parquet", index=False)
    h.to_parquet(PROCESSED / "mr_instruments.parquet", index=False)
    return {"proteins_tested": len(r), "instruments": len(h), "significant_q05": int((r["q"] < 0.05).sum()),
            "colocalising": int(r["colocalises"].fillna(False).astype(bool).sum())}


def load(store=None) -> int:
    """One mendelian_randomization evidence record per tested protein (significant or not)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    r = pd.read_parquet(PROCESSED / "mr_pqtl.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:MR_UKBPPP_UC", "Dataset",
                       "Mendelian randomisation: UKB-PPP plasma cis-pQTL (Open Targets 26.09) on UC GWAS (de Lange 2017)")
        for x in r.itertuples():
            node = store.resolve(x.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": x.gene, "exposure": "plasma protein level (UKB-PPP, per SD)", "outcome": "ulcerative colitis",
                   "method": x.method, "n_instruments": int(x.n_instruments), "beta": round(float(x.beta), 4),
                   "se": round(float(x.se), 4), "odds_ratio_per_sd": round(float(x.or_per_sd), 3),
                   "p": float(f"{x.p:.3g}"), "q": float(f"{x.q:.3g}"), "significant": bool(x.q < 0.05),
                   "direction": "higher protein raises UC risk" if x.beta > 0 else "higher protein lowers UC risk",
                   "lead_instrument": x.lead_instrument, "uses_tag_variant": bool(x.uses_tag_variant),
                   "colocalisation_h4": None if pd.isna(x.coloc_h4) else float(x.coloc_h4),
                   "colocalisation_clpp": None if pd.isna(x.coloc_clpp) else float(x.coloc_clpp),
                   "colocalises_with_ibd_gwas": None if x.colocalises is None or pd.isna(x.colocalises) else bool(x.colocalises),
                   "study": x.studyId}
            if hasattr(x, "beta_finngen") and pd.notna(x.beta_finngen):
                rec.update({"finngen_odds_ratio_per_sd": round(float(np.exp(x.beta_finngen)), 3),
                            "finngen_p": float(f"{x.p_finngen:.3g}"),
                            "replicates_in_finngen": bool(x.replicates_finngen),
                            "combined_odds_ratio_per_sd": round(float(np.exp(x.beta_combined)), 3),
                            "combined_q": float(f"{x.q_combined:.3g}")})
            if hasattr(x, "beta_liu2023") and pd.notna(x.beta_liu2023):
                rec.update({"liu2023_odds_ratio_per_sd": round(float(np.exp(x.beta_liu2023)), 3),
                            "liu2023_p": float(f"{x.p_liu2023:.3g}"), "replicates_in_liu2023": bool(x.replicates_liu2023)})
            ev = store.add_evidence(evidence_id("mendelian_randomization", x.gene), "mendelian_randomization",
                                    "DATASET:MR_UKBPPP_UC", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "causal_estimate_for", "DATASET:MR_UKBPPP_UC", ev)
            n += 1
    return n


if __name__ == "__main__":
    import sys
    if "download" in sys.argv:
        print(download())
    else:
        print(json.dumps(run(), indent=1))
        print(pd.read_parquet(PROCESSED / "mr_pqtl.parquet").head(25).to_string())
