"""Genetics-first target list: UC risk signal -> gene -> cell -> direction (added 2026-11-12).

The three target lists start from what changes in inflamed biopsies, which mixes cause and consequence. This list starts
from inherited risk variants, which precede the disease.

Risk signals (GRCh38, MHC excluded), in this order of preference:
1. Open Targets 26.09 fine-mapped credible sets (SuSiE and PICS, several independent signals per locus allowed) of the
   UC GWAS: Liu 2023 (GCST90446794), de Lange 2017 (GCST004133), FinnGen R12, Million Veteran Program (GCST90480319)
   and UK Biobank WGS (GCST90473823); then of the IBD GWAS (Liu 2023 GCST90292538, de Lange 2017 GCST004131, FinnGen
   R12), kept only when the lead variant is at least nominally associated with UC (p < 1e-4 in Lyu 2025 or Liu 2023).
2. Loci of the UC meta-analysis of Lyu et al. 2025 (medRxiv 10.1101/2025.09.23.25336483, a preprint that is NOT peer
   reviewed; European meta-analysis of MVP, FinnGen R12, Pan-UK Biobank and IIBDGC) that no credible set above covers:
   p < 5e-8 in at least two cohorts, 500 kb clumping, approximate fine-mapping (Wakefield Bayes factors, one causal
   variant). They are flagged as "preprint only".

For every signal, AlphaGenome predicts the effect of each credible-set variant on gene expression (alphagenome_uc.py):
Atlas (colon tissue, blood T cells, blood monocytes) for every variant, AlphaGenome 2 single-cell expression in colon
and blood cell types for the most likely variants. A prediction counts per variant, with the variant's posterior
probability shown next to it. Two tiers, calibrated on fine-mapped eQTL (alphagenome_uc.validate_*):
A gene is assigned to a variant in a context when it is the top-scoring protein-coding gene of that variant there,
its |score| is at least twice the second gene's and above the size threshold of the eQTL check (see thresholds(): on
fine-mapped eQTL this picks the measured eQTL gene in 93 % of tissue and 100 % of cell-type cases, with the right
direction in 98 % and 94 %). Tiers by the causal probability of the variants that carry the effect: confident (at least
half of the signal's posterior probability) and likely (10 to 50 %). Size alone was tried first and rejected: it let
about five genes per variant through and agreed with Mendelian randomisation only at chance.
"""
from __future__ import annotations

import gzip
import io
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from . import alphagenome_uc as ag
from . import mr
from .paths import PROCESSED, RAW

EXTRACT = RAW / "extra" / "ot_uc_ibd_gwas_credible_sets.parquet"
PARTS = EXTRACT.parent / "_gwas_cs_parts"
LYU = RAW / "gwas" / "lyu2025_uc"
UC_STUDIES = {"GCST90446794": "Liu 2023 UC", "GCST004133": "de Lange 2017 UC", "FINNGEN_R12_K11_UC_STRICT2": "FinnGen R12 UC",
              "GCST90480319": "MVP 2024 UC", "GCST90473823": "UK Biobank WGS 2025 UC"}
IBD_STUDIES = {"GCST90292538": "Liu 2023 IBD", "GCST004131": "de Lange 2017 IBD", "FINNGEN_R12_K11_IBD_STRICT": "FinnGen R12 IBD"}
PRIORITY = list(UC_STUDIES) + list(IBD_STUDIES)
P_GW, P_UC_NOMINAL, MERGE_KB, MIN_COHORTS = 5e-8, 1e-4, 250, 2
MHC = ag.MHC


# ------------------------------------------------------------------ Open Targets GWAS credible sets
def _filter_part(name: str, ids) -> pd.DataFrame:
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    cache = PARTS / name
    if cache.exists():
        return pd.read_parquet(cache)
    for attempt in range(5):                     # the FTP mirror drops long transfers now and then
        try:
            r = requests.get(f"{mr.OT}/credible_set/{name}", timeout=900)
            r.raise_for_status()
            break
        except (requests.ConnectionError, requests.exceptions.ChunkedEncodingError, requests.Timeout):
            if attempt == 4:
                raise
    tab = pq.read_table(io.BytesIO(r.content), columns=mr.COLS)
    t = tab.filter(pc.is_in(tab["studyId"], value_set=ids)).to_pandas()
    del tab
    t["locus"] = t["locus"].map(lambda xs: json.dumps(
        [[x["variantId"], round(float(x["posteriorProbability"] or 0), 4), x["beta"]]
         for x in (xs if xs is not None else []) if x["is95CredibleSet"]]))
    PARTS.mkdir(parents=True, exist_ok=True)
    t.to_parquet(cache, index=False)
    return t


def download(workers: int = 6) -> dict:
    """Stream the 200 credible-set parts and keep the UC and IBD GWAS credible sets (parts cached, resumable)."""
    import shutil
    from concurrent.futures import ThreadPoolExecutor

    import pyarrow as pa

    from .download import sha256
    if EXTRACT.exists():
        return {"rows": len(pd.read_parquet(EXTRACT, columns=["studyLocusId"]))}
    ids = pa.array(sorted(PRIORITY))
    parts = mr._parts("credible_set")
    keep = []
    with ThreadPoolExecutor(workers) as ex:
        for i, t in enumerate(ex.map(lambda n: _filter_part(n, ids), parts)):
            keep.append(t)
            if i % 20 == 0:
                print(f"part {i + 1}/{len(parts)}: {sum(len(k) for k in keep)} credible sets", flush=True)
    out = pd.concat(keep, ignore_index=True).sort_values("studyLocusId").reset_index(drop=True)
    out.to_parquet(EXTRACT, index=False)
    shutil.rmtree(PARTS, ignore_errors=True)
    mpath = RAW / "manifest.json"
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    manifest["extra/ot_uc_ibd_gwas_credible_sets.parquet"] = {
        "url": f"{mr.OT}/credible_set/ ({len(parts)} parts, filtered to the UC and IBD GWAS studies "
               f"{', '.join(PRIORITY)})",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256(EXTRACT), "bytes": EXTRACT.stat().st_size}
    mpath.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    return {"rows": len(out)}


# ------------------------------------------------------------------ Lyu 2025 (preprint) summary statistics
def lyu_significant(p_max: float = P_UC_NOMINAL) -> pd.DataFrame:
    """SNVs of the Lyu 2025 European UC meta-analysis with p < p_max, present in at least MIN_COHORTS cohorts.
    Columns: chrom, pos, ref, alt (1000 Genomes reference / alternative), beta for alt, se, p, rsid, n_cohorts."""
    out = PROCESSED / "lyu2025_uc_p1e-4.parquet"
    if out.exists() and p_max == P_UC_NOMINAL:
        return pd.read_parquet(out)
    rows = []
    for c in range(1, 23):
        with gzip.open(LYU / f"EUR_UC_chr{c}.txt.gz", "rt") as f:
            h = f.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(h)}
            for line in f:
                x = line.rstrip("\n").split("\t")
                try:
                    p = float(x[ix["P.value"]])
                except ValueError:
                    continue
                if p >= p_max:
                    continue
                _, pos, ref, alt = x[ix["MarkerName"]].split(":")
                if len(ref) != 1 or len(alt) != 1:
                    continue
                # MarkerName lists the two alleles alphabetically; the REF and ALT columns (1000 Genomes) say which
                # one is the reference. Without them the order stays alphabetical and effects() checks it.
                if x[ix["REF"]] in (ref, alt) and x[ix["ALT"]] in (ref, alt) and x[ix["REF"]] != x[ix["ALT"]]:
                    ref, alt = x[ix["REF"]], x[ix["ALT"]]
                d = x[ix["Direction"]]
                n = sum(ch in "+-" for ch in d)
                if n < MIN_COHORTS:
                    continue
                a1 = x[ix["Allele1"]].upper()
                beta = float(x[ix["Effect"]])
                rows.append({"chrom": str(c), "pos": int(pos), "ref": ref, "alt": alt,
                             "beta": beta if a1 == alt else -beta, "se": float(x[ix["StdErr"]]), "p": p,
                             "rsid": x[ix["SNP"]] if x[ix["SNP"]] != "NA" else None, "n_cohorts": n,
                             "het_i2": float(x[ix["HetISq"]])})
    g = pd.DataFrame(rows)
    g = g[~((g.chrom == MHC[0]) & g.pos.between(MHC[1], MHC[2]))].reset_index(drop=True)
    g["variant"] = g.chrom + "_" + g.pos.astype(str) + "_" + g.ref + "_" + g.alt
    if p_max == P_UC_NOMINAL:
        g.to_parquet(out, index=False)
    return g


# ------------------------------------------------------------------ risk signals
AUTOSOMES = {str(i) for i in range(1, 23)}


def protein_coding() -> set[str]:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "locus_group"], dtype=str)
    return set(hg.loc[hg.locus_group == "protein-coding gene", "symbol"])


def _explode(cs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in cs.itertuples():
        for vid, pip, beta in json.loads(r.locus):
            c, pos, ref, alt = vid.split("_")
            if len(ref) != 1 or len(alt) != 1 or c not in AUTOSOMES:
                continue
            if c == MHC[0] and MHC[1] <= int(pos) <= MHC[2]:
                continue
            rows.append({"cs": r.studyLocusId, "study": r.studyId, "method": r.finemappingMethod, "lead": r.variantId,
                         "lead_p": r.p, "variant": vid, "chrom": c, "pos": int(pos), "pip": pip,
                         "beta": beta if beta is not None else np.nan, "lead_beta": r.beta})
    return pd.DataFrame(rows)


def signals() -> pd.DataFrame:
    """One row per credible-set variant of every non-redundant UC risk signal (see the module docstring)."""
    cs = pd.read_parquet(EXTRACT)
    cs["p"] = cs["pValueMantissa"].astype(float) * 10.0 ** cs["pValueExponent"].astype(float)
    cs = cs[cs.p < P_GW]
    v = _explode(cs)
    lyu = lyu_significant()
    lyu_ix = lyu.set_index("variant")
    # IBD credible sets need a nominal UC association of one of their variants
    v["uc_nominal"] = v.variant.isin(lyu_ix.index)
    ok_ibd = set(v[v.uc_nominal].cs)
    v = v[v.study.isin(UC_STUDIES) | v.cs.isin(ok_ibd)]
    rank = {s: i for i, s in enumerate(PRIORITY)}
    order = (v.drop_duplicates("cs").assign(r=lambda d: d.study.map(rank)).sort_values(["r", "lead_p"]).cs.tolist())
    by_cs = {k: g for k, g in v.groupby("cs")}
    kept, seen, leads = [], set(), []          # leads: (chrom, pos, study)
    for k in order:
        g = by_cs[k]
        c0, p0 = g.lead.iloc[0].split("_")[:2]
        if seen & set(g.variant):
            continue
        if any(c == c0 and abs(p - int(p0)) < MERGE_KB * 1000 and s != g.study.iloc[0] for c, p, s in leads):
            continue
        kept.append(k)
        seen |= set(g.variant)
        leads.append((c0, int(p0), g.study.iloc[0]))
    out = v[v.cs.isin(kept)].copy()
    out["source"] = out.study.map({**UC_STUDIES, **IBD_STUDIES})
    out["ibd_signal"] = out.study.isin(IBD_STUDIES)
    out["preprint_only"] = False
    # Lyu 2025 loci that no credible set covers: approximate fine-mapping
    gw = lyu[lyu.p < P_GW].sort_values("p")
    new_leads = []
    for r in gw.itertuples():
        if any(c == r.chrom and abs(p - r.pos) < ag.CLUMP_KB * 1000 for c, p, _ in leads):
            continue
        if any(c == r.chrom and abs(p - r.pos) < ag.CLUMP_KB * 1000 for c, p in new_leads):
            continue
        new_leads.append((r.chrom, r.pos))
    extra = []
    for c, p in new_leads:
        g = lyu[(lyu.chrom == c) & ((lyu.pos - p).abs() < ag.CLUMP_KB * 1000)].copy()
        z, var = g.beta / g.se, g.se ** 2
        r_ = ag.PRIOR_SD ** 2 / (ag.PRIOR_SD ** 2 + var)
        labf = 0.5 * (np.log(1 - r_) + r_ * z ** 2)
        g["pip"] = np.exp(labf - labf.max()) / np.exp(labf - labf.max()).sum()
        g = g.sort_values("pip", ascending=False)
        g = g[g.pip.cumsum().shift(fill_value=0) < ag.CS]
        lead = g.sort_values("p").iloc[0]
        extra.append(pd.DataFrame({"cs": f"LYU2025_{c}_{p}", "study": "LYU2025_EUR_UC", "method": "ABF (one causal variant)",
                                   "lead": lead.variant, "lead_p": lead.p, "variant": g.variant, "chrom": g.chrom,
                                   "pos": g.pos, "pip": g.pip, "beta": g.beta, "lead_beta": lead.beta,
                                   "uc_nominal": True, "source": "Lyu 2025 UC (preprint)", "ibd_signal": False,
                                   "preprint_only": True}))
    out = pd.concat([out] + extra, ignore_index=True)
    out["beta"] = out.beta.fillna(out.variant.map(lyu_ix.beta)).fillna(out.lead_beta)
    out["risk_sign"] = np.where(out.beta >= 0, 1, -1)          # beta is for the alternative allele
    out["rsid"] = out.variant.map(lyu_ix.rsid)
    out["lead_rsid"] = out.lead.map(lyu_ix.rsid)
    out["uc_p_lyu2025"] = out.variant.map(lyu_ix.p)
    sig = {k: f"S{i + 1:03d}" for i, k in enumerate(
        out.drop_duplicates("cs").assign(c=lambda d: d.chrom.astype(int)).sort_values(["c", "pos"]).cs)}
    out["signal"] = out.cs.map(sig)
    out.to_parquet(PROCESSED / "genetic_signals.parquet", index=False)
    return out


# ------------------------------------------------------------------ AlphaGenome effects per variant
AG2_EXTRA = {   # immune conditions added for the genetics-first list (single-cell expression tracks observed)
    "UBERON:0001155_CL:0000623": ("colon NK cell", "T/NK"),
    "UBERON:0001155_CL:0000990": ("colon dendritic cell", "myeloid"),
    "UBERON:0000178_CL:0002677": ("blood regulatory T cell", "T/NK"),
}
CONDITIONS = {**ag.AG2_CONDITIONS, **AG2_EXTRA}
ATLAS_PIP_MIN, ATLAS_PER_SIGNAL = 0.01, 50
AG2_PIP_MIN, AG2_PER_SIGNAL = 0.1, 3
ATLAS_LINEAGE = {"colon tissue": "tissue", "blood T cells": "T/NK", "blood monocytes": "myeloid"}


RATIO = 2.0   # the assigned gene must score at least twice the second gene of the same variant and context


def _top_gene(d: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Rows of d (one per gene) that are the top-scoring protein-coding gene of their variant and context, with the
    ratio of their |score| to the second gene's (inf when there is no second gene)."""
    d = d[d.gene.isin(protein_coding())].dropna(subset=["score"]).copy()
    d["ab"] = d.score.abs()
    d = d.sort_values("ab", ascending=False)
    d["_r"] = d.groupby(by).cumcount()
    second = d[d["_r"] == 1][by + ["ab"]].rename(columns={"ab": "second"})
    top = d[d["_r"] == 0].drop(columns="_r").merge(second, on=by, how="left")
    top["ratio"] = top.ab / top.second.where(top.second > 0)
    top["ratio"] = top.ratio.fillna(np.inf)
    return top


def thresholds() -> dict:
    """The rule that assigns a gene to a variant, and how often it is right on fine-mapped eQTL.

    A first version called an effect confident from its size alone (the size at which the DIRECTION of a known eQTL is
    right). That says nothing about WHICH gene: of the other genes near an eQTL variant, 2.4 % pass the cell-type
    threshold, about five per variant over 15 cell types, and the list it gave agreed with Mendelian randomisation only
    at chance (3 of 7). The rule now also requires the gene to be the top-scoring gene of that variant in that context,
    at least RATIO times the second one. On the same eQTL, that gene is the measured eQTL gene in the share reported
    here (`right_gene`), and its direction is right in `right_direction`."""
    a = pd.read_parquet(PROCESSED / "alphagenome_validation_atlas.parquet")
    b = pd.read_parquet(PROCESSED / "alphagenome_validation_ag2.parquet")
    b = b[b.kind == "matched"]
    th = {"atlas_confident": float(a["abs"].quantile(.75)), "ag2_confident": float(b.sc_expression.abs().quantile(.5)),
          "ratio": RATIO}
    v = ag._finemapped("gtex-v10_ge_colon").sample(frac=1, random_state=1).head(300)     # the variants of validate_atlas
    s = ag.atlas_scores(v.variant.tolist()).rename(columns={"colon tissue": "score"})
    top = _top_gene(s[["variant", "gene", "score"]], ["variant"]).merge(
        a[["variant", "gene", "agree"]].rename(columns={"gene": "true_gene"}), on="variant")
    parts = [ag.ag2_scores(r.variant, r.condition).rename(columns={"sc_expression": "score"}) for r in b.itertuples()]
    s2 = pd.concat(parts, ignore_index=True)
    top2 = _top_gene(s2[["variant", "condition", "gene", "score"]], ["variant", "condition"]).merge(
        b.assign(agree=np.sign(b.sc_expression) == np.sign(b.beta))[["variant", "condition", "gene", "agree"]]
        .rename(columns={"gene": "true_gene"}), on=["variant", "condition"])

    def acc(x, thr):
        ok = x[(x.ab >= thr) & (x.ratio >= RATIO)]
        hit = ok[ok.gene == ok.true_gene]
        return {"pairs": int(len(ok)), "right_gene": round(float((ok.gene == ok.true_gene).mean()), 3),
                "right_direction": round(float(hit.agree.mean()), 3), "of_checked_eqtl": int(len(x))}

    th["check"] = {"atlas": acc(top, th["atlas_confident"]), "ag2": acc(top2, th["ag2_confident"]),
                   "top_gene_alone": {"atlas": round(float((top.gene == top.true_gene).mean()), 3),
                                      "ag2": round(float((top2.gene == top2.true_gene).mean()), 3)}}
    return th


def ag2_plan(s: pd.DataFrame) -> pd.DataFrame:
    return s[s.pip >= AG2_PIP_MIN].sort_values("pip", ascending=False).groupby("signal").head(AG2_PER_SIGNAL)


FLIPS = ag.CACHE / "flipped_variants.json"


def _swap(v: str) -> str:
    c, p, ref, alt = v.split("_")
    return f"{c}_{p}_{alt}_{ref}"


def _atlas_checked(variants: list[str]) -> tuple[pd.DataFrame, set[str]]:
    """Atlas scores; a variant whose first allele is not the genome reference (a few GWAS files list alleles the other
    way round) is queried with its alleles swapped and its score sign reversed. Returns the scores keyed on the
    original variant ids and the set of swapped variants (saved, so a rebuild sends the same queries)."""
    import re
    flipped = set(json.loads(FLIPS.read_text(encoding="utf-8"))) if FLIPS.exists() else set()
    while True:
        q = sorted(_swap(v) if v in flipped else v for v in variants)
        try:
            at = ag.atlas_scores(q)
            break
        except ValueError as exc:
            m = re.search(r"chr(\w+):(\d+):([ACGT])>([ACGT]) reference base does not match", str(exc))
            if not m:
                raise
            bad = "_".join(m.groups())
            orig = bad if bad in variants else _swap(bad)
            if orig in flipped:                       # neither allele is the reference: leave the variant out
                variants = [v for v in variants if v != orig]
                flipped.discard(orig)
            else:
                flipped.add(orig)
            FLIPS.write_text(json.dumps(sorted(flipped)), encoding="utf-8")
    back = {_swap(v): v for v in flipped}
    sw = at.variant.isin(back)
    for col in ag.ATLAS_GROUPS:
        at.loc[sw, col] = -at.loc[sw, col]
    at["variant"] = at.variant.map(lambda v: back.get(v, v))
    return at, flipped


def _ag2_checked(variant: str, condition: str, flipped: set[str]) -> pd.DataFrame:
    d = ag.ag2_scores(_swap(variant) if variant in flipped else variant, condition).copy()
    if variant in flipped:
        d["sc_expression"], d["rna_seq"], d["variant"] = -d.sc_expression, -d.rna_seq, variant
    return d


def effects(workers: int = 4) -> pd.DataFrame:
    """Risk-signed predicted effect of every scored variant on every gene and context that reaches a tier."""
    from concurrent.futures import ThreadPoolExecutor
    s = pd.read_parquet(PROCESSED / "genetic_signals.parquet")
    th = thresholds()
    a_sel = (s[s.pip >= ATLAS_PIP_MIN].sort_values(["signal", "pip"], ascending=[True, False])
             .groupby("signal").head(ATLAS_PER_SIGNAL))
    plan = ag2_plan(s)
    at, flipped = _atlas_checked(sorted(set(a_sel.variant) | set(plan.variant)))
    at = at.melt(id_vars=["variant", "gene"], value_vars=list(ag.ATLAS_GROUPS), var_name="context", value_name="score")
    at = _top_gene(at, ["variant", "context"])
    at = at[(at.ab >= th["atlas_confident"]) & (at.ratio >= RATIO)].copy()
    at["source"], at["lineage"], at["tier"] = "Atlas", at.context.map(ATLAS_LINEAGE), "confident"
    scored = set(at.variant)
    jobs = [(v, c) for v in sorted(set(plan.variant) & scored) for c in CONDITIONS]
    if jobs:
        _ag2_checked(*jobs[0], flipped)                           # create the client once, before the threads
    with ThreadPoolExecutor(workers) as ex:
        parts = list(ex.map(lambda j: _ag2_checked(*j, flipped), jobs))
    a2 = pd.concat(parts, ignore_index=True).rename(columns={"sc_expression": "score"})
    a2 = _top_gene(a2, ["variant", "condition"])
    a2 = a2[(a2.ab >= th["ag2_confident"]) & (a2.ratio >= RATIO)].copy()
    a2["context"] = a2.condition.map(lambda c: CONDITIONS[c][0])
    a2["lineage"] = a2.condition.map(lambda c: CONDITIONS[c][1])
    a2["source"], a2["tier"] = "AlphaGenome 2", "confident"
    cols = ["variant", "gene", "context", "score", "ratio", "source", "lineage", "tier"]
    e = pd.concat([at[cols], a2[cols]], ignore_index=True)
    e = e.merge(s[["signal", "variant", "pip", "risk_sign"]].drop_duplicates(), on="variant")
    e["effect"] = e.score * e.risk_sign
    e.to_parquet(PROCESSED / "genetic_effects.parquet", index=False)
    return e


# ------------------------------------------------------------------ the list
TIER_POINTS = {"confident": 3, "likely": 2}
PIP_CONFIDENT, PIP_LIKELY, GENES_PER_SIGNAL = 0.5, 0.1, 3


def _support(e: pd.DataFrame, tier: str) -> pd.DataFrame:
    """Per signal, gene and context: causal probability carried by the variants that have an effect of this tier in
    the majority direction, the PIP-weighted effect of those variants and the most likely one of them."""
    x = e[e.tier == tier].copy()
    if x.empty:
        return pd.DataFrame()
    x["sign"] = np.sign(x.effect)
    g = x.groupby(["signal", "gene", "context", "source", "lineage", "sign"], as_index=False).agg(
        pip_support=("pip", "sum"), effect=("effect", lambda v: float(np.average(v, weights=x.loc[v.index, "pip"]))))
    g = g.sort_values("pip_support", ascending=False).drop_duplicates(["signal", "gene", "context"])
    top = x.sort_values("pip", ascending=False).drop_duplicates(["signal", "gene", "context", "sign"])
    g = g.merge(top[["signal", "gene", "context", "sign", "variant", "pip"]].rename(
        columns={"variant": "top_variant", "pip": "top_pip"}), on=["signal", "gene", "context", "sign"])
    g["pip_support"] = g.pip_support.clip(upper=1.0)
    return g


def targets() -> dict:
    from .validation_shortlist import _lineage_share, known_genes, repurposing_genes
    e = pd.read_parquet(PROCESSED / "genetic_effects.parquet")
    s = pd.read_parquet(PROCESSED / "genetic_signals.parquet")
    c = _support(e, "confident")
    c["tier"] = np.where(c.pip_support >= PIP_CONFIDENT, "confident", np.where(c.pip_support >= PIP_LIKELY, "likely", ""))
    c = c[c.tier != ""]
    c["n_contexts"] = c.groupby(["signal", "gene"]).context.transform("nunique")
    c["tp"] = c.tier.map(TIER_POINTS)
    c["strength"] = c.pip_support * c.effect.abs()
    c = c.sort_values(["tp", "strength"], ascending=False)
    best = c.drop_duplicates(["signal", "gene"]).copy()
    ctx = (c.groupby(["signal", "gene"]).apply(
        lambda d: "; ".join(f"{r.context} {r.effect:+.3f}" for r in d.head(4).itertuples()), include_groups=False)
        .rename("contexts").reset_index())
    best = best.merge(ctx, on=["signal", "gene"], how="left")
    best["contexts"] = best.contexts.fillna(best.context + " " + best.effect.map("{:+.3f}".format))
    best = best.sort_values(["tp", "strength"], ascending=False).groupby("signal").head(GENES_PER_SIGNAL)
    n_sig = best.groupby("gene").signal.nunique()
    best = best.sort_values(["tp", "strength"], ascending=False).drop_duplicates("gene")
    best["n_signals"] = best.gene.map(n_sig)
    info = s.sort_values("pip", ascending=False).drop_duplicates("signal").set_index("signal")
    for col in ("source", "method", "lead", "lead_rsid", "lead_p", "ibd_signal", "preprint_only", "chrom"):
        best[col if col != "source" else "gwas"] = best.signal.map(info[col])
    best["lead_pos"] = best.lead.str.split("_").str[1].astype(int)
    best["cs_size"] = best.signal.map(s.groupby("signal").size())
    rs = s.dropna(subset=["rsid"]).drop_duplicates("variant").set_index("variant").rsid
    best["top_rsid"] = best.top_variant.map(rs).fillna(best.top_variant)
    best["lead_rsid"] = best.lead_rsid.fillna(best.lead)
    best["risk_allele"] = np.where(best.effect > 0, "raises", "lowers")
    best["hypothesis"] = np.where(best.effect > 0, "inhibit", "activate or restore")
    # ---- what the rest of the atlas says about each gene
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    best["meta_log2fc"] = best.gene.map(meta.meta_log2fc)
    best["meta_padj"] = best.gene.map(meta.meta_padj)
    changed = (best.meta_padj < 0.05) & (best.meta_log2fc.abs() >= 0.5)
    best["tissue_change"] = np.select([changed & (np.sign(best.meta_log2fc) == np.sign(best.effect)), changed],
                                      ["same direction", "opposite direction"], "no change")
    lin = _lineage_share()
    best["atlas_lineage"] = best.gene.map(lambda g: lin.get(g, (None, None))[0])
    best["atlas_share"] = best.gene.map(lambda g: lin.get(g, (None, None))[1])
    best["cell_match"] = (best.lineage == best.atlas_lineage) & best.atlas_lineage.notna()
    mrs = []
    for f in ("mr_pqtl", "mr_colon"):
        m = pd.read_parquet(PROCESSED / f"{f}.parquet")
        m = m[m.q < 0.05].sort_values("q").drop_duplicates("gene")
        mrs.append(m.assign(layer="plasma protein" if f == "mr_pqtl" else m.get("label", "gut expression"))[
            ["gene", "beta", "or_per_sd", "q", "layer"]])
    mrt = pd.concat(mrs).sort_values("q").drop_duplicates("gene").set_index("gene")
    best["mr_or"], best["mr_q"], best["mr_layer"] = best.gene.map(mrt.or_per_sd), best.gene.map(mrt.q), best.gene.map(mrt.layer)
    best["mr"] = np.select([best.mr_q.notna() & (np.sign(best.gene.map(mrt.beta)) == np.sign(best.effect)), best.mr_q.notna()],
                           ["agrees", "disagrees"], "not significant or not testable")
    di = pd.read_parquet(PROCESSED / "drug_indications.parquet").set_index("gene")
    best["drug"] = best.gene.map(di.category).fillna("no clinical drug")
    best["drug_names"] = best.gene.map(di.uc_drugs.where(di.uc_drugs.astype(str).str.len() > 0, di.other_drugs))
    lists = {"well known": set(known_genes()), "repurposing": set(repurposing_genes())}
    sl = PROCESSED / "validation_shortlist.parquet"
    if sl.exists():
        lists["first-in-class"] = set(pd.read_parquet(sl, columns=["gene"]).gene)
    try:
        from .evaluate import clinical_targets
        lists["UC clinical target"] = set(clinical_targets().gene)
    except (ImportError, FileNotFoundError):
        pass
    best["in_lists"] = best.gene.map(lambda g: ", ".join(k for k, v in lists.items() if g in v))
    best["specific"] = ~best.ibd_signal & ~best.preprint_only
    best["points"] = (best.tp + 2 * (best.mr == "agrees") - 2 * (best.mr == "disagrees")
                      + (best.tissue_change == "same direction") + best.cell_match.astype(int) + best.specific.astype(int))
    best["group"] = np.select([(best.tier == "confident") & (best.points >= 5), best.tier == "confident"], ["A", "B"], "C")
    best = best.sort_values(["group", "points", "strength"], ascending=[True, False, False]).reset_index(drop=True)
    best["rank"] = np.arange(1, len(best) + 1)
    best.drop(columns=["tp", "sign"]).to_parquet(PROCESSED / "genetic_targets.parquet", index=False)
    th = thresholds()
    plan = ag2_plan(s)
    summary = {"signals": int(s.signal.nunique()), "variants": int(len(s)),
               "signals_by_gwas": s.drop_duplicates("signal").source.value_counts().to_dict(),
               "signals_ibd_only": int(s.drop_duplicates("signal").ibd_signal.sum()),
               "signals_preprint_only": int(s.drop_duplicates("signal").preprint_only.sum()),
               "signals_with_gene": int(best.signal.nunique()),
               "signals_with_confident_gene": int(best[best.tier == "confident"].signal.nunique()),
               "genes": int(len(best)), "by_tier": best.tier.value_counts().to_dict(),
               "by_group": best.group.value_counts().to_dict(),
               "in_other_lists": best[best.in_lists != ""].gene.tolist(),
               "mr_agrees": best[best.mr == "agrees"].gene.tolist(), "mr_disagrees": best[best.mr == "disagrees"].gene.tolist(),
               "ag2_variants": int(plan.variant.nunique()), "ag2_predictions": int(plan.variant.nunique() * len(CONDITIONS)),
               "conditions": [v[0] for v in CONDITIONS.values()], "thresholds": th}
    (PROCESSED / "genetic_targets_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return summary


def run() -> dict:
    signals()
    effects()
    out = targets()
    ag.write_manifest()
    return out


def load(store=None) -> dict:
    """One genetic_target record per gene of the list (new evidence type)."""
    from datetime import date

    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "genetic_targets.parquet")
    n = 0
    with store.batch():
        if store.get_node("DATASET:GENETIC_TARGETS") is None:
            store.add_node("DATASET:GENETIC_TARGETS", "Dataset", "UC risk signals (Open Targets credible sets; Lyu 2025 "
                           "preprint) with the gene, cell and direction predicted by AlphaGenome (non-commercial use)")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "signal": r.signal, "gwas": r.gwas, "fine_mapping": r.method,
                   "lead_variant": r.lead_rsid, "lead_p": float(f"{r.lead_p:.3g}"),
                   "ibd_signal_with_nominal_uc_association": bool(r.ibd_signal), "preprint_only": bool(r.preprint_only),
                   "tier": r.tier, "causal_probability_of_variants_with_this_effect": round(float(r.pip_support), 3),
                   "most_likely_variant": r.top_rsid, "most_likely_variant_pip": round(float(r.top_pip), 3),
                   "context": r.context, "model": r.source, "effect_of_risk_allele": round(float(r.effect), 4),
                   "risk_allele": f"{r.risk_allele} expression", "therapeutic_hypothesis": r.hypothesis,
                   "all_contexts": r.contexts, "group": r.group, "points": int(r.points),
                   "tissue_change_vs_prediction": r.tissue_change, "mendelian_randomisation": r.mr,
                   "interpretation": "predicted from DNA sequence, not measured; tiers calibrated on fine-mapped eQTL; "
                                     "non-commercial research use"}
            ev = store.add_evidence(evidence_id("genetic_target", r.gene), "genetic_target", "DATASET:GENETIC_TARGETS",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "genetic_target_candidate_in", "DATASET:GENETIC_TARGETS", ev)
            n += 1
    return {"genetic_target_records": n}
