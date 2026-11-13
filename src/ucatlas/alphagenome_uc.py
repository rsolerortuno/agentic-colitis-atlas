"""AlphaGenome 2: how UC risk variants change gene expression in colon tissue and colon cell types (added 2026-11-12).

AlphaGenome (Google DeepMind; Avsec et al. 2026, Nature; the API is free for non-commercial use, and its outputs are
for non-commercial research only, not for clinical decisions) predicts, from DNA sequence, the effect of a variant on
RNA and chromatin tracks. Two parts of the API are used:
- AlphaGenome Atlas: pre-computed scores for every single-nucleotide variant, RNA-seq effect per gene in tissue tracks
  (colon tracks: colonic mucosa, transverse, sigmoid, descending and left colon) and blood immune cells;
- AlphaGenome 2 live predictions for colon cell types (composite conditions such as "colon fibroblast"), where the
  Atlas has no tracks.

Variants: the UC GWAS of Liu et al. 2023 (GRCh38), single-nucleotide variants with p < 5e-8 outside the MHC, grouped
into loci by 500 kb clumping (142 loci). The lead variant is rarely the causal one, so each locus is fine-mapped
approximately: Wakefield approximate Bayes factors (prior effect SD 0.2 on the log-odds scale, one causal variant per
locus), posterior inclusion probabilities (PIP) and the 95 % credible set. Per locus and gene, the predicted effect of
the risk allele is the PIP-weighted mean of the predicted effects of the credible-set variants, each signed so that it
refers to the allele that raises UC risk.

Validation before use (validate_atlas, validate_ag2): fine-mapped single-variant eQTL (PIP >= 0.9) are the truth.
GTEx colon eQTL against Atlas colon tracks: sign agreement 64 % overall, 82 % for the top half of |score|, 94 % for
the top quarter, 97 % for the top tenth (284 eQTL-gene pairs). Predictions are therefore used only above the top-quarter |score|
of that validation (CONFIDENT); weaker ones are reported as uncertain. Cell-type predictions are checked the same way
against IBDverse cell-type eQTL (gut biopsies of IBD patients): 100 % agreement in the top quarter in the matched
colon cell type (Spearman 0.60) against 75 % in a mismatched one (Spearman 0.26), so the cell type matters.

Result (2026-11-12): 7 of 142 loci have a confident effect, on 17 genes; among the atlas's genes of interest only
CCL20 (risk allele raises it in blood T cells and monocytes) and IL10 (risk allele lowers it in T cells, as the
protective MR estimate expects). The candidate genes were chosen for their change in inflamed biopsies and rarely sit
at UC risk loci. No gene of interest reaches the AlphaGenome 2 threshold in a colon cell type.

Correction (2026-11-13): that result called a gene from the size of its predicted change alone, which tests the
direction of a known eQTL and not which gene a variant acts on. It is withdrawn; genetic_targets.py assigns genes with
a rule checked for that question, and under it CCL20 and IL10 are not assigned. This module keeps the API access, the
cache and the validation.

The API key is read from the ALPHAGENOME_API_KEY environment variable, never from a file in the repository. Every API
answer is cached (data/raw/alphagenome/, with a sha256 manifest) so that run-all reproduces without the key.
Outputs: data/processed/alphagenome_uc_loci.parquet (credible sets), alphagenome_uc.parquet (locus x gene x track
group), alphagenome_validation.parquet, alphagenome_summary.json.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

GWAS = RAW / "gwas" / "GCST90446794_uc_liu2023.h.tsv.gz"
CACHE = RAW / "alphagenome"
P_GW, CLUMP_KB, PRIOR_SD, CS = 5e-8, 500, 0.2, 0.95
MHC = ("6", 25_000_000, 34_000_000)
# Atlas RNA-seq track groups (biosample names of the Atlas RNA_SEQ tracks)
ATLAS_GROUPS = {"colon tissue": r"colonic mucosa|transverse colon|sigmoid colon|descending colon|left colon",
                "blood T cells": r"T cell", "blood monocytes": r"CD14-positive monocyte"}
# AlphaGenome 2 colon cell-type conditions (UBERON colon + Cell Ontology) and the lineage of this atlas they stand for
AG2_CONDITIONS = {
    "UBERON:0001155_CL:0000584": ("colon enterocyte", "epithelial"),
    "UBERON:0001155_CL:0000034": ("colon stem cell", "epithelial"),
    "UBERON:0001155_CL:0000160": ("colon goblet cell", "epithelial"),
    "UBERON:0001155_CL:0000057": ("colon fibroblast", "stromal"),
    "UBERON:0001155_CL:0000115": ("colon endothelial cell", "endothelial"),
    "UBERON:0001155_CL:0000235": ("colon macrophage", "myeloid"),
    "UBERON:0001155_CL:0000576": ("colon monocyte", "myeloid"),
    "UBERON:0001155_CL:0000624": ("colon CD4 T cell", "T/NK"),
    "UBERON:0001155_CL:0000625": ("colon CD8 T cell", "T/NK"),
    "UBERON:0001155_CL:0000786": ("colon plasma cell", "plasma"),
    "UBERON:0001155_CL:0000236": ("colon B cell", "B"),
    "UBERON:0000178_CL:0000775": ("blood neutrophil", "neutrophil"),
}


# ------------------------------------------------------------------ variants
def uc_loci() -> pd.DataFrame:
    """Genome-wide significant SNVs of Liu 2023, 500 kb loci, approximate fine-mapping (Wakefield ABF)."""
    rows = []
    with gzip.open(GWAS, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            try:
                p = float(x[ix["p_value"]])
            except ValueError:
                continue
            if p >= P_GW:
                continue
            ea, oa = x[ix["effect_allele"]], x[ix["other_allele"]]
            if len(ea) != 1 or len(oa) != 1:
                continue
            c, pos = x[ix["chromosome"]], int(x[ix["base_pair_location"]])
            if c == MHC[0] and MHC[1] <= pos <= MHC[2]:
                continue
            ref, alt = x[ix["variant_id"]].split("_")[2:4]
            rows.append({"chrom": c, "pos": pos, "ref": ref, "alt": alt, "effect": ea, "beta": float(x[ix["beta"]]),
                         "se": float(x[ix["standard_error"]]), "p": p, "rsid": x[ix["rsid"]]})
    g = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    locus = np.full(len(g), -1)
    leads = []
    for i, r in enumerate(g.itertuples()):
        hit = next((k for k, (c, p0) in enumerate(leads) if c == r.chrom and abs(r.pos - p0) < CLUMP_KB * 1000), None)
        if hit is None:
            leads.append((r.chrom, r.pos))
            hit = len(leads) - 1
        locus[i] = hit
    g["locus"] = locus
    # Wakefield approximate Bayes factor per variant, PIP within locus (one causal variant)
    z = g.beta / g.se
    v = g.se ** 2
    r_ = PRIOR_SD ** 2 / (PRIOR_SD ** 2 + v)
    g["log_abf"] = 0.5 * (np.log(1 - r_) + r_ * z ** 2)
    g["pip"] = g.groupby("locus").log_abf.transform(lambda s: np.exp(s - s.max()) / np.exp(s - s.max()).sum())
    g = g.sort_values(["locus", "pip"], ascending=[True, False])
    g["in_cs"] = g.groupby("locus").pip.transform(lambda s: s.cumsum().shift(fill_value=0) < CS)
    # sign that turns an alt-vs-ref effect into the effect of the risk-raising allele
    risk_is_alt = ((g.effect == g.alt) & (g.beta > 0)) | ((g.effect != g.alt) & (g.beta < 0))
    g["risk_sign"] = np.where(risk_is_alt, 1, -1)
    g["lead"] = g.groupby("locus").p.transform("min") == g.p
    g["variant"] = g.chrom + "_" + g.pos.astype(str) + "_" + g.ref + "_" + g.alt
    g.to_parquet(PROCESSED / "alphagenome_uc_loci.parquet", index=False)
    return g


# ------------------------------------------------------------------ API access with cache
def _key() -> str | None:
    return os.environ.get("ALPHAGENOME_API_KEY")


def _cache_path(kind: str, payload: str):
    return CACHE / kind / (hashlib.sha1(payload.encode()).hexdigest()[:20] + ".parquet")


def _variant(v: str):
    from alphagenome.data import genome
    c, p, ref, alt = v.split("_")
    return genome.Variant(chromosome=f"chr{c}", position=int(p), reference_bases=ref, alternate_bases=alt)


def atlas_scores(variants: list[str], chunk: int = 200) -> pd.DataFrame:
    """Atlas RNA-seq scores (alt vs ref) per variant and gene, averaged over each track group (cached per chunk)."""
    out = []
    client = meta = None
    for i in range(0, len(variants), chunk):
        part = variants[i:i + chunk]
        path = _cache_path("atlas", json.dumps(part))
        if not path.exists():
            if _key() is None:
                raise RuntimeError("AlphaGenome Atlas answer not cached and ALPHAGENOME_API_KEY not set")
            from alphagenome.atlas import atlas
            if client is None:
                client = atlas.create(_key())
                meta = client.scorer_metadata()["RNA_SEQ"].track_metadata
            ad = client.query_variants([_variant(v) for v in part], requested_scorers=["RNA_SEQ"],
                                       progress_bar=False)["RNA_SEQ"]
            obs = ad.obs.reset_index(drop=True)
            X = np.asarray(ad.X, dtype=float)
            rows = pd.DataFrame({"variant": obs.variant.astype(str).str.replace("chr", "", regex=False)
                                 .str.replace(":", "_", regex=False).str.replace(">", "_", regex=False),
                                 "gene": obs.gene_name.astype(str)})
            for name, pat in ATLAS_GROUPS.items():
                idx = np.where(meta.biosample_name.str.contains(pat, case=False, regex=True).values)[0]
                rows[name] = X[:, idx].mean(axis=1) if len(idx) else np.nan
            path.parent.mkdir(parents=True, exist_ok=True)
            rows.to_parquet(path, index=False)
        out.append(pd.read_parquet(path))
    return pd.concat(out, ignore_index=True)


def ag2_scores(variant: str, condition: str) -> pd.DataFrame:
    """AlphaGenome 2 live prediction for one variant in one condition: per gene, the RNA-seq LFC and the single-cell
    expression score, mean over tracks (cached)."""
    path = _cache_path("ag2", variant + "|" + condition)
    if not path.exists():
        if _key() is None:
            raise RuntimeError("AlphaGenome 2 answer not cached and ALPHAGENOME_API_KEY not set")
        from alphagenome.models.v2 import dna_client as dc, variant_scorers as vs
        global _AG2
        if _AG2 is None:
            _AG2 = dc.create(_key())
        v = _variant(variant)
        res = _AG2.score_variant(interval=v.reference_interval.resize(dc.SEQUENCE_LENGTH_1MB), variant=v,
                                 condition=condition,
                                 variant_scorers=[vs.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
                                                  vs.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION"]])
        cols = {}
        for ad, name in zip(res, ("rna_seq", "sc_expression")):
            cols[name] = pd.Series(np.asarray(ad.X, dtype=float).mean(axis=1),
                                   index=ad.obs.gene_name.astype(str).values).groupby(level=0).mean()   # repeated names
        d = pd.DataFrame(cols).reset_index().rename(columns={"index": "gene"})
        d = d.groupby("gene", as_index=False).mean()
        d["variant"], d["condition"] = variant, condition
        path.parent.mkdir(parents=True, exist_ok=True)
        d.to_parquet(path, index=False)
    return pd.read_parquet(path)


_AG2 = None


# ------------------------------------------------------------------ validation against fine-mapped eQTL
IBDVERSE_TO_AG2 = {   # IBDverse cell-type eQTL study prefix -> AlphaGenome 2 condition (matched) and a mismatched one
    "Colonocyte_RBFOX1+": "UBERON:0001155_CL:0000584", "Colonocyte_CAECAM7+_KRT20+_IFI27+": "UBERON:0001155_CL:0000584",
    "Colonocyte_stem_cell_LGR5+": "UBERON:0001155_CL:0000034", "Colon_Goblet_MUC2+_SPDEF+": "UBERON:0001155_CL:0000160",
    "PlasmaB_CD38+++": "UBERON:0001155_CL:0000786", "CD4+_memory": "UBERON:0001155_CL:0000624",
    "CD8+_TRM_TGCR2+": "UBERON:0001155_CL:0000625", "Macrophages_CD163++": "UBERON:0001155_CL:0000235",
    "Classical_monocytes": "UBERON:0001155_CL:0000576", "Naive_B_cells": "UBERON:0001155_CL:0000236"}
EPITHELIAL = {"UBERON:0001155_CL:0000584", "UBERON:0001155_CL:0000034", "UBERON:0001155_CL:0000160"}


def _finemapped(prefix: str, pip: float = 0.9) -> pd.DataFrame:
    """Single-variant eQTL credible sets (top PIP >= pip, SNV) from the Open Targets gut eQTL table."""
    d = pd.read_parquet(RAW / "extra" / "ot_gut_eqtl_credible_sets.parquet")
    d = d[d.studyId.str.startswith(prefix)]
    rows = []
    for r in d.itertuples():
        top = max(json.loads(r.locus), key=lambda x: x[1])
        c, p, ref, alt = top[0].split("_")
        if top[1] >= pip and len(ref) == 1 and len(alt) == 1:
            rows.append({"variant": top[0], "pip": top[1], "beta": top[2], "study": r.studyId,
                         "gene_id": r.studyId.rsplit("_", 1)[1].upper()})
    return pd.DataFrame(rows).drop_duplicates(["variant", "gene_id"])


def _gene_names() -> dict:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    return dict(zip(hg.ensembl_gene_id, hg.symbol))


def validate_atlas(n: int = 300) -> dict:
    """GTEx colon single-variant eQTL vs Atlas colon tracks: sign agreement by |score| quantile and the threshold."""
    v = _finemapped("gtex-v10_ge_colon").sample(frac=1, random_state=1).head(n)
    names = _gene_names()
    v["gene"] = v.gene_id.map(names)
    s = atlas_scores(v.variant.tolist())
    m = v.merge(s[["variant", "gene", "colon tissue"]], on=["variant", "gene"]).rename(columns={"colon tissue": "score"})
    m["agree"], m["abs"] = np.sign(m.score) == np.sign(m.beta), m.score.abs()
    out = {"pairs": int(len(m)), "agreement_all": round(float(m.agree.mean()), 3),
           "threshold_top_quarter": float(m["abs"].quantile(0.75)), "spearman": round(float(m.score.corr(m.beta, method="spearman")), 3)}
    for q in (0.5, 0.75, 0.9):
        out[f"agreement_top_{int(round((1 - q) * 100))}pct"] = round(float(m[m["abs"] >= m["abs"].quantile(q)].agree.mean()), 3)
    m.assign(check="atlas_gtex_colon").to_parquet(PROCESSED / "alphagenome_validation_atlas.parquet", index=False)
    return out


def validate_ag2(per_type: int = 25) -> dict:
    """IBDverse cell-type single-variant eQTL vs AlphaGenome 2 in the matched colon cell type and in a mismatched one."""
    names = _gene_names()
    rows = []
    for prefix, cond in IBDVERSE_TO_AG2.items():
        v = _finemapped(f"OTAR2057_IBDverse_ge_{prefix}_")
        if v.empty:
            continue
        v = v.sample(frac=1, random_state=1).head(per_type)
        other = "UBERON:0001155_CL:0000624" if cond in EPITHELIAL else "UBERON:0001155_CL:0000584"
        for r in v.itertuples():
            gene = names.get(r.gene_id)
            for kind, c in (("matched", cond), ("mismatched", other)):
                s = ag2_scores(r.variant, c)
                hit = s[s.gene == gene]
                if not hit.empty:
                    rows.append({"cell_type": prefix, "condition": c, "kind": kind, "variant": r.variant, "gene": gene,
                                 "beta": r.beta, "rna_seq": float(hit.rna_seq.iloc[0]),
                                 "sc_expression": float(hit.sc_expression.iloc[0])})
    m = pd.DataFrame(rows)
    m.to_parquet(PROCESSED / "alphagenome_validation_ag2.parquet", index=False)
    out = {}
    for score in ("sc_expression", "rna_seq"):
        for kind, d in m.groupby("kind"):
            d = d.assign(agree=np.sign(d[score]) == np.sign(d.beta), abs=d[score].abs())
            top = d[d["abs"] >= d["abs"].quantile(0.75)]
            out[f"{score}_{kind}"] = {"pairs": int(len(d)), "agreement_all": round(float(d.agree.mean()), 3),
                                      "agreement_top_quarter": round(float(top.agree.mean()), 3),
                                      "threshold_top_quarter": float(d["abs"].quantile(0.75)),
                                      "spearman": round(float(d[score].corr(d.beta, method="spearman")), 3)}
    return out


# ------------------------------------------------------------------ UC loci: predicted effect of the risk allele
AG2_PIP_MIN, AG2_PER_LOCUS = 0.05, 5


def interest_genes() -> set[str]:
    """Genes of the three target lists, the Candidates top 100 and the UC clinical targets."""
    out = set()
    for f in ("validation_known", "validation_shortlist", "validation_repurposing"):
        p = PROCESSED / f"{f}.parquet"
        if p.exists():
            out |= set(pd.read_parquet(p, columns=["gene"]).gene)
    try:
        from .runs import current_run
        out |= {r["gene"] for r in json.loads((current_run() / "convergence.json").read_text(encoding="utf-8"))["ranking"][:100]}
        from .evaluate import clinical_targets
        out |= set(clinical_targets().gene)
    except (FileNotFoundError, KeyError):
        pass
    return out


def _weighted(d: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """PIP-weighted mean of the risk-signed scores per locus and gene, with the variant that contributes most."""
    w = d.pip / d.groupby(["locus", "gene"]).pip.transform("sum")
    out = d[["locus", "gene"]].copy()
    for c in cols:
        out[c] = d[c] * d.risk_sign * w
    agg = out.groupby(["locus", "gene"], as_index=False)[cols].sum()
    top = d.assign(_c=d.pip * d[cols].abs().max(axis=1)).sort_values("_c").groupby(["locus", "gene"]).tail(1)
    return agg.merge(top[["locus", "gene", "variant", "rsid", "pip"]].rename(
        columns={"variant": "top_variant", "rsid": "top_rsid", "pip": "top_pip"}), on=["locus", "gene"])


def run() -> dict:
    g = uc_loci()
    cs = g[g.in_cs]
    val_atlas = validate_atlas()
    val_ag2 = validate_ag2()
    thr_atlas = val_atlas["threshold_top_quarter"]
    thr_ag2 = val_ag2["sc_expression_matched"]["threshold_top_quarter"]
    s = atlas_scores(cs.variant.tolist())
    m = s.merge(cs[["variant", "locus", "pip", "risk_sign", "rsid"]], on="variant")
    groups = list(ATLAS_GROUPS)
    atlas = _weighted(m, groups).melt(id_vars=["locus", "gene", "top_variant", "top_rsid", "top_pip"], value_vars=groups,
                                      var_name="context", value_name="effect")
    atlas["source"], atlas["lineage"] = "Atlas", atlas.context.map({"colon tissue": "tissue", "blood T cells": "T/NK",
                                                                     "blood monocytes": "myeloid"})
    atlas["confident"] = atlas.effect.abs() >= thr_atlas
    # AlphaGenome 2 in colon cell types: loci near genes of interest, top credible-set variants
    genes = interest_genes()
    near = set(m.loc[m.gene.isin(genes), "locus"])
    sel = cs[cs.locus.isin(near) & (cs.pip >= AG2_PIP_MIN)].sort_values("pip", ascending=False).groupby("locus").head(AG2_PER_LOCUS)
    rows = []
    for r in sel.itertuples():
        for cond in AG2_CONDITIONS:
            x = ag2_scores(r.variant, cond)
            x = x[x.gene.isin(genes)]
            rows.append(x.assign(locus=r.locus, pip=r.pip, risk_sign=r.risk_sign, rsid=r.rsid))
    ag = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    parts = [atlas]
    if not ag.empty:
        for cond, d in ag.groupby("condition"):
            e = _weighted(d, ["sc_expression"]).rename(columns={"sc_expression": "effect"})
            e["context"], e["lineage"] = AG2_CONDITIONS[cond]
            e["source"], e["confident"] = "AlphaGenome 2", e.effect.abs() >= thr_ag2
            parts.append(e)
    out = pd.concat(parts, ignore_index=True)
    lead = g[g.lead].drop_duplicates("locus").set_index("locus")
    out["lead_rsid"] = out.locus.map(lead.rsid)
    out["lead_p"] = out.locus.map(lead.p)
    out["locus_name"] = out.locus.map(lambda l: f"chr{lead.at[l, 'chrom']}:{lead.at[l, 'pos']:,}")
    out["of_interest"] = out.gene.isin(genes)
    out["direction"] = np.where(out.effect > 0, "risk allele raises expression", "risk allele lowers expression")
    out.to_parquet(PROCESSED / "alphagenome_uc.parquet", index=False)
    conf = out[out.confident]
    gi = conf[conf.of_interest]
    s_ = {"loci": int(g.locus.nunique()), "credible_set_variants": int(len(cs)),
          "validation_atlas": val_atlas, "validation_ag2": val_ag2,
          "loci_with_confident_effect": int(conf.locus.nunique()),
          "genes_with_confident_effect": int(conf.gene.nunique()),
          "genes_of_interest_with_confident_effect": sorted(gi.gene.unique().tolist()),
          "ag2_variants": int(len(sel)), "ag2_predictions": int(len(sel) * len(AG2_CONDITIONS))}
    (PROCESSED / "alphagenome_summary.json").write_text(json.dumps(s_, indent=1), encoding="utf-8")
    write_manifest()
    return s_


def load(store=None) -> dict:
    """One alphagenome_variant_effect record per gene with a confident prediction (new evidence type)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "alphagenome_uc.parquet")
    d = d[d.confident]
    n = 0
    with store.batch():
        if store.get_node("DATASET:ALPHAGENOME") is None:
            store.add_node("DATASET:ALPHAGENOME", "Dataset", "AlphaGenome 2 and AlphaGenome Atlas variant-effect "
                                                             "predictions for UC GWAS loci (non-commercial use)")
        for gene, x in d.groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            x = x.reindex(x.effect.abs().sort_values(ascending=False).index).head(8)
            rec = {"gene": gene, "model": "AlphaGenome (Atlas RNA-seq; AlphaGenome 2 single-cell expression)",
                   "variants": "Liu 2023 UC GWAS, approximate fine-mapping (95 % credible sets)",
                   "effects": [{"locus": r.locus_name, "lead": r.lead_rsid, "context": r.context, "source": r.source,
                                "effect_of_risk_allele": round(float(r.effect), 4), "direction": r.direction,
                                "top_variant": r.top_rsid, "top_pip": round(float(r.top_pip), 3)} for r in x.itertuples()],
                   "interpretation": "predicted from DNA sequence, not measured; confident = above the top-quarter "
                                     "|score| of the validation against fine-mapped eQTL; non-commercial research use"}
            ev = store.add_evidence(evidence_id("alphagenome_variant_effect", gene), "alphagenome_variant_effect",
                                    "DATASET:ALPHAGENOME", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "predicted_regulated_by_uc_variant_in", "DATASET:ALPHAGENOME", ev)
            n += 1
    return {"alphagenome_variant_effect_records": n}


def write_manifest() -> dict:
    """sha256 of every cached AlphaGenome answer (data/raw/alphagenome/manifest.json)."""
    files = {}
    for f in sorted(list(CACHE.rglob("*.parquet")) + list(CACHE.glob("flipped_variants.json"))):
        files[f.relative_to(CACHE).as_posix()] = {"sha256": hashlib.sha256(f.read_bytes()).hexdigest(),
                                                  "bytes": f.stat().st_size}
    m = {"source": "AlphaGenome API (alphagenome 2.0.0: Atlas RNA_SEQ scores and AlphaGenome 2 live predictions); "
                   "outputs for non-commercial research only (AlphaGenome Terms of Service)",
         "checked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "files": files}
    (CACHE / "manifest.json").write_text(json.dumps(m, indent=1, sort_keys=True), encoding="utf-8")
    return {"files": len(files)}
