"""Analyses added on 2026-11-10 so that the report covers every view of the kidney atlas, adapted to UC.

All deterministic, from data/processed; designs fixed before looking at the results (evals/history.md):

- sensitivity(): rank stability of the convergence score and the two profiles under 1,000 random re-weightings
  (each weight multiplied by a log-uniform factor in [0.5, 2]; profiles.sensitivity).
- subgroups(): does the answer change with the patient? (kidney: tumour genotype and ccA/ccB)
  (a) age: genes up (padj < 0.05, log2FC >= 1) in both paediatric RNA-seq cohorts (PROTECT, RISK) but in at most
      one of the five adult cohorts, and the reverse (up in >= 4 of 5 adult cohorts, in neither paediatric cohort);
  (b) site: MSCCR inflamed vs control biopsies, rectum and proximal colon separately (patient means, Welch t-test,
      BH); a site-specific gene is up (padj < 0.05, log2FC >= 1) at one site and not significant at the other;
  (c) histological activity axis: per MSCCR biopsy, the neutrophil lineage score (composition markers) by endoscopic
      Mayo subscore, the UC counterpart of the ccA/ccB score.
- biologics(): suitability for an antibody or other biologic (kidney: ADC and radioligand suitability):
      surface or secreted (Open Targets antibody tractability buckets: UniProt location high confidence 1, GO cellular
      component high confidence or UniProt location medium 0.7, signal peptide or transmembrane 0.5) x 0.30
      + increase in inflamed mucosa (meta log2FC / 3, capped) x 0.25
      + low expression in vital organs (1 - HPA consensus nTPM / 100, floored at 0) x 0.25
      + share of patients carrying it (inflamed MSCCR patients above the 95th percentile of controls; 50 % = full) x 0.20,
  over genes up in UC (meta padj < 0.05, log2FC >= 1, measured in >= 2 cohorts). Controls: antibody targets of UC
  drugs (TNF, ITGA4, ITGB7, MADCAM1, IL12B, IL23A, TNFSF15).
- immune(): immune microenvironment. MSCCR patients split into thirds of an immune score (mean of the T/NK, myeloid,
  neutrophil, B and plasma lineage scores); lineage scores per third. Immune-target score per gene (0-1 each):
      immune share of transcripts in inflamed single cells (SCP3755) x 0.3 + coupling: max Spearman rho with the T/NK,
      myeloid or neutrophil score across patients / 0.6 x 0.25 + up in UC (meta log2FC / 2) x 0.2
      + predicted tractability x 0.15 + response: higher in ustekinumab non-responders (t / 4) x 0.1.
  Controls: immune targets of UC drugs (JAK1, TYK2, IL23A, IL12B, ITGB7, S1PR1, TNF).
- novelty(): an ordering score for the discovery pool (kidney: novelty score), used only to order the table and the
  cards; status filters are unchanged.
- candidate_evidence(): per candidate, which evidence layers are measured, meet the criterion, or contradict it.
- score_sensitivity(): rank stability of the biologic, immune, novelty and pocket scores under random re-weighting.
Outputs: data/processed/report_layers.json and a few parquet tables (see each function).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED, RAW

PAED = ["GSE109142", "GSE117993"]
ADULT = ["GSE193677", "GSE206285", "GSE92415", "GSE73661", "GSE16879"]
BIOLOGIC_CONTROLS = ["TNF", "ITGA4", "ITGB7", "MADCAM1", "IL12B", "IL23A", "TNFSF15"]
IMMUNE_CONTROLS = ["JAK1", "TYK2", "IL23A", "IL12B", "ITGB7", "S1PR1", "TNF"]
IMMUNE_LINEAGES = ["T/NK", "myeloid", "neutrophil", "B", "plasma"]


def _json_default(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return None if np.isnan(x) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    raise TypeError(type(x))


# ------------------------------------------------------------------ rank stability
def sensitivity() -> dict:
    from .config import load
    from .profiles import sensitivity as sens
    cfg = load()
    t = pd.read_parquet(PROCESSED / "convergence_table.parquet")
    out = {}
    weights = {"convergence": cfg["convergence"]["stream_weights"], **cfg.get("profiles", {})}
    for name, w in weights.items():
        if not all(f"c_{k}" in t for k in w):
            continue
        d = sens(t, w)
        base = (sum(v * t[f"c_{k}"] for k, v in w.items())).rank(ascending=False, method="first").astype(int)
        d["base_rank"] = d["gene"].map(dict(zip(t["gene"], base)))
        d = d.sort_values("base_rank")
        d.to_parquet(PROCESSED / f"sensitivity_{name}.parquet", index=False)
        out[name] = {"top10_stable": int((d.head(10)["p_top10"] >= 0.8).sum())}
    return out


# ------------------------------------------------------------------ subgroups
def _patient_de(region_mask: pd.Series, x: pd.DataFrame, s: pd.DataFrame) -> pd.DataFrame:
    s = s[region_mask]
    s = s[s["group"].isin(["UC_inflamed", "healthy"])]
    pm = x[s["sample"]].T.groupby(s["patient"].to_numpy() + "|" + s["group"].to_numpy()).mean().T
    infl = np.array([c.endswith("|UC_inflamed") for c in pm.columns])
    a, b = pm.loc[:, infl].to_numpy(float), pm.loc[:, ~infl].to_numpy(float)
    t = stats.ttest_ind(a, b, axis=1, equal_var=False)
    d = pd.DataFrame({"gene": pm.index, "log2fc": a.mean(axis=1) - b.mean(axis=1), "p": t.pvalue,
                      "n_inflamed": infl.sum(), "n_control": (~infl).sum()}).dropna()
    d["padj"] = multipletests(d["p"], method="fdr_bh")[1]
    return d


def subgroups() -> dict:
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    de = de[de.contrast.isin(["uc_vs_control", "inflamed_vs_control"])]
    up = de.assign(up=(de.padj < 0.05) & (de.log2fc >= 1)).pivot_table(index="gene", columns="cohort", values="up",
                                                                       aggfunc="max").fillna(False)
    lfc = de.pivot_table(index="gene", columns="cohort", values="log2fc")
    paed_up = up[PAED].all(axis=1) & (up[ADULT].sum(axis=1) <= 1)
    adult_up = (up[ADULT].sum(axis=1) >= 4) & ~up[PAED].any(axis=1)
    counts = {"paediatric": {c: int(up[c].sum()) for c in PAED}, "adult": {c: int(up[c].sum()) for c in ADULT}}
    # site within MSCCR
    x = pd.read_parquet(PROCESSED / "expr_GSE193677.parquet")
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    s = s[(s.cohort == "GSE193677") & (s.timepoint == "baseline") & s["sample"].isin(x.columns)]
    proximal = {"RightColon", "Cecum", "Transverse", "Colon"}
    rec = _patient_de(s["region"] == "Rectum", x, s)
    prox = _patient_de(s["region"].isin(proximal), x, s)
    site = rec.merge(prox, on="gene", suffixes=("_rectum", "_proximal"))
    site.to_parquet(PROCESSED / "subgroup_site.parquet", index=False)
    sig_r = (site.padj_rectum < 0.05) & (site.log2fc_rectum >= 1)
    sig_p = (site.padj_proximal < 0.05) & (site.log2fc_proximal >= 1)
    rect_only = site[sig_r & (site.padj_proximal >= 0.05)]
    prox_only = site[sig_p & (site.padj_rectum >= 0.05)]
    # activity axis: neutrophil score by Mayo
    mk = json.loads((PROCESSED / "composition_markers.json").read_text())
    s2 = s[s.group.isin(["UC_inflamed", "UC_noninflamed", "healthy"])].copy()
    sub = x[s2["sample"]].reindex(mk["neutrophil"]).dropna()
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    s2["neutrophil_score"] = z.mean(axis=0).reindex(s2["sample"]).to_numpy()
    s2["mayo_group"] = np.where(s2.group == "healthy", "control", s2["mayo"].map(
        lambda m: None if pd.isna(m) else f"Mayo {int(m)}"))
    act = s2.dropna(subset=["mayo_group", "neutrophil_score"])[["sample", "patient", "mayo_group", "neutrophil_score"]]
    act.to_parquet(PROCESSED / "subgroup_activity_axis.parquet", index=False)
    uc = act[act.mayo_group != "control"]
    rho = stats.spearmanr(uc["mayo_group"].str[-1].astype(int), uc["neutrophil_score"])
    tract = _tractable()
    clin = _clinical_genes()

    def rows(genes, extra):
        out = []
        for g in genes:
            r = {"gene": g, "predicted_tractable": bool(tract.get(g, False)), "clinical_drug": g in clin}
            r.update(extra(g))
            out.append(r)
        return out
    paed_tab = lfc[paed_up].assign(m=lfc.loc[paed_up, PAED].mean(axis=1)).sort_values("m", ascending=False).head(12)
    adult_tab = lfc[adult_up].assign(m=lfc.loc[adult_up, ADULT].mean(axis=1)).sort_values("m", ascending=False).head(12)
    out = {
        "counts": counts, "site_counts": {"rectum": int(sig_r.sum()), "proximal": int(sig_p.sum()),
                                          "rectum_only": len(rect_only), "proximal_only": len(prox_only),
                                          "n_rectum": [int(rec.n_inflamed.iloc[0]), int(rec.n_control.iloc[0])],
                                          "n_proximal": [int(prox.n_inflamed.iloc[0]), int(prox.n_control.iloc[0])],
                                          "rho_log2fc": round(float(stats.spearmanr(site.log2fc_rectum, site.log2fc_proximal).statistic), 3)},
        "age_counts": {"paediatric_only": int(paed_up.sum()), "adult_only": int(adult_up.sum())},
        "activity_axis": {"spearman_rho": round(float(rho.statistic), 3), "p": float(f"{rho.pvalue:.3g}"), "n": len(uc)},
        "tables": {
            "Paediatric UC (PROTECT and RISK)": rows(paed_tab.index, lambda g: {"a": round(float(lfc.at[g, "GSE109142"]), 2),
                                                                             "b": round(float(lfc.at[g, "GSE117993"]), 2),
                                                                             "c": round(float(lfc.loc[g, ADULT].mean()), 2)}),
            "Adult UC (five adult cohorts)": rows(adult_tab.index, lambda g: {"a": round(float(lfc.loc[g, ADULT].mean()), 2),
                                                                           "b": round(float(lfc.at[g, "GSE193677"]), 2),
                                                                           "c": round(float(lfc.loc[g, PAED].mean()), 2)}),
            "Rectum only (MSCCR)": rows(rect_only.sort_values("log2fc_rectum", ascending=False).gene.head(12),
                                        lambda g: {"a": round(float(site.set_index("gene").at[g, "log2fc_rectum"]), 2),
                                                   "b": float(f"{site.set_index('gene').at[g, 'padj_rectum']:.2g}"),
                                                   "c": round(float(site.set_index("gene").at[g, "log2fc_proximal"]), 2)}),
            "Proximal colon only (MSCCR)": rows(prox_only.sort_values("log2fc_proximal", ascending=False).gene.head(12),
                                                lambda g: {"a": round(float(site.set_index("gene").at[g, "log2fc_proximal"]), 2),
                                                           "b": float(f"{site.set_index('gene').at[g, 'padj_proximal']:.2g}"),
                                                           "c": round(float(site.set_index("gene").at[g, "log2fc_rectum"]), 2)}),
        },
        "table_columns": {"Paediatric UC (PROTECT and RISK)": ["PROTECT log2FC", "RISK log2FC", "adult mean"],
                          "Adult UC (five adult cohorts)": ["adult mean log2FC", "MSCCR log2FC", "paediatric mean"],
                          "Rectum only (MSCCR)": ["rectum log2FC", "rectum padj", "proximal log2FC"],
                          "Proximal colon only (MSCCR)": ["proximal log2FC", "proximal padj", "rectum log2FC"]},
    }
    return out


# ------------------------------------------------------------------ helpers
def _tractable() -> dict[str, bool]:
    t = pd.read_parquet(PROCESSED / "convergence_table.parquet", columns=["gene", "c_tractability_predicted"])
    return dict(zip(t.gene, t.c_tractability_predicted > 0))


def _clinical_genes() -> set[str]:
    t = pd.read_parquet(PROCESSED / "convergence_table.parquet", columns=["gene", "c_tractability"])
    return set(t.loc[t.c_tractability >= 1, "gene"])


def _hgnc_ensembl() -> dict[str, str]:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    return dict(zip(hg.ensembl_gene_id, hg.symbol))


def _vital() -> pd.Series:
    from .discovery import vital_expression
    return vital_expression()


def _patients_lineages() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    from .composition import _patient_matrix, _zscores
    pm, infl = _patient_matrix()
    mk = json.loads((PROCESSED / "composition_markers.json").read_text())
    sc = pd.DataFrame({lin: _zscores(pm, g).mean(axis=0) for lin, g in mk.items()})
    return pm, sc, infl.to_numpy().astype(bool)


# ------------------------------------------------------------------ biologics
def biologics() -> dict:
    meta = pd.read_parquet(PROCESSED / "meta.parquet")
    up = meta[(meta.meta_padj < 0.05) & (meta.meta_log2fc >= 1) & (meta.k >= 2)]
    genes = list(dict.fromkeys(up.gene.tolist() + BIOLOGIC_CONTROLS))
    parts = sorted((RAW / "extra").glob("ot_target_tractability_*.parquet"))
    tr = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    tr = tr[(tr.modality == "AB") & tr.value.astype(bool)]
    sym = _hgnc_ensembl()
    tr["gene"] = tr.targetId.map(sym)
    loc_score = {"UniProt loc high conf": 1.0, "GO CC high conf": 0.7, "UniProt loc med conf": 0.7,
                 "UniProt SigP or TMHMM": 0.5}
    surf = tr.assign(s=tr.category.map(loc_score).fillna(0.0)).groupby("gene")["s"].max()
    vital = _vital()
    pm, _, infl = _patients_lineages()
    rows = []
    for g in genes:
        m = meta.set_index("gene")
        lf = float(m.at[g, "meta_log2fc"]) if g in m.index else np.nan
        carry = np.nan
        if g in pm.index:
            v = pm.loc[g].to_numpy(float)
            thr = np.quantile(v[~infl], 0.95)
            carry = float((v[infl] > thr).mean())
        vit = float(vital.get(g)) if g in vital.index else np.nan
        comp = {"surface": float(surf.get(g, 0.0) or 0.0),
                "increase": float(np.clip(lf / 3, 0, 1)) if pd.notna(lf) else 0.0,
                "vital_low": float(np.clip(1 - vit / 100, 0, 1)) if pd.notna(vit) else 0.5,
                "carriers": float(np.clip(carry / 0.5, 0, 1)) if pd.notna(carry) else 0.0}
        w = {"surface": 0.30, "increase": 0.25, "vital_low": 0.25, "carriers": 0.20}
        rows.append({"gene": g, "suitability": round(sum(w[k] * comp[k] for k in w), 3), **{f"c_{k}": round(v, 3) for k, v in comp.items()},
                     "meta_log2fc": None if pd.isna(lf) else round(lf, 2), "vital_max_ntpm": None if pd.isna(vit) else round(vit, 1),
                     "patients_carrying_pct": None if pd.isna(carry) else round(100 * carry, 1),
                     "control": g in BIOLOGIC_CONTROLS, "in_universe": g in set(up.gene)})
    d = pd.DataFrame(rows).sort_values("suitability", ascending=False).reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    d.to_parquet(PROCESSED / "biologics.parquet", index=False)
    return {"candidates": int(d.in_universe.sum()), "surface": int((d.c_surface > 0).sum()),
            "controls": {g: {"rank": int(d.loc[d.gene == g, "rank"].iloc[0]), "in_universe": bool(d.loc[d.gene == g, "in_universe"].iloc[0])}
                         for g in BIOLOGIC_CONTROLS if (d.gene == g).any()},
            "weights": {"surface": 0.30, "increase": 0.25, "vital_low": 0.25, "carriers": 0.20}}


# ------------------------------------------------------------------ immune microenvironment
def immune() -> dict:
    pm, sc, infl = _patients_lineages()
    lins = [x for x in IMMUNE_LINEAGES if x in sc]
    score = sc[lins].mean(axis=1)
    cls = pd.qcut(score, 3, labels=["low", "intermediate", "high"])
    by = sc[lins].groupby(cls.to_numpy(), observed=True).mean()
    classes = {str(k): int(v) for k, v in cls.value_counts().items()}
    infl_by = pd.Series(infl, index=sc.index).groupby(cls.to_numpy(), observed=True).mean()
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet", columns=["gene", "health", "compartment", "share"])
    sh = sh[sh.health == "inflamed"]
    immune_comp = {"B Cell", "Plasma", "T Cell and NK", "Macrophage and DC", "Neutrophil", "Mast Cell", "Eosinophil"}
    ishare = sh[sh.compartment.isin(immune_comp)].groupby("gene")["share"].sum()
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    tract = _tractable()
    rg = pd.read_parquet(PROCESSED / "response_genes.parquet")
    ust = rg[rg.arm == "ustekinumab"].set_index("gene")["t"]
    cand = meta[(meta.meta_padj < 0.05) & (meta.k >= 2)].index.intersection(pm.index)
    cand = list(dict.fromkeys(list(cand) + [g for g in IMMUNE_CONTROLS if g in pm.index]))
    X = pm.loc[cand].to_numpy(float)
    rk = pd.DataFrame(X).rank(axis=1).to_numpy()
    rows = []
    lin_r = {lin: stats.rankdata(sc[lin].to_numpy()) for lin in ("T/NK", "myeloid", "neutrophil") if lin in sc}

    def rho(a, b):
        a, b = a - a.mean(), b - b.mean()
        return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))
    for i, g in enumerate(cand):
        rhos = {lin: rho(rk[i], r) for lin, r in lin_r.items()}
        best = max(rhos.values())
        comp = {"immune_share": float(np.clip(ishare.get(g, 0.0), 0, 1)),
                "coupling": float(np.clip(best / 0.6, 0, 1)),
                "up": float(np.clip(meta.at[g, "meta_log2fc"] / 2, 0, 1)) if g in meta.index else 0.0,
                "tractable": float(tract.get(g, False)),
                "response": float(np.clip(ust.get(g, 0.0) / 4, 0, 1))}
        w = {"immune_share": 0.3, "coupling": 0.25, "up": 0.2, "tractable": 0.15, "response": 0.1}
        rows.append({"gene": g, "immune_target_score": round(sum(w[k] * comp[k] for k in w), 3),
                     **{f"c_{k}": round(v, 3) for k, v in comp.items()},
                     **{f"rho_{k.replace('/', '_')}": round(v, 3) for k, v in rhos.items()},
                     "meta_log2fc": round(float(meta.at[g, "meta_log2fc"]), 2) if g in meta.index else None,
                     "control": g in IMMUNE_CONTROLS})
    d = pd.DataFrame(rows).sort_values("immune_target_score", ascending=False).reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    d.to_parquet(PROCESSED / "immune_targets.parquet", index=False)
    pd.DataFrame({"patient": sc.index, "immune_score": score.to_numpy(), "class": cls.astype(str).to_numpy(),
                  "inflamed": infl, **{lin: sc[lin].to_numpy() for lin in lins}}).to_parquet(
        PROCESSED / "immune_classes.parquet", index=False)
    ctrl = {g: int(d.loc[d.gene == g, "rank"].iloc[0]) for g in IMMUNE_CONTROLS if (d.gene == g).any()}
    return {"patients": int(len(sc)), "classes": classes, "inflamed_share_by_class": {str(k): round(float(v), 3) for k, v in infl_by.items()},
            "lineage_means_by_class": {str(k): {lin: round(float(v), 3) for lin, v in r.items()} for k, r in by.iterrows()},
            "universe": len(d), "controls": ctrl, "controls_in_top100": int(sum(1 for r in ctrl.values() if r <= 100)),
            "weights": {"immune_share": 0.3, "coupling": 0.25, "up": 0.2, "tractable": 0.15, "response": 0.1}}


# ------------------------------------------------------------------ discovery ordering
NOVELTY_WEIGHTS = {"increase": 0.3, "tissue_change": 0.2, "understudied": 0.2, "druggable": 0.15, "safety": 0.15}


def novelty() -> dict:
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    d["c_increase"] = (d.meta_log2fc.fillna(0) / 3).clip(0, 1)
    d["c_tissue_change"] = d.retained.fillna(0).clip(0, 1)
    d["c_understudied"] = (1 - np.log10(d.uc_papers.fillna(100) + 1) / np.log10(101)).clip(0, 1)
    d["c_druggable"] = d.tractable.map(lambda x: 1.0 if len(x) and ("SM" in list(x)) else (0.7 if len(x) else 0.0))
    d["c_safety"] = (1 - d.vital_max_ntpm.fillna(50) / 100).clip(0, 1)
    d["novelty_score"] = sum(w * d[f"c_{k}"] for k, w in NOVELTY_WEIGHTS.items()).round(3)
    lin_track = {"epithelial": "epithelial", "stromal": "stromal", "endothelial": "stromal"}
    d["track"] = d.lineage.map(lambda x: lin_track.get(x, "immune"))
    d["rank_in_track"] = d.groupby("track")["novelty_score"].rank(ascending=False, method="first").astype(int)
    d.to_parquet(PROCESSED / "discovery_scored.parquet", index=False)
    return {"by_track": d.track.value_counts().to_dict(),
            "fic_by_track": d[d.status == "first-in-class opportunity"].track.value_counts().to_dict(),
            "weights": NOVELTY_WEIGHTS}


# ------------------------------------------------------------------ evidence coverage per candidate
def candidate_evidence(genes: list[str]) -> list[dict]:
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    paired = de[de.contrast == "inflamed_vs_noninflamed_paired"].set_index("gene")
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    act = pd.read_parquet(PROCESSED / "activity.parquet")
    act = act[act.cohort == "GSE193677"].set_index("gene")
    mr = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").set_index("gene")
    mrg = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    mrc = mrg[mrg.layer == "colon_tissue"].set_index("gene")
    vis = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())
    vis = {t["gene"]: t for t in vis["targets"]}
    disc = pd.read_parquet(PROCESSED / "discovery.parquet").set_index("gene")
    conv = pd.read_parquet(PROCESSED / "convergence_table.parquet").set_index("gene")
    rg = pd.read_parquet(PROCESSED / "response_genes.parquet")
    rg = rg[rg.arm.isin(["ustekinumab", "placebo (UNIFI)"])].pivot_table(index="gene", columns="arm", values="padj")
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet", columns=["gene", "health", "compartment", "share"])
    sh = sh[sh.health == "inflamed"].groupby(["gene", "compartment"])["share"].sum()
    st = json.loads((PROCESSED / "structure_summary.json").read_text()) if (PROCESSED / "structure_summary.json").exists() else {"results": []}
    st = {r["gene"]: r for r in st["results"]}
    cards = []
    for g in genes:
        obs = []

        def add(name, status, detail, source, **fields):
            obs.append({"name": name, "status": status, "detail": detail, "source": source,
                        "fields": {k: v for k, v in fields.items() if v is not None}})
        if g in meta.index:
            m = meta.loc[g]
            ok = m.meta_padj < 0.05 and abs(m.meta_log2fc) >= 1 and m.k >= 2
            add("Meta-analysis, eight cohorts", "meets" if ok else "measured_not_meeting",
                "pooled change in UC vs control" + ("" if m.k >= 2 else " (one cohort only: no pooled test)"), "meta.parquet",
                log2fc=round(float(m.meta_log2fc), 2), padj=float(f"{m.meta_padj:.2g}") if pd.notna(m.meta_padj) else None,
                cohorts=int(m.k), significant_up=int(m.n_sig_up), significant_down=int(m.n_sig_down))
        else:
            add("Meta-analysis, eight cohorts", "missing", "not measured", "meta.parquet")
        if g in paired.index:
            p = paired.loc[g]
            add("Inflamed vs non-inflamed, same patient", "meets" if p.padj < 0.05 and abs(p.log2fc) >= 0.5 else "measured_not_meeting",
                "local inflammation rather than disease-wide", "de_results.parquet", log2fc=round(float(p.log2fc), 2),
                padj=float(f"{p.padj:.2g}"))
        if g in comp.index:
            c = comp.loc[g]
            add("Composition adjustment", {"intrinsic": "meets", "partly_composition": "partly", "composition": "contradicts"}[c.verdict],
                f"share of the bulk change kept after the cell mix is accounted for ({c.verdict.replace('_', ' ')})",
                "composition_effect.parquet", retained=round(float(c.retained), 2), adjusted_log2fc=round(float(c.adjusted_log2fc), 2))
        if g in act.index:
            a = act.loc[g]
            add("Disease activity (Mayo)", "meets" if a.padj < 0.05 and abs(a.rho) >= 0.2 else "measured_not_meeting",
                "correlation with the endoscopic Mayo subscore (MSCCR)", "activity.parquet", rho=round(float(a.rho), 2),
                padj=float(f"{a.padj:.2g}"))
        if g in sh.index.get_level_values(0):
            s = sh.loc[g].sort_values(ascending=False)
            add("Single-cell source (inflamed)", "measured", f"largest share in {s.index[0]}", "sc3_celltype_share.parquet",
                top=s.index[0], share=round(float(s.iloc[0]), 2))
        if g in vis:
            v = vis[g]
            add("Tissue (Visium, deconvolved)", "measured", f"in-situ expresser {v['top_in_situ_expresser']}",
                "spatial_visium_summary.json", contribution=v["in_situ_contribution"],
                neighbours=", ".join(e["cell_group"] for e in v["enriched_neighbors"][:3]) or None)
        else:
            add("Tissue (Visium, deconvolved)", "missing", "not analysed or too few spots", "spatial_visium_summary.json")
        for name, tab, src, unit in (("Genetics: plasma protein MR", mr, "mr_pqtl.parquet", "protein"),
                                     ("Genetics: colon expression MR", mrc, "mr_colon.parquet", "expression")):
            if g in tab.index:
                r = tab.loc[g]
                status = "meets" if r.q < 0.05 else "measured_not_meeting"
                add(name, status, f"genetically higher {unit} {'raises' if r.beta > 0 else 'lowers'} UC risk" if r.q < 0.05
                    else "no significant causal estimate", src, odds_ratio=round(float(r.or_per_sd), 2), q=float(f"{r.q:.2g}"))
            else:
                add(name, "missing", "no genetic instrument (not tested)", src)
        if g in conv.index:
            c = conv.loc[g]
            add("Open Targets genetic association", "meets" if c.ot_genetic >= 0.5 else "measured_not_meeting",
                "germline association score with UC", "convergence_table.parquet", score=round(float(c.ot_genetic), 2))
            add("Tractability", "meets" if c.c_tractability_predicted > 0 else "measured_not_meeting",
                "predicted small-molecule or antibody tractability", "convergence_table.parquet",
                clinical_drug=bool(c.c_tractability >= 1))
        if g in rg.index:
            r = rg.loc[g]
            add("Treatment response (UNIFI)", "measured", "differs between ustekinumab non-responders and responders"
                if r.get("ustekinumab", 1) < 0.05 else "no difference by response", "response_genes.parquet",
                ustekinumab_padj=float(f"{r.get('ustekinumab', np.nan):.2g}"),
                placebo_padj=float(f"{r.get('placebo (UNIFI)', np.nan):.2g}"))
        if g in st and st[g].get("modeled"):
            add("Structure (AlphaFold pockets)", "measured", "best pocket score of the predicted monomer",
                "structure_summary.json", pocket_score=st[g]["best_pocket_score"], plddt=st[g]["mean_plddt"])
        if g in disc.index:
            add("Discovery screen", "meets" if disc.at[g, "status"] == "first-in-class opportunity" else "measured",
                disc.at[g, "status"], "discovery.parquet", uc_papers=int(disc.at[g, "uc_papers"]))
        measured = sum(o["status"] != "missing" for o in obs)
        mods = []
        if g in conv.index:
            c = conv.loc[g]
            up = g in meta.index and meta.at[g, "meta_log2fc"] > 0
            if c.c_tractability_predicted > 0 and up:
                mods.append({"name": "Inhibitor (small molecule or antibody)", "status": "plausible" if g not in comp.index
                             or comp.at[g, "verdict"] != "composition" else "caution",
                             "rationale": "up in UC and predicted tractable" + ("; but the increase is mostly cell influx"
                                                                                if g in comp.index and comp.at[g, "verdict"] == "composition" else "")})
            if g in mr.index and mr.at[g, "q"] < 0.05 and mr.at[g, "beta"] < 0:
                mods.append({"name": "Agonist or protein replacement", "status": "genetics favours activation",
                             "rationale": "genetically higher plasma protein lowers UC risk"})
            if g in meta.index and abs(meta.at[g, "meta_log2fc"]) >= 1 and meta.at[g, "k"] >= 4:
                mods.append({"name": "Tissue biomarker", "status": "plausible",
                             "rationale": f"changes in {int(meta.at[g, 'k'])} cohorts"})
        cards.append({"gene": g, "measured": measured, "total": len(obs), "observations": obs, "modalities": mods})
    return cards


# ------------------------------------------------------------------ stability of the newer scores
def score_sensitivity() -> dict:
    """Rank stability of the biologic, immune, novelty and pocket scores under the same 1,000 random re-weightings
    as the convergence score (each weight x log-uniform [0.5, 2]). Outputs sensitivity_<score>.parquet (top 200 by
    the configured weights) and, per score, how many of its configured top 10 stay in the top 10 in >= 80 % of
    re-weightings and the median Spearman correlation between the configured and re-weighted scores."""
    from .profiles import sensitivity as sens
    from .structure import WEIGHTS as POCKET
    tables = {
        "biologics": (pd.read_parquet(PROCESSED / "biologics.parquet"),
                      {"surface": 0.30, "increase": 0.25, "vital_low": 0.25, "carriers": 0.20}),
        "immune": (pd.read_parquet(PROCESSED / "immune_targets.parquet"),
                   {"immune_share": 0.3, "coupling": 0.25, "up": 0.2, "tractable": 0.15, "response": 0.1}),
        "novelty": (pd.read_parquet(PROCESSED / "discovery_scored.parquet"), NOVELTY_WEIGHTS),
    }
    pk = PROCESSED / "structure_pockets.parquet"
    if pk.exists():
        d = pd.read_parquet(pk)
        d = d[d["rank"] == 1].rename(columns={f"comp_{k}": f"c_{k}" for k in POCKET})
        tables["pockets"] = (d, POCKET)
    out = {}
    rng = np.random.default_rng(11)
    for name, (t, w) in tables.items():
        t = t.dropna(subset=[f"c_{k}" for k in w]).reset_index(drop=True)
        r = sens(t, w)
        base = sum(v * t[f"c_{k}"] for k, v in w.items())
        r["base_rank"] = r["gene"].map(dict(zip(t["gene"], base.rank(ascending=False, method="first").astype(int))))
        r = r.sort_values("base_rank")
        r.to_parquet(PROCESSED / f"sensitivity_{name}.parquet", index=False)
        C = t[[f"c_{k}" for k in w]].to_numpy(float)
        w0 = np.array(list(w.values()))
        rhos = [stats.spearmanr(C @ w0, C @ (w0 * np.exp(rng.uniform(np.log(0.5), np.log(2), len(w0))))).statistic
                for _ in range(200)]
        out[name] = {"top10_stable": int((r.head(10)["p_top10"] >= 0.8).sum()), "median_rho": round(float(np.median(rhos)), 3),
                     "n": int(len(t))}
    return out


def run(genes: list[str] | None = None) -> dict:
    out = {"sensitivity": sensitivity(), "subgroups": subgroups(), "biologics": biologics(), "immune": immune(),
           "novelty": novelty()}
    out["score_sensitivity"] = score_sensitivity()
    (PROCESSED / "report_layers.json").write_text(json.dumps(out, indent=1, default=_json_default), encoding="utf-8")
    return out


if __name__ == "__main__":
    r = run()
    print(json.dumps({k: v for k, v in r.items() if k != "subgroups"}, indent=1, default=_json_default))
    print(json.dumps({k: v for k, v in r["subgroups"].items() if k != "tables"}, indent=1, default=_json_default))
