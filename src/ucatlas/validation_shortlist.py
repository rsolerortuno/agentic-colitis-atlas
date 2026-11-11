"""Shortlist for experimental validation (added 2026-11-10, revised 2026-11-11): one row per first-in-class gene, with
a traffic light per criterion and a proposed first experiment. Everything comes from tables already in data/processed;
nothing is re-estimated.

Criteria (green / amber / red; grey = not measured):
- reproducible: cohorts with a significant increase (green >= 6 of 8, amber 4-5, red < 4), and the covariate-adjusted
  MSCCR change holds (if it does not, at most amber);
- protein: the protein also rises in the independent proteomics study (green: replicates; red: measured, does not
  replicate; grey: not quantified);
- tissue change: two views of "is it more RNA per cell, or more cells of the type that makes it?". (a) the bulk
  composition verdict of the discovery screen (intrinsic or partly composition, lineage level); (b) the change of the
  gene inside its own single-cell clusters (Smillie 2019, inflamed vs healthy, mean expression per cell in the three
  clusters that hold most of its RNA, weighted by that share; clusters with >= 50 cells in both states). (b) sees
  shifts inside a lineage that (a) cannot (crypt expansion, loss of mature enterocytes). Green: per-cell change
  >= 1 log2 with an intrinsic verdict, or >= 1.5 with a partly-composition verdict; red: per-cell change < 0.3 (the
  single cells explain the rise by cell numbers); amber otherwise. Calibrated on markers: DUOX2 +7.1, MMP3 +4.9,
  CHI3L1 +4.1 (state changes), TNF +0.2, IL1B -0.5, OSM -0.6 (rise by immune influx), AQP8 -2.7;
- cell source: share of the gene's transcripts in its dominant lineage in inflamed single cells (Smillie 2019, the
  lineage the discovery screen uses; green >= 0.5, amber 0.3-0.5, red < 0.3); a lineage not replicated in the second
  atlas (Kong 2023) caps it at amber. A regeneration flag also caps it at amber: an epithelial gene with >= 76 % of its
  epithelial RNA in the crypt zone (stem, TA, cycling TA, secretory TA, enterocyte progenitors) that rises per cell
  there. 76 % is halfway between EPCAM (58 %, any epithelial gene) and the crypt markers LGR5, OLFM4, MKI67 and REG1A
  (92-97 %). Such a gene may belong to the repair programme of the crypt, and blocking it could slow healing;
- causal genetics: a significant Mendelian-randomisation estimate in plasma or colon in the direction that supports
  inhibition (genetically higher level raises risk) is green, the opposite direction red, tested and null amber (it
  scores 0: a null test is not support), no instrument grey;
- druggable: an antibody route or strong small-molecule evidence in
  Open Targets (high-quality ligand, high-quality pocket, druggable family, clinical precedent) is green; a pocket
  score >= 0.6 alone, a weak bucket (structure with ligand) or a predicted membrane or signal peptide only is amber
  (the pocket score ranked the UC small-molecule targets 4th to 25th of 45 proteins: it does not discriminate alone).
  Antibody route = an extracellular part, from UniProt topology: an extracellular topological domain, a GPI anchor or
  secretion. "Membrane" or "Cell membrane" as a location is not enough (OSBPL3, EHD3 and RFTN1 sit on the cytoplasmic
  face); a protein marked secreted but also resident in the ER membrane (ERO1A) counts as amber. Checked on the UC
  clinical targets: ITGA2, S1PR1, CXCR2, ICAM1, IL2RA, ITGB7 (extracellular domain), IL23A, IL12B, CCL11, CXCL10
  (secreted) pass; JAK3 and PDE4B (small-molecule targets) do not;
- safety: highest vital-organ RNA (> 50 nTPM), loss-of-function constraint (LOEUF < 0.35), known safety events, and a
  Mendelian disease caused by loss of the gene (UniProt); green with none, amber with one flag, red with two or more;
- UC literature: what the papers about the gene in UC or IBD show (shortlist_literature.py: Europe PMC records read by
  a classifier subagent, every call backed by a verbatim sentence). Green: functional evidence of a pathogenic role
  (inhibition should help) and none protective; red: functional evidence of a protective role only (inhibition may
  harm); amber: a specific expression, biomarker or genetic result, or functional evidence both ways; grey: no
  paper about the gene, or only named in gene lists or reviews. Clinical targets in the calibration are not reviewed (grey). Only papers in a colitis or IBD
  setting count; functional results in other inflammatory settings are shown as context;
- mouse: whether the gene has ever been perturbed in mouse colitis (information, not a score).
Points: green 2, amber 1 (0 for a null genetic test), grey 0, red -1; sorted by points, then greens, then meta log2FC.
Calibration: the same rules applied to the UC clinical targets that are in the discovery table
(validation_shortlist_calibration.parquet). They are not expected to score high on every light (several rise by
immune influx, which is how an anti-inflammatory works); the comparison shows what each light separates.
Proposed experiment: by dominant lineage (epithelial: patient-derived colon organoids; stromal: primary colonic
fibroblasts; immune: lamina propria immune cells) and by route (antibody, small molecule or genetic), followed by a
mouse colitis model when none exists.
Outputs: data/processed/validation_shortlist.parquet, validation_shortlist_calibration.parquet,
validation_shortlist_summary.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

ORDER = {"green": 2, "amber": 1, "grey": 0, "red": -1}
CRYPT = {"Stem", "TA 1", "TA 2", "Cycling TA", "Secretory TA", "Enterocyte Progenitors"}
CRYPT_SHARE = 0.76
STRONG_SM = ("High-Quality Ligand", "High-Quality Pocket", "Druggable Family", "Approved Drug", "Advanced Clinical",
             "Phase 1 Clinical")
CRITERIA = ["reproducible", "protein", "tissue", "cell", "genetics", "druggable", "safety", "literature", "mouse"]


def _lineage_share() -> dict[str, tuple[str, float]]:
    # the same source as the lineage of the discovery screen (Smillie 2019, inflamed), so the two always agree
    from .sc_replication import shares_2019
    lin = shares_2019("inflamed").sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    return {g: (r.lineage, float(r.share)) for g, r in lin.iterrows()}


class _SingleCell:
    """Per-cell change inside the gene's own clusters and the crypt-zone share (Smillie 2019)."""

    def __init__(self):
        s = pd.read_parquet(PROCESSED / "sc_celltype_share.parquet")
        w = s.pivot_table(index=["gene", "cluster"], columns="health", values=["mean_cp10k", "n_cells", "share"])
        w.columns = [f"{a}_{b}" for a, b in w.columns]
        self.w = w.reset_index()
        x = s[(s.health == "inflamed") & (s.compartment == "epithelial")]
        self.crypt = (x[x.cluster.isin(CRYPT)].groupby("gene").cp10k_sum.sum() / x.groupby("gene").cp10k_sum.sum())

    def per_cell(self, g: str, k: int = 3) -> tuple[float, str]:
        d = self.w[self.w.gene == g].dropna(subset=["mean_cp10k_healthy", "mean_cp10k_inflamed"])
        d = d[(d.n_cells_healthy >= 50) & (d.n_cells_inflamed >= 50)].nlargest(k, "share_inflamed")
        if d.empty:
            return np.nan, ""
        wt = d.share_inflamed / d.share_inflamed.sum()
        lfc = np.log2((d.mean_cp10k_inflamed + 0.01) / (d.mean_cp10k_healthy + 0.01))
        return float((wt * lfc).sum()), ", ".join(f"{c} {v:+.1f}" for c, v in zip(d.cluster, lfc))


def _topology(symbol: str) -> dict:
    from .external import UNIPROT_URL, _cached, _request
    params = {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
              "fields": "accession,ft_topo_dom,ft_transmem,ft_lipid,cc_subcellular_location", "format": "json", "size": 1}
    d = _cached("uniprot", json.dumps(params, sort_keys=True),
                lambda: _request("GET", UNIPROT_URL, params=params).json())["data"]
    if not d.get("results"):
        return {}
    r = d["results"][0]
    f = r.get("features", [])
    locs = [x["location"]["value"] for c in r.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION"
            for x in c.get("subcellularLocations", []) if "location" in x]
    return {"extracellular_domain": any(x["type"] == "Topological domain" and "xtracellular" in x.get("description", "")
                                        for x in f),
            "gpi": any(x["type"] == "Lipidation" and "GPI" in x.get("description", "") for x in f),
            "secreted": any(x.startswith("Secreted") for x in locs),
            "er_membrane": any(x.startswith("Endoplasmic reticulum membrane") for x in locs)}


def antibody_route(symbol: str) -> tuple[str, str]:
    """('full' | 'partial' | 'none', reason) from UniProt topology."""
    t = _topology(symbol)
    if t.get("extracellular_domain"):
        return "full", "extracellular domain"
    if t.get("gpi"):
        return "full", "GPI-anchored at the surface"
    if t.get("secreted"):
        return ("partial", "secreted, but also ER membrane") if t.get("er_membrane") else ("full", "secreted")
    return "none", ""


def _sm_labels():
    from . import external
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    ens = dict(zip(hg.symbol, hg.ensembl_gene_id))

    def get(g):
        t = external.ot_target(ens[g]) if g in ens else None
        return [x.split(": ", 1)[1] for x in (t or {}).get("tractability_labels", []) if x.startswith("SM: ")]
    return get


def _experiment(lineage: str, surface: bool, sm: bool, mouse_papers: int | None, regen: bool) -> str:
    if lineage == "epithelial":
        system = "knock-down or CRISPR knockout in patient-derived colon organoids stimulated with TNF and IFN-γ " \
                 "(read-outs: barrier integrity, chemokine release, epithelial death)"
        if regen:
            system += ", with a wound-healing or organoid regrowth assay, since the gene rises in regenerating crypt cells"
    elif lineage in ("stromal", "endothelial"):
        system = "silencing in primary colonic fibroblasts from inflamed UC stimulated with OSM and IL-1β " \
                 "(read-outs: IL-6, CXCL1, MMP release; neutrophil attraction in co-culture)"
    else:
        system = "knock-down in lamina propria immune cells or blood-derived equivalents " \
                 "(read-outs: activation and cytokine profile)"
    route = ("a blocking antibody can be tested directly" if surface else
             "small-molecule tool compounds may exist (Open Targets ligand or pocket evidence)" if sm else "genetic tools first")
    mouse = (" Then a knockout or blockade in DSS or T-cell transfer colitis (never tested in mouse colitis)."
             if mouse_papers == 0 else f" Mouse colitis data exist ({mouse_papers} paper{'s' if mouse_papers > 1 else ''}): check them first."
             if mouse_papers else "")
    return f"{system}; {route}.{mouse}"


def _load():
    t = {"adj": pd.read_parquet(PROCESSED / "de_adjusted.parquet").set_index("gene"),
         "prot": pd.read_parquet(PROCESSED / "proteomics.parquet").set_index("gene"),
         "lin": _lineage_share(), "sc": _SingleCell(), "sm": _sm_labels(),
         "saf": pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene"),
         "ml": pd.read_parquet(PROCESSED / "mouse_literature.parquet").set_index("gene")["papers"].to_dict()}
    rep_p = PROCESSED / "sc_replication.parquet"
    t["screp"] = pd.read_parquet(rep_p).set_index("gene")["same_lineage"].to_dict() if rep_p.exists() else {}
    vis = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())["targets"]
    t["vis"] = {v["gene"]: v for v in vis}
    mr_rows = {}
    for f_, lab in (("mr_pqtl.parquet", "plasma"), ("mr_colon.parquet", "colon")):
        r = pd.read_parquet(PROCESSED / f_)
        if "layer" in r:
            r = r[r.layer == "colon_tissue"]
        for x in r.itertuples():
            mr_rows.setdefault(x.gene, []).append((lab, float(x.beta), float(x.q)))
    t["mr"] = mr_rows
    lit_p = PROCESSED / "shortlist_literature.parquet"
    t["lit"] = pd.read_parquet(lit_p).set_index("gene") if lit_p.exists() else pd.DataFrame()
    pk = pd.read_parquet(PROCESSED / "structure_pockets.parquet")
    t["pocket"] = pk[pk["rank"] == 1].set_index("gene")["pocket_score"].to_dict()
    return t


def _row(x, t) -> dict:
    g = x.gene
    c = {}
    rep = "green" if x.cohorts_up >= 6 else "amber" if x.cohorts_up >= 4 else "red"
    holds = bool(t["adj"].at[g, "holds"]) if g in t["adj"].index else None
    if holds is False and rep == "green":
        rep = "amber"
    c["reproducible"] = (rep, f"up in {int(x.cohorts_up)}/8 cohorts" +
                         ("" if holds is None else f"; adjusted: {'holds' if holds else 'does not hold'}"))
    prot = t["prot"]
    if g in prot.index and pd.notna(prot.at[g, "padj"]):
        r_ = prot.at[g, "replicates"]
        c["protein"] = ("green" if r_ == "yes" else "red" if r_ == "no" else "amber",
                        f"log2 ratio {prot.at[g, 'log2_ratio']:.2f}, padj {prot.at[g, 'padj']:.1e}")
    else:
        c["protein"] = ("grey", "not quantified")
    pc, pc_txt = t["sc"].per_cell(g)
    comp = x.composition or "not assessed"
    if np.isnan(pc):
        tis = "amber" if comp in ("intrinsic", "partly_composition") else "red"
    elif pc < 0.3:
        tis = "red"
    elif (comp == "intrinsic" and pc >= 1) or (comp == "partly_composition" and pc >= 1.5):
        tis = "green"
    else:
        tis = "amber"
    ret = f"{float(x.retained):.0%} kept" if pd.notna(x.retained) else ""
    c["tissue"] = (tis, f"bulk: {comp.replace('_', ' ')}{', ' + ret if ret else ''}; per cell: " +
                   (f"{pc:+.1f} log2 ({pc_txt})" if not np.isnan(pc) else "no matched clusters"))
    lname, lshare = t["lin"].get(g, (x.lineage, np.nan))
    v = t["vis"].get(g, {})
    vtxt = f"; in situ {v['top_in_situ_expresser'].replace('_', ' ')}" if v.get("top_in_situ_expresser") else ""
    col = "grey" if np.isnan(lshare) else "green" if lshare >= 0.5 else "amber" if lshare >= 0.3 else "red"
    sr = t["screp"]
    rtxt = "" if g not in sr else "; same lineage in the second atlas" if sr[g] else "; not replicated in the second atlas"
    if col == "green" and sr.get(g) is False:
        col = "amber"
    crypt = float(t["sc"].crypt.get(g, np.nan))
    regen = bool(lname == "epithelial" and crypt >= CRYPT_SHARE and (pc if not np.isnan(pc) else 0) >= 1)
    if regen and col == "green":
        col = "amber"
    ctxt = f"; {crypt:.0%} of epithelial RNA in the crypt zone, rising there: possible repair programme" if regen else ""
    c["cell"] = (col, (f"{lname} {lshare:.0%} of transcripts{rtxt}{ctxt}{vtxt}" if not np.isnan(lshare)
                       else "not in the single-cell atlas"))
    mrs = t["mr"].get(g, [])
    sig = [m for m in mrs if m[2] < 0.05]
    if sig:
        c["genetics"] = ("green" if any(m[1] > 0 for m in sig) else "red",
                         "; ".join(f"{m[0]} OR {np.exp(m[1]):.2f} (q {m[2]:.1e})" for m in sig))
    elif mrs:
        c["genetics"] = ("amber", "tested, no causal effect (scores 0)")
    else:
        c["genetics"] = ("grey", "no genetic instrument")
    ab, ab_txt = antibody_route(g)
    surface = ab == "full"
    sm = t["sm"](g)
    strong = [s for s in sm if s in STRONG_SM]
    weak = [s for s in sm if s not in STRONG_SM]
    pscore = t["pocket"].get(g)
    mods = list(x.tractable) if x.tractable is not None else []
    drug = "green" if surface or strong else "amber" if ab == "partial" or (pscore or 0) >= 0.6 or sm or mods else "red"
    parts = ([f"antibody route: {ab_txt}"] if ab != "none" else ["no extracellular part"]) + (["Open Targets: " + ", ".join(strong)] if strong else []) + \
            ([f"pocket {pscore:.2f}"] if pscore is not None else []) + ([f"weak: {', '.join(weak)}"] if weak else [])
    c["druggable"] = (drug, "; ".join(parts) or ("predicted " + "/".join(mods) if mods else "no route"))
    flags = []
    if pd.notna(x.vital_max_ntpm) and x.vital_max_ntpm > 50:
        flags.append(f"vital organs {x.vital_max_ntpm:.0f} nTPM")
    saf = t["saf"]
    if g in saf.index:
        if saf.at[g, "constrained"] is True:
            flags.append(f"constrained (LOEUF {saf.at[g, 'loeuf']:.2f})")
        if saf.at[g, "n_safety_events"] > 0:
            flags.append("known safety events")
        if "mendelian_loss_disease" in saf.columns and saf.at[g, "mendelian_loss_disease"]:
            flags.append("loss causes " + "; ".join(saf.at[g, "mendelian_diseases"]))
    vital = f"vital organs ≤ {x.vital_max_ntpm:.0f} nTPM" if pd.notna(x.vital_max_ntpm) else "vital organs not measured"
    c["safety"] = ("green" if not flags else "amber" if len(flags) == 1 else "red",
                   "; ".join(flags) if flags else f"{vital}, not constrained, no Mendelian disease")
    lit = t["lit"]
    c["literature"] = ((lit.at[g, "light"], lit.at[g, "verdict"])
                       if g in lit.index else ("grey", "not reviewed"))
    papers = t["ml"].get(g)
    papers = int(papers) if papers is not None and pd.notna(papers) else None
    mouse_lit = lit.at[g, "mouse_colitis"] if g in lit.index and "mouse_colitis" in lit.columns else ""
    if mouse_lit:   # the literature review found a perturbation in mouse colitis that the title search missed
        papers = max(papers or 0, mouse_lit.count("PMID"))
    c["mouse"] = ("grey", "not searched" if papers is None else "never tested in mouse colitis" if papers == 0
                  else f"tested in mouse colitis: {mouse_lit}" if mouse_lit else f"{papers} mouse colitis papers")
    pts = {k: (0 if k == "genetics" and v[0] == "amber" else ORDER[v[0]]) for k, v in c.items() if k != "mouse"}
    return {"gene": g, "lineage": x.lineage, "track": x.track, "status": x.status,
            "meta_log2fc": round(float(x.meta_log2fc), 2), "uc_papers": int(x.uc_papers),
            "criteria": json.dumps(c, ensure_ascii=False), "per_cell_log2fc": None if np.isnan(pc) else round(pc, 2),
            "crypt_share": None if np.isnan(crypt) else round(crypt, 3), "regeneration_flag": regen,
            "greens": sum(v[0] == "green" for k, v in c.items() if k != "mouse"),
            "reds": sum(v[0] == "red" for k, v in c.items() if k != "mouse"), "score": int(sum(pts.values())),
            "experiment": _experiment(x.lineage, surface, bool(strong) or (pscore or 0) >= 0.6, papers, regen)}


def _rank(rows: list[dict]) -> pd.DataFrame:
    out = pd.DataFrame(rows).sort_values(["score", "greens", "meta_log2fc"], ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def run() -> dict:
    from .evaluate import clinical_targets
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    t = _load()
    out = _rank([_row(x, t) for x in d[d.status == "first-in-class opportunity"].itertuples()])
    out.to_parquet(PROCESSED / "validation_shortlist.parquet", index=False)
    cal = _rank([_row(x, t) for x in d[d.gene.isin(set(clinical_targets()["gene"]))].itertuples()])
    cal.to_parquet(PROCESSED / "validation_shortlist_calibration.parquet", index=False)

    def greens(df):
        return {k: round(float(df.criteria.map(lambda s: json.loads(s)[k][0] == "green").mean()), 2) for k in CRITERIA[:-1]}
    s = {"genes": len(out), "top10": out.gene.head(10).tolist(), "no_red": int((out.reds == 0).sum()),
         "regeneration_flag": out.gene[out.regeneration_flag].tolist(),
         "mendelian_flag": out.gene[out.criteria.str.contains("loss causes")].tolist(),
         "green_share_first_in_class": greens(out), "green_share_clinical_targets": greens(cal),
         "clinical_targets": len(cal), "median_points_first_in_class": float(out.score.median()),
         "median_points_clinical_targets": float(cal.score.median())}
    (PROCESSED / "validation_shortlist_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, ensure_ascii=False))
    for f in ("validation_shortlist.parquet", "validation_shortlist_calibration.parquet"):
        d = pd.read_parquet(PROCESSED / f)
        print("\n", f)
        for r in d.itertuples():
            c = json.loads(r.criteria)
            print(f"{r.rank:>2} {r.gene:<11}{r.lineage:<11}{r.score:>3} pc={r.per_cell_log2fc} " +
                  " ".join(f"{k[:4]}={v[0][:2]}" for k, v in c.items() if k != "mouse"))
