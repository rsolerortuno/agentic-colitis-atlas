"""M7: deterministic, self-contained HTML report from runs/<run>/findings.json (Jinja2 + inline Plotly + Cytoscape.js
+ 3Dmol.js). Layout follows the kidney atlas; every view is rebuilt from data/processed.

The LLM only wrote the narrative paragraphs and agent findings inside findings.json; layout, tables, figures,
citations and provenance are produced here. The output opens without a server or network.
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from . import report_figures as F
from .claim_display import statement_html
from .graph import GraphStore, evidence_id
from .paths import PROCESSED, RAW, REPORTS, ROOT
from .report_help import HELP
from .report_text import NAV, NAV_TEXT, PLAIN

ASSETS = Path(__file__).parent / "report_assets"
COHORT = {"GSE193677": "MSCCR (RNA-seq)", "GSE206285": "UNIFI", "GSE92415": "PURSUIT", "GSE16879": "GSE16879",
          "GSE73661": "GSE73661", "GSE109142": "PROTECT (RNA-seq)", "GSE117993": "RISK (RNA-seq)",
          "GSE165512": "GSE165512 (RNA-seq)"}


def anchor(ev_id: str) -> str:
    return "ev-" + "".join(c if c.isalnum() else "-" for c in ev_id)


def collect_ids(f: dict, extra_ids: list[str] | None = None) -> list[str]:
    ids: list[str] = []
    for p in f["narrative"].get("summary", []):
        ids += p["evidence_ids"]
    for paras in f["narrative"].get("sections", {}).values():
        for p in paras:
            ids += p["evidence_ids"]
    for t in f.get("targets", []):
        ids += t["evidence_ids"]
    for block in f.get("agent_findings", {}).values():
        for x in block.get("findings", []):
            ids += x["evidence_ids"]
    for r in f["convergence"]["ranking"][:25]:
        for v in r["evidence_ids"].values():
            ids += v
    for rows in (f["convergence"].get("profiles") or {}).values():
        for r in rows[:25]:
            for v in r["evidence_ids"].values():
                ids += v
    return list(dict.fromkeys(ids + (extra_ids or [])))


SKIP = {"evidence_id", "type", "source", "date", "sentence", "url", "labels", "result_sha256", "interpretation"}


def evidence_view(e: dict) -> dict:
    t, link, label = e["type"], e.get("url"), e["source"]
    if t in ("literature_relation", "provenance", "pubtator_annotation") and e.get("pmid"):
        link, label = f"https://pubmed.ncbi.nlm.nih.gov/{e['pmid']}/", f"PMID {e['pmid']}"
    if t == "omics_de":
        detail = (f"{e['gene']} · {COHORT.get(e['cohort'], e['cohort'])} · {e['contrast_text']}: log2FC {e['log2fc']}, "
                  f"padj {e['padj']} (n {e['n_disease']} vs {e['n_control']})")
    elif t == "expression_meta":
        detail = (f"{e['gene']} · meta-analysis of {e['cohorts']} cohorts: log2FC {e['log2fc']}, padj {e['padj']}, "
                  f"I² {e['i2']}; significant up in {e['cohorts_significant_up']}, down in {e['cohorts_significant_down']}")
    elif t == "mendelian_randomization":
        coloc = e.get("colocalises_with_ibd_gwas")
        detail = (f"{e['gene']} · plasma protein (UKB-PPP cis-pQTL) on UC risk: OR per SD {e['odds_ratio_per_sd']}, "
                  f"p {e['p']}, q {e['q']}, {e['n_instruments']} instrument(s), {e['method']}"
                  + ("" if coloc is None else f"; colocalises with IBD GWAS: {'yes' if coloc else 'no'}")
                  + (f"; FinnGen OR {e['finngen_odds_ratio_per_sd']} (p {e['finngen_p']}), combined q {e['combined_q']}"
                     if "finngen_p" in e else ""))
    elif t == "mendelian_randomization_gut":
        coloc = e.get("colocalises_with_ibd_gwas")
        detail = (f"{e['gene']} · {e['exposure']} ({e['tissue_or_cell']}) on UC risk: OR per SD {e['odds_ratio_per_sd']}, "
                  f"p {e['p']}, q {e['q']}, {e['n_instruments']} instrument(s)"
                  + ("" if coloc is None else f"; colocalises with IBD GWAS: {'yes' if coloc else 'no'}")
                  + (f"; FinnGen OR {e['finngen_odds_ratio_per_sd']} (p {e['finngen_p']}), combined q {e['combined_q']}"
                     if "finngen_p" in e else ""))
    elif t == "spatial_visium":
        detail = (f"{e['gene']} · Visium (GSE189184): co-varies most with {e['compartment']} (rho {e['compartment_rho']}, "
                  f"first in {e['sections_compartment_first']} of {e['uc_sections']} UC sections); UC vs healthy "
                  f"log2 ratio {e['log2_ratio_uc_vs_healthy']}")
    elif t == "spatial_visium_cell":
        nbs = ", ".join(f"{x['cell_group']} (log2 {x['log2_enrichment']}, padj {x['padj']})" for x in e["enriched_neighbors"]) or "none"
        detail = (f"{e['gene']} · deconvolved Visium (GSE189184, {e['uc_spots']} UC spots): in-situ expresser "
                  f"{e['top_in_situ_expresser']} ({e['in_situ_contribution']} of the signal), composition-corrected "
                  f"source {e['composition_corrected_source']}; neighbours of {e['n_high_spots']} high spots: {nbs}")
    elif t == "protein_structure":
        detail = (f"{e['gene']} ({e['accession']}) AlphaFold {e['model']} v{e['version']}: mean pLDDT {e['mean_plddt']}, "
                  f"{e['n_cavities']} cavities, best pocket score {e['best_pocket_score']}; ligand-site recall "
                  f"{e.get('ligand_site_recall')}")
    elif t == "proteomics_de":
        detail = (f"{e['gene']} · protein in UC colon mucosa ({e['study']}): log2 UC/healthy {e['log2_ratio']}, padj "
                  f"{e['padj']}; RNA meta log2FC {e.get('rna_meta_log2fc')}; replicates the RNA change: "
                  f"{e.get('replicates_rna_change') or 'not tested'}")
    elif t == "mouse_colitis":
        detail = (f"{e['gene']} · mouse models (MGI): {e['direction']}; spontaneous colitis "
                  f"{e['spontaneous_colitis_genotypes']}, more severe induced colitis "
                  f"{e['more_severe_induced_colitis_genotypes']}, protected {e['protected_from_induced_colitis_genotypes']} "
                  f"genotype(s); PMIDs {', '.join(e.get('pmids', [])[:4])}")
    elif t == "spatial_cosmx":
        det = "; ".join(f"{x['panel']} panel: best {x['best_lineage']} {x['pct_cells_best']} % vs {x['pct_cells_other_median']} % "
                        f"of other cells ({'above background' if x['above_background'] else 'background level'})"
                        for x in e.get("detection", []))
        detail = f"{e['gene']} · CosMx single cells, inflamed UC rectum (GSE312415): {det}"
    elif t == "pubtator_relations":
        uc = ", ".join(f"{k} {v}" for k, v in e.get("with_uc_or_ibd", {}).items()) or "none"
        detail = (f"{e['gene']} · PubTator3 machine-extracted relations: {e['relations']} in {e['publications']} "
                  f"publications; with UC/IBD/colitis: {uc}")
    elif t == "literature_relation":
        detail = f"“{e['sentence']}” ({e['subject_text']} → {e['object_text']})"
    else:
        detail = "; ".join(f"{k} {v}" for k, v in e.items() if k not in SKIP and v not in (None, "", []))[:400]
    return {"id": e["evidence_id"], "anchor": anchor(e["evidence_id"]), "type": t, "label": label, "link": link,
            "detail": detail}


def network(findings: dict, store: GraphStore) -> dict:
    genes = [r["gene"] for r in findings["convergence"]["ranking"][:30]]
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    nodes, edges, seen = [], [], set()
    uc = "MONDO:0005101"
    nodes.append({"data": {"id": uc, "label": "UC", "kind": "disease"}})
    for g in genes:
        n = store.resolve(g)
        if not n:
            continue
        nodes.append({"data": {"id": n["id"], "label": g, "kind": "gene"}})
        seen.add(n["id"])
        if g in meta.index and meta.at[g, "meta_padj"] < 0.05:
            edges.append({"data": {"source": n["id"], "target": uc,
                                   "kind": "upregulated_in" if meta.at[g, "meta_log2fc"] > 0 else "downregulated_in",
                                   "ev": evidence_id("expression_meta", g)}})
    drugs_per = {}
    for nid in [x["data"]["id"] for x in nodes if x["data"]["kind"] == "gene"]:
        for c in store.claims_for(nid, limit=400):
            other = c["object"] if c["subject"] == nid else c["subject"]
            if c["predicate"] == "interacts_with" and other in seen and c["subject"] < c["object"]:
                edges.append({"data": {"source": c["subject"], "target": c["object"], "kind": "interaction",
                                       "ev": c["evidence_id"]}})
            if c["predicate"] == "targets" and c["object"] == nid and drugs_per.get(nid, 0) < 3:
                d = store.get_node(c["subject"])
                if d and d["id"] not in seen:
                    nodes.append({"data": {"id": d["id"], "label": d["name"].title(), "kind": "drug"}})
                    seen.add(d["id"])
                drugs_per[nid] = drugs_per.get(nid, 0) + 1
                edges.append({"data": {"source": c["subject"], "target": nid, "kind": "drug", "ev": c["evidence_id"]}})
    keep, out = set(), []
    for e in edges:
        k = (e["data"]["source"], e["data"]["target"], e["data"]["kind"])
        if k not in keep:
            keep.add(k)
            out.append({"data": {**e["data"], "id": f"e{len(out)}"}})
    return {"nodes": nodes, "edges": out}


def _num(x, nd=2):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return round(float(x), nd)


# ------------------------------------------------------------------ tables
def cell_tables(genes: list[str], store: GraphStore) -> dict:
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    inf = sh[sh.health == "inflamed"]
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    de = pd.read_parquet(PROCESSED / "sc_celltype_de.parquet")
    rows = []
    for g in genes:
        d = inf[inf.gene == g]
        if d.empty:
            continue
        top = d.sort_values("mean_cp10k", ascending=False).iloc[0]
        fam = d.groupby("compartment")["share"].sum().sort_values(ascending=False)
        v = d.mean_cp10k.to_numpy(float)
        tau = float((1 - v / v.max()).sum() / (len(v) - 1)) if v.max() > 0 and len(v) > 1 else None
        x = de[(de.gene == g)].sort_values("padj")
        ivh = None if x.empty else {"log2fc": _num(x.iloc[0].log2fc), "padj": float(f"{x.iloc[0].padj:.2g}"),
                                    "cluster": x.iloc[0].cluster}
        rows.append({"gene": g, "top_label": top.cluster, "top_group": top.compartment,
                     "top_pct": _num(top.pct_expressing, 1), "top_share": f"{fam.iloc[0]:.0%} {fam.index[0]}",
                     "neutrophil_share": _num(fam.get("Neutrophil", 0.0)), "tau": _num(tau),
                     "retained": _num(comp.at[g, "retained"]) if g in comp.index else None,
                     "verdict": comp.at[g, "verdict"] if g in comp.index else "",
                     "ivh": ivh, "evidence_ids": [evidence_id("sc_summary", g), evidence_id("composition_effect", g)]})
    rep = pd.read_parquet(PROCESSED / "sc_replication.parquet").set_index("gene")
    sc_rep_rows = [{"gene": g, "lineage_2019": rep.at[g, "lineage_2019"], "share_2019": _num(rep.at[g, "share_2019"]),
                    "lineage_2023": rep.at[g, "lineage_2023"], "share_2023": _num(rep.at[g, "share_2023"])}
                   for g in genes if g in rep.index][:20]
    vis = json.loads((PROCESSED / "spatial_visium_summary.json").read_text(encoding="utf-8"))
    tv = {t["gene"]: t for t in vis["targets"]}
    comp_v = pd.read_parquet(PROCESSED / "spatial_visium.parquet").set_index("gene")
    spatial_rows = [{"gene": g, "top_expresser": (tv[g]["top_in_situ_expresser"] or "").replace("_", " "),
                     "contribution": tv[g]["in_situ_contribution"],
                     "corrected": (tv[g]["top_corrected_source"] or "").replace("_", " "),
                     "compartment": comp_v.at[g, "compartment"] if g in comp_v.index else None,
                     "n_high": tv[g]["n_high_spots"], "neighbors": tv[g]["enriched_neighbors"][:4],
                     "evidence_ids": [evidence_id("spatial_visium_cell", g), evidence_id("spatial_visium", g)]}
                    for g in genes if g in tv]
    att = pd.read_parquet(PROCESSED / "spatial_attribution.parquet")
    att = att[(att.state == "inflamed") & att.gene.isin(genes)]
    xen = att.groupby(["gene", "lineage"]).agg(share=("share", "mean"), pct=("pct_expressing", "mean")).reset_index()
    xen_rows = [{"gene": g, "lineage": d.sort_values("share").iloc[-1].lineage, "share": _num(d.share.max()),
                 "pct": _num(d.sort_values("share").iloc[-1].pct, 1)} for g, d in xen.groupby("gene")]
    sec = pd.read_parquet(PROCESSED / "spatial_sections.parquet")
    cols = [c for c in ("gsm", "study", "state", "cells", "median_counts_per_cell", "assigned_fraction") if c in sec]
    cells3 = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet", columns=["health", "cluster", "n_cells"]).drop_duplicates()
    return {"cell_rows": rows, "sc_rep_rows": sc_rep_rows, "spatial_rows": spatial_rows, "xenium_rows": xen_rows,
            "xenium_sections": sec[cols].round(2).to_dict("records"), "visium": vis,
            "sc3_cells": int(cells3.n_cells.sum())}


def genetics_tables(genes: list[str]) -> dict:
    from .evaluate import clinical_targets
    clin = set(clinical_targets()["gene"])
    p = pd.read_parquet(PROCESSED / "mr_pqtl.parquet")
    c = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    pi = p.set_index("gene")
    sig = p[p.q < 0.05]
    out = {"plasma": {"tested": len(p), "significant": len(sig),
                      "colocalising": int(sig.colocalises.fillna(False).astype(bool).sum()),
                      "clinical_tested": int(p.gene.isin(clin).sum()),
                      "clinical_significant": sorted(set(sig.gene) & clin)}}
    for layer, key in (("colon_tissue", "colon"), ("gut_cell", "cell")):
        d = c[c.layer == layer]
        s = d[d.q < 0.05]
        out[key] = {"tested": len(d), "significant": len(s), "colocalising": int(s.colocalises.fillna(False).astype(bool).sum())}
    col = c[(c.layer == "colon_tissue") & (c.q < 0.05)].sort_values("p").head(25)

    def coloc(x):
        return "yes" if x is True else ("no" if x is False else "-")
    def fg(x):
        if pd.isna(getattr(x, "p_finngen", np.nan)):
            return "—"
        return f"OR {np.exp(x.beta_finngen):.2f}, p {x.p_finngen:.1e}" + (" ✓" if x.replicates_finngen is True else "")
    for key, d in (("plasma", p), ("colon", c[c.layer == "colon_tissue"]), ("cell", c[c.layer == "gut_cell"])):
        s_ = d[d.q < 0.05]
        out[key]["finngen_replicated"] = int((s_.get("replicates_finngen", pd.Series(dtype=object)) == True).sum())  # noqa: E712
        out[key]["combined_significant"] = int((d.get("q_combined", pd.Series(dtype=float)) < 0.05).sum())
    out["plasma_finngen"] = {x.gene: fg(x) for x in p.itertuples()}
    out["colon_rows"] = [{"gene": x.gene, "label": x.label, "or": _num(x.or_per_sd), "q": f"{x.q:.1e}", "finngen": fg(x),
                          "n": int(x.n_instruments), "coloc": coloc(x.colocalises),
                          "plasma": f"OR {pi.at[x.gene, 'or_per_sd']:.2f}, q {pi.at[x.gene, 'q']:.1e}" if x.gene in pi.index else None}
                         for x in col.itertuples()]
    ci = c.set_index(["gene", "layer"])

    def fmt(g, layer):
        if (g, layer) not in ci.index:
            return "not tested"
        r = ci.loc[(g, layer)]
        lab = f" in {r['label'].replace('_', ' ')}" if layer == "gut_cell" else ""
        return f"OR {r['or_per_sd']:.2f}{lab}, q {r['q']:.1e}" + (" (colocalises)" if r["colocalises"] is True else "")
    out["candidates"] = [{"gene": g,
                          "plasma": (f"OR {pi.at[g, 'or_per_sd']:.2f}, q {pi.at[g, 'q']:.1e}" + (" (colocalises)" if pi.at[g, "colocalises"] is True else ""))
                          if g in pi.index else "not tested",
                          "colon": fmt(g, "colon_tissue"), "cell": fmt(g, "gut_cell")} for g in genes[:20]]
    return out


def structure_tables(report_genes: list[str]) -> dict:
    p = PROCESSED / "structure_summary.json"
    if not p.exists():
        return {"rows": [], "controls": [], "config": {"weights": {}, "pocket_volume": [0, 0], "controls": [], "min_plddt_keep": 50},
                "modeled": 0}
    st = json.loads(p.read_text(encoding="utf-8"))
    rows = []
    for r in st["results"]:
        if not r["modeled"] or not r["pockets"]:
            continue
        b = r["pockets"][0]
        v = r["ligand_site_validation"]
        rows.append({"gene": r["gene"], "control": r["control"], "uniprot": r["uniprot"], "mean_plddt": _num(r["mean_plddt"], 1),
                     "rank": b["rank"], "volume_a3": b["volume_a3"], "max_depth_a": b["max_depth_a"],
                     "avg_hydropathy": b["avg_hydropathy"], "pocket_score": b["pocket_score"],
                     "recall": v.get("recall"), "recall_top": v.get("recall_in_top_pockets"),
                     "evidence_ids": [evidence_id("protein_structure", r["gene"])]})
    rows.sort(key=lambda r: -r["pocket_score"])
    controls = [{"gene": r["gene"], "rank": i + 1, "of": len(rows)} for i, r in enumerate(rows) if r["control"]]
    return {"rows": rows, "controls": controls, "config": st["config"], "modeled": st["modeled"], "summary": st}


def discovery_tables(layers: dict) -> dict:
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    st = {r["gene"]: r for r in json.loads((PROCESSED / "structure_summary.json").read_text())["results"]} \
        if (PROCESSED / "structure_summary.json").exists() else {}
    names = {}
    ext = json.loads((PROCESSED / "external_targets.json").read_text(encoding="utf-8"))
    for g, e in ext.items():
        names[g] = ((e.get("open_targets") or {}).get("name") or "")
    fic = d[d.status == "first-in-class opportunity"].sort_values("novelty_score", ascending=False)
    cards = []
    for r in fic.groupby("track").head(4).sort_values("novelty_score", ascending=False).itertuples():
        why = (f"Up {r.meta_log2fc:.1f} log2 units in UC in {int(r.cohorts_up)} cohorts; {r.retained:.0%} of the rise survives "
               f"the composition adjustment; made mostly by {r.lineage} cells; {int(r.uc_papers)} UC papers; tractable as "
               f"{', '.join(r.tractable)}; at most {r.vital_max_ntpm:.0f} nTPM in a vital organ.")
        pk = st.get(r.gene, {}).get("best_pocket_score")
        cards.append({"gene": r.gene, "name": names.get(r.gene), "track": r.track, "rank_in_track": int(r.rank_in_track),
                      "novelty": round(float(r.novelty_score), 2), "papers": int(r.uc_papers), "retained": float(r.retained),
                      "pocket": pk, "why": why, "evidence_ids": [evidence_id("discovery_profile", r.gene)]})
    top = d.sort_values("novelty_score", ascending=False).groupby("track").head(40)
    rows = [{"gene": r.gene, "track": r.track, "lineage": r.lineage, "neutrophil_share": float(r.neutrophil_share or 0),
             "novelty": round(float(r.novelty_score), 2), "log2fc": round(float(r.meta_log2fc), 2),
             "retained": _num(r.retained), "papers": int(r.uc_papers), "clinical": bool(r.clinical_drug),
             "tractable": ", ".join(r.tractable), "vital": _num(r.vital_max_ntpm, 1), "status": r.status,
             "evidence_ids": [evidence_id("discovery_profile", r.gene)]} for r in top.itertuples()]
    return {"pool": len(d), "by_track": d.track.value_counts().to_dict(), "fic_by_track": fic.track.value_counts().to_dict(),
            "status": d.status.value_counts().to_dict(), "cards": cards, "rows": rows, "weights": layers["novelty"]["weights"]}


def treatment_tables() -> dict:
    r = pd.read_parquet(PROCESSED / "response_genes.parquet")
    n = r.groupby("arm")[["n_nonresponders", "n_responders"]].first().sum(axis=1)
    top = {}
    for arm in ("ustekinumab", "placebo (UNIFI)", "infliximab (ACT1)", "placebo (ACT1)", "etrolizumab", "golimumab", "vedolizumab"):
        d = r[r.arm == arm].sort_values("p").head(8)
        top[arm] = [{"gene": x.gene, "t": _num(x.t), "q": f"{x.padj:.1e}",
                     "direction": "higher in non-responders" if x.t > 0 else "higher in responders"} for x in d.itertuples()]
    it = pd.read_parquet(PROCESSED / "response_interaction.parquet").sort_values("p_interaction").head(10)
    return {"n": {"ustekinumab": int(n.get("ustekinumab", 0)), "placebo_unifi": int(n.get("placebo (UNIFI)", 0)),
                  "infliximab_act1": int(n.get("infliximab (ACT1)", 0)), "placebo_act1": int(n.get("placebo (ACT1)", 0)),
                  "etrolizumab": int(n.get("etrolizumab", 0))},
            "significant": {a: int((d.padj < 0.05).sum()) for a, d in r.groupby("arm")},
            "interaction_significant": int((pd.read_parquet(PROCESSED / "response_interaction.parquet").padj_interaction < 0.05).sum()),
            "top": top, "check_genes": F.CHECK_GENES[:5],
            "interaction": [{"gene": x.gene, "drug": x.drug, "lor": _num(x.log_or_interaction), "q": f"{x.padj_interaction:.2f}"}
                            for x in it.itertuples()]}


def confirm_rows(genes: list[str]) -> list[dict]:
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    de = de[de.contrast.isin(["inflamed_vs_control", "uc_vs_control"])]
    w = de.pivot_table(index="gene", columns="cohort", values="log2fc")
    q = de.pivot_table(index="gene", columns="cohort", values="padj")
    out = []
    for g in genes[:20]:
        if g not in w.index:
            continue
        m = w.at[g, "GSE193677"]
        rep = [c for c in ("GSE109142", "GSE117993") if pd.notna(q.at[g, c]) and q.at[g, c] < 0.05 and np.sign(w.at[g, c]) == np.sign(m)]
        out.append({"gene": g, "msccr": _num(m), "protect": _num(w.at[g, "GSE109142"]), "protect_q": f"{q.at[g, 'GSE109142']:.1e}" if pd.notna(q.at[g, "GSE109142"]) else "—",
                    "risk": _num(w.at[g, "GSE117993"]), "risk_q": f"{q.at[g, 'GSE117993']:.1e}" if pd.notna(q.at[g, "GSE117993"]) else "—",
                    "replicates": "both" if len(rep) == 2 else ("one" if rep else "no")})
    return out


def agents_vs_ranking(f: dict) -> list[dict]:
    t = pd.read_parquet(PROCESSED / "convergence_table.parquet").set_index("gene")
    rank = {r["gene"]: r["rank"] for r in f["convergence"]["ranking"]}
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    mr = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").set_index("gene")
    out = []
    for i, x in enumerate(f.get("targets", []), 1):
        g = x["gene"]
        out.append({"gene": g, "category": x["category"].replace("_", " "), "position": i, "rank": rank.get(g),
                    "log2fc": _num(t.at[g, "meta_log2fc"]) if g in t.index else None,
                    "retained": _num(comp.at[g, "retained"]) if g in comp.index else None,
                    "mr_q": float(f"{mr.at[g, 'q']:.2g}") if g in mr.index else None,
                    "target": _num(t.at[g, "score_target"]) if g in t.index else None,
                    "biomarker": _num(t.at[g, "score_biomarker"]) if g in t.index else None})
    return out


def datasets() -> list[dict]:
    rows = []
    for mp in [RAW / "manifest.json"] + sorted((RAW / "spatial").glob("*/manifest.json")):
        m = json.loads(mp.read_text(encoding="utf-8"))
        for k, v in sorted(m.items()):
            url = v.get("url") or ""
            rows.append({"file": k if mp.parent == RAW else f"{mp.parent.name}/{k}", "url": url if url.startswith("http") else None,
                         "host": url.split("/")[2] if url.startswith("http") else url[:60], "bytes": v.get("bytes"),
                         "retrieved": (v.get("retrieved_utc") or "")[:10], "sha256": v.get("sha256")})
    return rows


def badges(f: dict, cell: dict) -> dict:
    stage = ["PHASE_1", "PHASE_2", "PHASE_3", "PHASE_4", "APPROVAL"]
    mr = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").set_index("gene")
    out: dict[str, dict] = {}
    for r in cell["cell_rows"]:
        out.setdefault(r["gene"], {}).update(cell_label=r["top_group"], cell_share=None, verdict=r["verdict"])
    for r in f["convergence"]["ranking"]:
        d = out.setdefault(r["gene"], {})
        d["rank"] = r["rank"]
        known = [x for x in r["drugs"] if x["max_stage"] in stage]
        d["drug"] = max(known, key=lambda x: stage.index(x["max_stage"]))["name"].title() if known else None
    for g, d in out.items():
        if g in mr.index and mr.at[g, "q"] < 0.05:
            d["mr"] = f"MR OR {mr.at[g, 'or_per_sd']:.2f}"
    return out


# ------------------------------------------------------------------ render
def render(findings_path: Path, out: Path | None = None) -> Path:
    from . import report_extra as X
    from .protein_viewer import build as build_protein_viewer, gzip_b64, script_json
    f = json.loads(findings_path.read_text(encoding="utf-8"))
    f["models"] = {k: str(v).replace("�", "·") for k, v in f["models"].items()}
    store = GraphStore()
    ranking = f["convergence"]["ranking"]
    top = [r["gene"] for r in ranking[:25]]
    targets = [t["gene"] for t in f.get("targets", [])]
    report_genes = list(dict.fromkeys(targets + top))
    hl = set(targets + top[:10])
    layers = json.loads((PROCESSED / "report_layers.json").read_text(encoding="utf-8"))
    sens = PROCESSED / "sensitivity_convergence.parquet"
    if sens.exists():
        s = pd.read_parquet(sens).set_index("gene")
        for r in ranking[:25]:
            if r["gene"] in s.index:
                x = s.loc[r["gene"]]
                r["rank_interval_90"] = {"rank_p05": int(x.rank_p05), "rank_p95": int(x.rank_p95), "p_top10": round(float(x.p_top10), 3)}
    sw = f["convergence"]["config"]["stream_weights"]
    cell = cell_tables(report_genes[:30], store)
    G = genetics_tables(report_genes)
    S = structure_tables(report_genes)
    D = discovery_tables(layers)
    T = treatment_tables()
    disc_genes = [c["gene"] for c in D["cards"]]
    from .spatial_deconv import CELL_GROUPS  # noqa: F401
    cohort_html, _, _ = F.cohort_scatter("GSE193677", "GSE206285", hl)
    figs = {
        "stream_stack": F.stream_stack(ranking, sw), "component_heatmap": F.component_heatmap(ranking, list(sw)),
        "volcano_ic": F.volcano("inflamed_vs_control", hl, "Inflamed UC vs control colon (MSCCR, RNA-seq)",
                                "log2 fold change (UC − control)"),
        "volcano_paired": F.volcano("inflamed_vs_noninflamed_paired", hl, "Inflamed vs non-inflamed mucosa, same patients (MSCCR)",
                                    "log2 fold change (inflamed − non-inflamed)"),
        "cohort_agree": cohort_html, "replication": F.replication_heatmap(report_genes[:30]),
        "activity": F.activity_bars(report_genes[:20]),
        "dotplot": F.dotplot(report_genes[:28] + [g for g in disc_genes if g not in report_genes][:8]),
        "share_stack": F.share_stack(report_genes[:30]), "neutrophil_retained": F.neutrophil_vs_retained(report_genes[:25]),
        "map_uc": F.visium_map("GSM5695840", "UC resection (C2): estimated cell shares and genes"),
        "map_healthy": F.visium_map("GSM5695843", "Healthy colon (C5): estimated cell shares and genes"),
        "expressers": F.expresser_stack(report_genes[:30]), "neighbors": F.neighbor_heatmap(report_genes[:30]),
        "pathways": F.pathway_heatmap(), "tf_bars": F.tf_bars(),
        "mr": F.mr_volcano("mr_pqtl.parquet", "Plasma protein (UKB-PPP cis-pQTL) on UC risk"),
        "mr_colon": F.mr_volcano("mr_colon.parquet", "Colon expression (GTEx v10, CEDAR eQTL) on UC risk", layer="colon_tissue"),
        "mr_agree": F.genetics_agreement(), "mr_cells": F.mr_celltypes(report_genes),
        "pocket_bubble": F.pocket_bubble(report_genes, disc_genes) if S["rows"] else None,
        "pocket_components": F.pocket_components(S["config"]["weights"]) if S["rows"] else None,
        "ligand_validation": F.ligand_validation(S["summary"]) if S["rows"] else None,
        "control_ranking": F.control_ranking() if S["rows"] else None,
        "discovery_scatter": F.discovery_scatter(), "discovery_bars": F.discovery_bars(D["weights"]),
        "safety": F.safety_scatter(),
        "treatment_controls": F.treatment_controls(), "treatment_arms": F.treatment_arms(hl),
        "response": F.response_heatmap(), "subgroup_counts": F.subgroup_counts(layers["subgroups"]),
        "activity_axis": F.activity_axis(), "biologics": F.biologics_figure(), "immune_classes": F.immune_classes(layers["immune"]),
        "composition_shift": F.composition_shift(report_genes[:22]),
        "confirm_protect": F.cohort_scatter("GSE193677", "GSE109142", hl)[0],
        "confirm_risk": F.cohort_scatter("GSE193677", "GSE117993", hl)[0],
        "rank_conv": F.rank_intervals("convergence"), "rank_target": F.rank_intervals("target"),
        "rank_biomarker": F.rank_intervals("biomarker"),
        "drug_stage": F.drug_stage(ranking),
    }
    figs["sc_replication"], cell["sc_rep"] = F.sc_replication_figure()
    figs.update({"proteomics": F.proteomics_scatter(hl), "finngen": F.finngen_scatter(), "pubtator": F.pubtator_chart()})
    third = {}
    if (PROCESSED / "proteomics_summary.json").exists():
        pr = pd.read_parquet(PROCESSED / "proteomics.parquet").set_index("gene")
        third["proteomics"] = json.loads((PROCESSED / "proteomics_summary.json").read_text())
        third["proteomics_rows"] = [{"gene": g, "rna": _num(pr.at[g, "meta_log2fc"]), "protein": _num(pr.at[g, "log2_ratio"]),
                                     "padj": f"{pr.at[g, 'padj']:.1e}" if pd.notna(pr.at[g, "padj"]) else "—",
                                     "replicates": pr.at[g, "replicates"] or "not tested",
                                     "evidence_ids": [evidence_id("proteomics_de", g)]}
                                    for g in report_genes[:30] if g in pr.index]
    if (PROCESSED / "mouse_colitis_summary.json").exists():
        ms = json.loads((PROCESSED / "mouse_colitis_summary.json").read_text())
        mc = pd.read_parquet(PROCESSED / "mouse_colitis.parquet").set_index("gene")
        third["mouse"] = ms
        figs["mouse"] = F.mouse_rates(ms)
        third["mouse_rows"] = [{"gene": g, "direction": mc.at[g, "direction"], "spont": int(mc.at[g, "spontaneous_genotypes"]),
                                "worse": int(mc.at[g, "worse_induced_genotypes"]), "protected": int(mc.at[g, "protected_genotypes"]),
                                "evidence_ids": [evidence_id("mouse_colitis", g)]}
                               for g in list(dict.fromkeys(report_genes + sorted(ms["clinical_targets_with_phenotype"]))) if g in mc.index]
        third["mouse_candidates_without"] = [g for g in report_genes[:25] if g not in mc.index]
    if (PROCESSED / "cosmx_detection.parquet").exists():
        cd = pd.read_parquet(PROCESSED / "cosmx_detection.parquet")
        cs = pd.read_parquet(PROCESSED / "cosmx_sections.parquet")
        third["cosmx_sections"] = cs.round(3).to_dict("records")
        cd6 = cd[cd.panel == "6k"].set_index("gene")
        third["cosmx_rows"] = [{"gene": g, "best": cd6.at[g, "best_lineage"], "pct_best": cd6.at[g, "pct_best"],
                                "pct_other": cd6.at[g, "pct_others_median"], "spec": cd6.at[g, "specificity"],
                                "above": bool(cd6.at[g, "above_background"]), "evidence_ids": [evidence_id("spatial_cosmx", g)]}
                               for g in report_genes[:30] if g in cd6.index]
        third["cosmx_detected"] = cd6[cd6.above_background].sort_values("specificity", ascending=False).index.tolist()[:15]
        third["cosmx_counts"] = {"tested": int(len(cd6)), "above": int(cd6.above_background.sum())}
    lit_rel = PROCESSED.parents[0] / "literature" / "pubtator_relations.parquet"
    if lit_rel.exists():
        pubs = pd.read_parquet(PROCESSED.parents[0] / "literature" / "publications.parquet")
        third["literature"] = {"publications": len(pubs), "pubtator_relations": len(pd.read_parquet(lit_rel))}
    lit1, lit2 = F.literature_charts()
    figs.update({"lit_relations": lit1, "lit_years": lit2})
    evals = json.loads((ROOT / "evals" / "results.json").read_text(encoding="utf-8"))
    for k in ("clinical_targets", "temporal_wide"):
        for v in evals.get(k, {}).values():
            v.setdefault("median_rank", None)
    figs["outcomes"] = F.outcome_auroc(evals)
    extra = X.tables(report_genes)
    from .report_layers import candidate_evidence
    cards = candidate_evidence(list(dict.fromkeys(targets + top[:5])))
    bio = pd.read_parquet(PROCESSED / "biologics.parquet")
    bio_top = pd.concat([bio[bio.in_universe].head(15), bio[bio.control]]).drop_duplicates("gene")
    imm = pd.read_parquet(PROCESSED / "immune_targets.parquet")
    imm_top = pd.concat([imm.head(15), imm[imm.control]]).drop_duplicates("gene")
    rep = pd.read_parquet(PROCESSED / "replication.parquet").to_dict("records")
    ab = json.loads((ROOT / "evals" / "agent_benchmark.json").read_text(encoding="utf-8"))
    team_mean = {"cost": float(np.mean([r["cost_usd"] for r in ab["team"]])), "minutes": float(np.mean([r["minutes"] for r in ab["team"]]))}
    vpath = findings_path.parent / "verification_final.json"
    verification = json.loads(vpath.read_text(encoding="utf-8")) if vpath.exists() else None
    sess_p = findings_path.parent / "session.json"
    sess = json.loads(sess_p.read_text(encoding="utf-8")) if sess_p.exists() else {}
    cov_p = findings_path.parent / "coverage.json"
    cov = json.loads(cov_p.read_text(encoding="utf-8")) if cov_p.exists() else {}
    session = {"run": findings_path.parent.name, "statements": verification["statements"] if verification else "—",
               "unsupported": verification["unsupported"] if verification else "—",
               "cost": f"{sess.get('total_cost_usd', 0):.2f}", "minutes": round(sess.get("duration_ms", 0) / 60000, 1),
               "coverage": f"{cov.get('key_findings', '—')} key findings were all cited"}
    sb_p = ROOT / "evals" / "verifier_sabotage.json"
    sabotage = "not run"
    if sb_p.exists():
        sb = json.loads(sb_p.read_text(encoding="utf-8"))["corruptions"]
        sabotage = "; ".join(f"{k} {v['caught']}/{v['applied']}" for k, v in sb.items())
    graph_stats = json.loads((PROCESSED / "graph_stats.json").read_text())
    samples = pd.read_parquet(PROCESSED / "samples.parquet")
    xen_cells = int(pd.read_parquet(PROCESSED / "spatial_sections.parquet")["cells"].sum())
    kpis = {"universe": f["convergence"]["universe_size"], "claims": graph_stats["claims_total"], "samples": len(samples),
            "cells": cell["sc3_cells"] + int(pd.read_parquet(PROCESSED / "sc_cells.parquet")["n_cells"].sum()),
            "spots": cell["visium"]["spots"], "xenium_cells": xen_cells,
            "mr_genes": len(set(pd.read_parquet(PROCESSED / "mr_pqtl.parquet").gene) | set(pd.read_parquet(PROCESSED / "mr_colon.parquet").gene)),
            "proteins": S["modeled"], "fic": D["status"].get("first-in-class opportunity", 0), "pool": D["pool"]}
    # table evidence must be numbered too
    extra_ids = []
    blocks = [cell["cell_rows"], cell["spatial_rows"], S["rows"], D["cards"], D["rows"],
              third.get("proteomics_rows", []), third.get("mouse_rows", []), third.get("cosmx_rows", [])]
    for block in blocks:
        for r in block:
            extra_ids += r["evidence_ids"]
    present = {e["evidence_id"] for e in store.get_evidence(list(dict.fromkeys(extra_ids)))}
    for block in blocks:
        for r in block:   # a table row cites only records that exist (not every gene has every layer)
            r["evidence_ids"] = [i for i in r["evidence_ids"] if i in present]
    extra_ids = [i for i in extra_ids if i in present]
    ids = collect_ids(f, extra_ids)
    records = {e["evidence_id"]: e for e in store.get_evidence(ids)}
    ev_views = [evidence_view(records[i]) for i in ids if i in records]
    number = {v["id"]: k + 1 for k, v in enumerate(ev_views)}
    env = Environment(loader=FileSystemLoader(ASSETS), autoescape=select_autoescape(["html", "j2"]))

    def cite(ev_ids):
        parts = [f'<a class="cite" href="#{anchor(i)}" title="{html.escape(i)}">{number[i]}</a>' for i in ev_ids if i in number]
        return "<sup>[" + ",".join(parts) + "]</sup>" if parts else ""

    env.globals["H"] = HELP
    env.globals["PLAIN"] = PLAIN
    env.filters["cite"] = cite
    env.filters["statement_html"] = lambda st: statement_html(st, records)
    env.filters["evattr"] = lambda ev_ids: " ".join(ev_ids)
    env.filters["fmt"] = lambda x: "—" if x is None or (isinstance(x, float) and math.isnan(x)) else x
    viewer_p = PROCESSED / "spatial_viewer.json"
    html_out = env.get_template("report.html.j2").render(
        f=f, top25=ranking[:25], sw=sw, evidence=ev_views, figs=figs, nav=NAV, nav_text=NAV_TEXT, kpis=kpis,
        badges=badges(f, cell), cell=cell, third=third,
        score_sens=layers.get("score_sensitivity", {}), G=G, S=S, D=D, T=T, layers=layers, extra=extra, replication=rep,
        cohort_names={**COHORT, "GSE109142": "PROTECT", "GSE117993": "RISK"},
        activity_caption=F.activity_summary(report_genes[:20]), candidate_evidence=cards,
        bio_top=bio_top.to_dict("records"), imm_top=imm_top.to_dict("records"),
        confirm_rows=confirm_rows(report_genes), profiles=f["convergence"].get("profiles") or {},
        profile_weights=f["convergence"].get("profile_weights") or {}, evals=evals,
        baselines=json.loads((ROOT / "evals" / "baselines.json").read_text(encoding="utf-8")),
        agents_vs_ranking=agents_vs_ranking(f), AB=ab, team_mean=team_mean, session=session, sabotage=sabotage,
        verification=verification, graph_stats=graph_stats, datasets=datasets(),
        css=(ASSETS / "style.css").read_text(encoding="utf-8"), plotly_js=get_plotlyjs(),
        cytoscape_js=(ASSETS / "cytoscape.min.js").read_text(encoding="utf-8"),
        protein_library_js=(ASSETS / "3Dmol-min.js").read_text(encoding="utf-8"),
        protein_viewer_js=(ASSETS / "protein_viewer.js").read_text(encoding="utf-8"),
        protein_b64=gzip_b64(script_json(build_protein_viewer())),
        viewer_b64=gzip_b64(viewer_p.read_text(encoding="utf-8") if viewer_p.exists() else "null"),
        network_json=json.dumps(network(f, store)))
    out = out or REPORTS / "uc_report.html"
    out.write_text(html_out, encoding="utf-8")
    return out
