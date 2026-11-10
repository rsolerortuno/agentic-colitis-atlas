"""M7: deterministic, self-contained HTML report from runs/<run>/findings.json (Jinja2 + inline Plotly + Cytoscape.js).

The LLM only wrote the narrative paragraphs inside findings.json; layout, tables, figures, citations and provenance
are produced here. The output opens without a server or network.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from .claim_display import statement_html
from .graph import GraphStore
from .paths import PROCESSED, REPORTS

ASSETS = Path(__file__).parent / "report_assets"
COHORT = {"GSE193677": "MSCCR (RNA-seq)", "GSE206285": "UNIFI", "GSE92415": "PURSUIT", "GSE16879": "GSE16879",
          "GSE73661": "GSE73661", "GSE109142": "PROTECT (RNA-seq)", "GSE117993": "RISK (RNA-seq)",
          "GSE165512": "GSE165512 (RNA-seq)"}
NAV = [("Start", [("Overview", "overview", "The question, the answer in brief and the prioritised candidates.")]),
       ("Findings", [("Candidates", "candidates", "The convergence ranking, target and biomarker profiles."),
                     ("Discovery", "discovery", "New, untargeted candidates that pass composition, novelty and safety filters.")]),
       ("Evidence", [("Expression", "expression", "What changes in inflamed colon, how it replicates across eight cohorts."),
                     ("Cells & space", "cells", "Which cell makes each gene, neutrophils, composition and the tissue view."),
                     ("Pathways", "pathways", "Immune pathways, transcription factors and cell-cell circuits."),
                     ("Genetics", "genetics", "Does genetically higher plasma protein change UC risk? (Mendelian randomisation)")]),
       ("Clinical use", [("Treatment response", "treatment", "What separates responders to four biologics, and placebo."),
                         ("Literature & drugs", "literature", "What the papers and the drug databases say."),
                         ("Network", "network", "Candidate genes, drugs and interactions as a graph.")]),
       ("Trust and methods", [("Validation", "validation", "Recovery of known UC targets and known-answer questions."),
                              ("Methods & evidence", "methods", "Datasets, verification, limits and every cited record.")])]


def anchor(ev_id: str) -> str:
    return "ev-" + "".join(c if c.isalnum() else "-" for c in ev_id)


def collect_ids(f: dict) -> list[str]:
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
                ids += v[:6]
    return list(dict.fromkeys(ids))


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
                  + ("" if coloc is None else f"; colocalises with IBD GWAS: {'yes' if coloc else 'no'}"))
    elif t == "spatial_visium":
        detail = (f"{e['gene']} · Visium (GSE189184): co-varies most with {e['compartment']} (rho {e['compartment_rho']}, "
                  f"first in {e['sections_compartment_first']} of {e['uc_sections']} UC sections); UC vs healthy "
                  f"log2 ratio {e['log2_ratio_uc_vs_healthy']}")
    elif t == "literature_relation":
        detail = f"“{e['sentence']}” ({e['subject_text']} → {e['object_text']})"
    else:
        detail = "; ".join(f"{k} {v}" for k, v in e.items() if k not in SKIP and v not in (None, "", []))[:400]
    return {"id": e["evidence_id"], "anchor": anchor(e["evidence_id"]), "type": t, "label": label, "link": link,
            "detail": detail}


# ------------------------------------------------------------------ figures
def _fig(fig: go.Figure, h=420) -> str:
    fig.update_layout(height=h, margin=dict(l=50, r=20, t=40, b=45), paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", font=dict(size=12, color="#4d4c49"))
    fig.update_xaxes(gridcolor="rgba(137,135,129,.22)", zeroline=False)
    fig.update_yaxes(gridcolor="rgba(137,135,129,.22)", zeroline=False)
    return fig.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False, "responsive": True})


def volcano(de: pd.DataFrame, label: list[str]) -> str:
    d = de[(de.cohort == "GSE193677") & (de.contrast == "inflamed_vs_control")].dropna(subset=["padj"])
    y = -np.log10(d["padj"].clip(lower=1e-300))
    sig = (d["padj"] < 0.05) & (d["log2fc"].abs() >= 1)
    col = np.where(~sig, "#c9c8c0", np.where(d["log2fc"] > 0, "#e34948", "#2a78d6"))
    fig = go.Figure(go.Scattergl(x=d["log2fc"], y=y, mode="markers", marker=dict(size=4, color=col),
                                 text=d["gene"], hovertemplate="%{text}<br>log2FC %{x:.2f}<br>-log10 padj %{y:.1f}<extra></extra>"))
    lab = d[d["gene"].isin(label)]
    fig.add_trace(go.Scatter(x=lab["log2fc"], y=-np.log10(lab["padj"].clip(lower=1e-300)), mode="markers+text",
                             text=lab["gene"], textposition="top center", marker=dict(size=7, color="#1c1c1b"),
                             hoverinfo="skip"))
    fig.update_layout(showlegend=False, title="Inflamed UC vs control colon (MSCCR, RNA-seq)",
                      xaxis_title="log2 fold change", yaxis_title="-log10 adjusted p")
    return _fig(fig, 460)


def replication_heatmap(de: pd.DataFrame, genes: list[str]) -> str:
    d = de[de.contrast.isin(["inflamed_vs_control", "uc_vs_control"]) & de.gene.isin(genes)]
    m = d.pivot_table(index="gene", columns="cohort", values="log2fc").reindex(index=genes[::-1], columns=list(COHORT))
    fig = go.Figure(go.Heatmap(z=m.values, x=[COHORT[c] for c in m.columns], y=m.index, colorscale="RdBu_r", zmid=0,
                               colorbar=dict(title="log2FC"),
                               hovertemplate="%{y} · %{x}<br>log2FC %{z:.2f}<extra></extra>"))
    fig.update_layout(title="Top candidates across the eight cohorts (UC vs control)")
    return _fig(fig, max(380, 18 * len(genes) + 120))


def score_bars(ranking: list[dict], weights: dict) -> str:
    rows = ranking[:25][::-1]
    fig = go.Figure()
    colors = ["#2a78d6", "#1baf7a", "#eb6834", "#6250d6", "#eda100", "#e87ba4"]
    for k, c in zip(weights, colors):
        fig.add_trace(go.Bar(y=[r["gene"] for r in rows], x=[weights[k] * r["components"][k] for r in rows],
                             name=k, orientation="h", marker_color=c))
    fig.update_layout(barmode="stack", title="Convergence score by evidence stream (top 25)", xaxis_title="weighted score",
                      legend=dict(orientation="h", y=-0.12))
    return _fig(fig, 640)


def network(findings: dict, store: GraphStore) -> dict:
    genes = [r["gene"] for r in findings["convergence"]["ranking"][:30]]
    nodes, edges, seen = [], [], set()
    for g in genes:
        n = store.resolve(g)
        if not n:
            continue
        nodes.append({"data": {"id": n["id"], "label": g, "kind": "gene"}})
        seen.add(n["id"])
    for nid in list(seen):
        for c in store.claims_for(nid, limit=400):
            other = c["object"] if c["subject"] == nid else c["subject"]
            if c["predicate"] == "interacts_with" and other in seen and c["subject"] < c["object"]:
                edges.append({"data": {"source": c["subject"], "target": c["object"], "kind": "interaction",
                                       "ev": c["evidence_id"]}})
            if c["predicate"] == "targets" and c["object"] == nid:
                d = store.get_node(c["subject"])
                if d and d["id"] not in seen:
                    nodes.append({"data": {"id": d["id"], "label": d["name"], "kind": "drug"}})
                    seen.add(d["id"])
                edges.append({"data": {"source": c["subject"], "target": nid, "kind": "drug", "ev": c["evidence_id"]}})
    return {"nodes": nodes, "edges": edges}


def render(findings_path: Path, out: Path | None = None) -> Path:
    f = json.loads(findings_path.read_text(encoding="utf-8"))
    store = GraphStore()
    ids = collect_ids(f)
    records = {e["evidence_id"]: e for e in store.get_evidence(ids)}
    ev_views = [evidence_view(records[i]) for i in ids if i in records]
    number = {v["id"]: k + 1 for k, v in enumerate(ev_views)}
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    ranking = f["convergence"]["ranking"]
    top = [r["gene"] for r in ranking[:25]]
    targets = [t["gene"] for t in f.get("targets", [])]
    figs = {"volcano": volcano(de, list(dict.fromkeys(targets + top[:10]))),
            "replication": replication_heatmap(de, list(dict.fromkeys(targets + top))[:30]),
            "scores": score_bars(ranking, f["convergence"]["config"]["stream_weights"])}
    from . import report_extra as X
    figs.update({"cells": X.cells_figure(list(dict.fromkeys(targets + top[:20]))[:25]),
                 "pathways": X.pathway_heatmap(), "response": X.response_heatmap()})
    if (PROCESSED / "mr_pqtl.parquet").exists():
        figs["mr"] = X.mr_figure()
    extra = X.tables(list(dict.fromkeys(targets + top)))
    rep = pd.read_parquet(PROCESSED / "replication.parquet").to_dict("records") if (PROCESSED / "replication.parquet").exists() else []
    vpath = findings_path.parent / "verification_final.json"

    env = Environment(loader=FileSystemLoader(ASSETS), autoescape=select_autoescape(["html", "j2"]))

    def cite(ev_ids):
        parts = [f'<a class="cite" href="#{anchor(i)}" title="{html.escape(i)}">{number[i]}</a>' for i in ev_ids if i in number]
        return "<sup>[" + ",".join(parts) + "]</sup>" if parts else ""

    env.filters["cite"] = cite
    env.filters["statement_html"] = lambda st: statement_html(st, records)
    env.filters["evattr"] = lambda ev_ids: " ".join(ev_ids)
    html_out = env.get_template("report.html.j2").render(
        f=f, ranking=ranking, evidence=ev_views, figs=figs, nav=NAV, replication=rep, extra=extra,
        profiles=f["convergence"].get("profiles") or {},
        verification=json.loads(vpath.read_text(encoding="utf-8")) if vpath.exists() else None,
        graph_stats=json.loads((PROCESSED / "graph_stats.json").read_text()),
        css=(ASSETS / "style.css").read_text(encoding="utf-8"), plotly_js=get_plotlyjs(),
        cytoscape_js=(ASSETS / "cytoscape.min.js").read_text(encoding="utf-8"),
        network_json=json.dumps(network(f, store)))
    out = out or REPORTS / "uc_report.html"
    out.write_text(html_out, encoding="utf-8")
    return out
