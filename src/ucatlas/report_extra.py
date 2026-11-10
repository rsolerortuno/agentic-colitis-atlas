"""M18: figures and tables of the later milestones for the report (cells, pathways, treatment, discovery, validation).
All values come from data/processed and evals; nothing is typed by hand."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .paths import PROCESSED, ROOT


def _fig(fig: go.Figure, h=420) -> str:
    from .render import _fig as base
    return base(fig, h)


def cells_figure(genes: list[str]) -> str:
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    sh = sh[(sh["health"] == "inflamed") & sh["gene"].isin(genes)]
    neu = sh[sh["compartment"] == "Neutrophil"].groupby("gene")["share"].sum().reindex(genes).fillna(0)
    ret = comp["retained"].reindex(genes).fillna(np.nan)
    fig = go.Figure()
    fig.add_trace(go.Bar(y=genes[::-1], x=neu[genes[::-1]], orientation="h", name="share of transcripts in neutrophils",
                         marker_color="#eb6834"))
    fig.add_trace(go.Bar(y=genes[::-1], x=ret[genes[::-1]].clip(0, 1), orientation="h",
                         name="share of the bulk change kept after composition adjustment", marker_color="#2a78d6"))
    fig.update_layout(barmode="group", title="Neutrophil share vs tissue change (top candidates)", xaxis_title="fraction",
                      legend=dict(orientation="h", y=-0.12))
    return _fig(fig, max(420, 22 * len(genes) + 140))


def mr_figure() -> str:
    """Mendelian randomisation: log odds ratio per SD of plasma protein vs -log10 p, every tested protein."""
    r = pd.read_parquet(PROCESSED / "mr_pqtl.parquet")
    y = -np.log10(r["p"].clip(lower=1e-300))
    sig = r["q"] < 0.05
    col = np.where(~sig, "#c9c8c0", np.where(r["beta"] > 0, "#e34948", "#2a78d6"))
    fig = go.Figure(go.Scattergl(x=r["beta"], y=y, mode="markers", marker=dict(size=5, color=col), text=r["gene"],
                                 hovertemplate="%{text}<br>log OR per SD %{x:.2f}<br>-log10 p %{y:.1f}<extra></extra>"))
    lab = r[sig].nsmallest(15, "p")
    fig.add_trace(go.Scatter(x=lab["beta"], y=-np.log10(lab["p"].clip(lower=1e-300)), mode="markers+text", text=lab["gene"],
                             textposition="top center", marker=dict(size=7, color="#1c1c1b"), hoverinfo="skip"))
    fig.update_layout(showlegend=False, title=f"Plasma protein (UKB-PPP cis-pQTL) on UC risk: {len(r):,} proteins",
                      xaxis_title="log odds ratio of UC per SD of protein", yaxis_title="-log10 p")
    return _fig(fig, 460)


def pathway_heatmap() -> str:
    d = pd.read_parquet(PROCESSED / "pathway_activity.parquet")
    m = d.pivot_table(index="source", columns="cohort", values="t")
    fig = go.Figure(go.Heatmap(z=m.values, x=list(m.columns), y=list(m.index), colorscale="RdBu_r", zmid=0,
                               colorbar=dict(title="t"), hovertemplate="%{y} · %{x}<br>t %{z:.1f}<extra></extra>"))
    fig.update_layout(title="Pathway activity, inflamed / UC vs control (PROGENy)")
    return _fig(fig, 460)


def response_heatmap() -> str:
    d = pd.read_parquet(PROCESSED / "response_pathways.parquet")
    d = d[d["kind"] == "pathway"]
    m = d.pivot_table(index="name", columns="arm", values="t_nr_vs_r")
    fig = go.Figure(go.Heatmap(z=m.values, x=list(m.columns), y=list(m.index), colorscale="RdBu_r", zmid=0,
                               colorbar=dict(title="t (NR − R)"),
                               hovertemplate="%{y} · %{x}<br>t %{z:.1f}<extra></extra>"))
    fig.update_layout(title="Pathway activity in non-responders vs responders, by arm")
    return _fig(fig, 460)


def tables(genes: list[str] | None = None) -> dict:
    out = {}
    c = pd.read_parquet(PROCESSED / "communication.parquet")
    c = c[c.get("inflamed", 0) > 0].head(15)
    out["communication"] = [{"ligand": r.ligand, "receptor": r.receptor, "sender": r.sender, "receiver": r.receiver,
                             "score": round(float(r.inflamed), 1), "log2_change": round(float(r.log2_change), 2)}
                            for r in c.itertuples()]
    osm = pd.read_parquet(PROCESSED / "communication.parquet")
    osm = osm[(osm.ligand == "OSM") & osm.receptor.str.contains("OSMR")].head(3)
    out["osm"] = [{"ligand": r.ligand, "receptor": r.receptor, "sender": r.sender, "receiver": r.receiver,
                   "score": round(float(r.inflamed), 1)} for r in osm.itertuples()]
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    out["discovery_counts"] = d["status"].value_counts().to_dict()
    f = d[d.status == "first-in-class opportunity"].sort_values("meta_log2fc", ascending=False).head(20)
    out["first_in_class"] = [{"gene": r.gene, "log2fc": round(r.meta_log2fc, 2), "lineage": r.lineage,
                              "papers": int(r.uc_papers), "tractable": ", ".join(r.tractable),
                              "vital": r.vital_max_ntpm, "retained": round(float(r.retained), 2)} for r in f.itertuples()]
    ev = json.loads((ROOT / "evals" / "results.json").read_text(encoding="utf-8"))
    out["eval"] = [{"score": k, **{x: v[x] for x in ("auroc", "auroc_ci95", "recall_at_100", "median_rank")}}
                   for k, v in ev["clinical_targets"].items()]
    out["n_gold"] = len(ev["gold_sets"]["uc_clinical_targets"])
    qa = json.loads((ROOT / "evals" / "qa_answers.json").read_text(encoding="utf-8"))
    out["qa"] = f"{qa['correct']}/{qa['total']}"
    rg = pd.read_parquet(PROCESSED / "response_genes.parquet")
    out["response_counts"] = rg.groupby("arm").apply(lambda x: int((x.padj < 0.05).sum())).to_dict()
    ab = ROOT / "evals" / "agent_benchmark.json"
    if ab.exists():
        b = json.loads(ab.read_text(encoding="utf-8"))
        out["agent_benchmark"] = [{"arm": arm, **{k: r[k] for k in ("run", "statements", "unsupported", "checklist_score",
                                                                     "checklist_max", "evidence_cited_share", "cost_usd",
                                                                     "minutes", "turns")}}
                                  for arm in ("team", "single") for r in b[arm]]
    sec = pd.read_parquet(PROCESSED / "spatial_sections.parquet")
    cols = [c for c in ("gsm", "study", "state", "cells", "median_counts_per_cell", "assigned_fraction") if c in sec]
    out["spatial_sections"] = sec[cols].round(2).to_dict("records")
    vp = PROCESSED / "spatial_visium.parquet"
    if vp.exists():
        v = pd.read_parquet(vp).set_index("gene")
        keep = [g for g in (genes or []) if g in v.index][:20]
        out["visium"] = [{"gene": g, "compartment": v.at[g, "compartment"], "rho": v.at[g, "compartment_rho"],
                          "second": v.at[g, "second_compartment"], "first_in": int(v.at[g, "sections_compartment_first"]),
                          "sections": int(v.at[g, "uc_sections"]), "log2_ratio": v.at[g, "log2_ratio_uc_vs_healthy"]}
                         for g in keep]
        vs = pd.read_parquet(PROCESSED / "spatial_visium_sections.parquet")
        out["visium_sections"] = {"uc": int((vs.state == "UC_inflamed").sum()), "healthy": int((vs.state == "healthy").sum()),
                                  "spots": int(vs.spots.sum())}
    mp = PROCESSED / "mr_pqtl.parquet"
    if mp.exists():
        m = pd.read_parquet(mp)
        mi = m.set_index("gene")
        from .evaluate import clinical_targets
        clin = set(clinical_targets()["gene"])
        sig = m[m["q"] < 0.05]
        out["mr_counts"] = {"tested": len(m), "significant": len(sig),
                            "colocalising": int(sig["colocalises"].fillna(False).astype(bool).sum()),
                            "clinical_tested": int(m["gene"].isin(clin).sum()),
                            "clinical_significant": sorted(set(sig["gene"]) & clin)}
        out["mr"] = [{"gene": x.gene, "or": round(float(x.or_per_sd), 2),
                      "ci": [round(float(np.exp(x.beta - 1.96 * x.se)), 2), round(float(np.exp(x.beta + 1.96 * x.se)), 2)],
                      "p": f"{x.p:.1e}", "q": f"{x.q:.1e}", "n": int(x.n_instruments),
                      "coloc": "yes" if x.colocalises is True else ("no" if x.colocalises is False else "-"),
                      "clinical": x.gene in clin} for x in sig.sort_values("p").head(25).itertuples()]
        out["mr_candidates"] = [{"gene": g, "tested": g in mi.index,
                                 "or": round(float(mi.at[g, "or_per_sd"]), 2) if g in mi.index else None,
                                 "q": f"{mi.at[g, 'q']:.1e}" if g in mi.index else None}
                                for g in (genes or [])[:15]]
    if "temporal_wide" in ev:
        out["temporal_wide"] = [{"score": k, **{x: v[x] for x in ("n_gold", "auroc", "auroc_ci95", "recall_at_100")}}
                                for k, v in ev["temporal_wide"].items()]
        out["temporal_cutoff"] = ev["gold_sets"].get("cutoff_wide")
    return out
