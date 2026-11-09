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


def tables() -> dict:
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
    sec = pd.read_parquet(PROCESSED / "spatial_sections.parquet")
    out["spatial_sections"] = sec[["gsm", "state", "cells", "assigned_fraction"]].round(2).to_dict("records")
    return out
