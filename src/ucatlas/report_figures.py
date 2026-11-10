"""Plotly figures of the report (deterministic; built from data/processed, never from LLM text).

Same colour rules as the kidney atlas: categorical hues in fixed slot order, one-hue sequential ramp for magnitude,
blue <-> red diverging pair with a grey midpoint for polarity, scatter charts with at most three categories. Figures
are transparent; the page re-themes fonts and grids from the CSS tokens.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .paths import LITERATURE, PROCESSED

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQS = [[i / 6, c] for i, c in enumerate(SEQ)]
DIV = [[0, "#184f95"], [0.25, "#6da7ec"], [0.5, "#f0efec"], [0.75, "#ee8f8e"], [1, "#c03130"]]
HEAT = [[0, "#fff5d6"], [0.35, "#f6c177"], [0.7, "#e0603a"], [1, "#9e1b1b"]]
GRAY = "#898781"
TRACK = {"epithelial": CAT[0], "stromal": CAT[1], "immune": CAT[2]}
STREAM = {"expression": CAT[0], "replication": CAT[6], "inflammation": CAT[1], "genetics": CAT[3],
          "literature": CAT[2], "tractability": CAT[4]}
FONT = dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12)
CFG = {"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["lasso2d", "select2d"]}
COHORT = {"GSE193677": "MSCCR", "GSE206285": "UNIFI", "GSE92415": "PURSUIT", "GSE16879": "Leuven IFX",
          "GSE73661": "Leuven VDZ", "GSE109142": "PROTECT", "GSE117993": "RISK", "GSE165512": "GSE165512"}
COMPARTMENT_ORDER = ["Enterocytes", "Goblet", "Stem Cycling", "Tuft", "Enteroendocrine", "Stromal", "Endothelial",
                     "Pericyte", "Glial", "Neutrophil", "Macrophage and DC", "Mast Cell", "Eosinophil", "T Cell and NK",
                     "B Cell", "Plasma", "Erythrocyte"]
SHARE_GROUPS = [("epithelial", "Epithelium", CAT[0], ["Enterocytes", "Goblet", "Stem Cycling", "Tuft", "Enteroendocrine"]),
                ("stromal", "Stroma and vessels", CAT[1], ["Stromal", "Endothelial", "Pericyte", "Glial"]),
                ("neutrophil", "Neutrophil", CAT[7], ["Neutrophil"]),
                ("myeloid", "Macrophage and DC", CAT[2], ["Macrophage and DC", "Mast Cell", "Eosinophil"]),
                ("lymphoid", "T, NK, B and plasma", CAT[3], ["T Cell and NK", "B Cell", "Plasma"])]


def _html(fig: go.Figure, height: int, title: str | None = None, top: int | None = None, **layout) -> str:
    fig.update_layout(height=height, template="plotly_white", font=FONT, paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", margin=dict(l=60, r=20, t=top or (44 if title else 16), b=48),
                      title=dict(text=title, x=0, font=dict(size=14)) if title else None,
                      legend=dict(orientation="h", y=-0.18, font=dict(size=11)), hoverlabel=dict(font=dict(size=12)),
                      **layout)
    fig.update_xaxes(gridcolor="rgba(137,135,129,.22)", zerolinecolor="rgba(137,135,129,.45)", automargin=True)
    fig.update_yaxes(gridcolor="rgba(137,135,129,.22)", zerolinecolor="rgba(137,135,129,.45)", automargin=True)
    return fig.to_html(full_html=False, include_plotlyjs=False, config=CFG)


def _p(name):
    return pd.read_parquet(PROCESSED / name)


# ------------------------------------------------------------------ convergence
def stream_stack(ranking: list[dict], weights: dict) -> str:
    rows = ranking[:25][::-1]
    fig = go.Figure()
    for k, w in weights.items():
        fig.add_trace(go.Bar(y=[r["gene"] for r in rows], x=[w * r["components"][k] for r in rows], name=k, orientation="h",
                             marker_color=STREAM.get(k, GRAY), marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.2f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted contribution to the convergence score")
    return _html(fig, 640, "What drives each candidate's score? (top 25)")


def component_heatmap(ranking: list[dict], keys: list[str]) -> str:
    rows = ranking[:25]
    fig = go.Figure(go.Heatmap(z=[[r["components"][k] for k in keys] for r in rows], x=keys, y=[r["gene"] for r in rows],
                               colorscale=SEQS, zmin=0, zmax=1, colorbar=dict(title="component", thickness=12),
                               hovertemplate="%{y} · %{x}: %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed")
    return _html(fig, 640, "Evidence streams per candidate (0 – 1)")


# ------------------------------------------------------------------ expression
def volcano(contrast: str, highlight: set[str], title: str, xlab: str) -> str:
    d = _p("de_results.parquet")
    d = d[(d.cohort == "GSE193677") & (d.contrast == contrast)].dropna(subset=["padj"]).copy()
    d["mlog"] = -np.log10(d["padj"].clip(lower=1e-300))
    bg, hi = d[~d.gene.isin(highlight)], d[d.gene.isin(highlight)]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=bg.log2fc, y=bg.mlog, mode="markers", name="all genes",
                               marker=dict(size=3, color="rgba(137,135,129,.35)"), text=bg.gene,
                               hovertemplate="%{text}<br>log2FC %{x:.2f} · -log10 padj %{y:.1f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=hi.log2fc, y=hi.mlog, mode="markers+text", name="report candidates and top 10",
                             text=hi.gene, textposition="top center", textfont=dict(size=10),
                             marker=dict(size=8, color=CAT[0], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{text}<br>log2FC %{x:.2f} · -log10 padj %{y:.1f}<extra></extra>"))
    fig.update_xaxes(title=xlab)
    fig.update_yaxes(title="−log10 adjusted p")
    return _html(fig, 430, title)


def cohort_scatter(a: str, b: str, highlight: set[str], title: str | None = None) -> tuple[str, float, int]:
    de = _p("de_results.parquet")
    pick = lambda c: de[(de.cohort == c) & de.contrast.isin(["inflamed_vs_control", "uc_vs_control"])].set_index("gene")
    j = pick(a)[["log2fc", "padj"]].join(pick(b)[["log2fc", "padj"]], rsuffix="_b", how="inner").dropna()
    rho = float(j.log2fc.corr(j.log2fc_b, method="spearman"))
    s = j.sample(min(len(j), 6000), random_state=1)
    fig = go.Figure(go.Scattergl(x=s.log2fc, y=s.log2fc_b, mode="markers", name="genes",
                                 marker=dict(color=GRAY, size=3, opacity=.35), text=s.index,
                                 hovertemplate="%{text}: " + COHORT[a] + " %{x:.2f}, " + COHORT[b] + " %{y:.2f}<extra></extra>"))
    h = j[j.index.isin(highlight)]
    fig.add_trace(go.Scatter(x=h.log2fc, y=h.log2fc_b, mode="markers+text", name="report candidates", text=h.index,
                             textposition="top center", textfont=dict(size=9), marker=dict(color=CAT[0], size=8),
                             hovertemplate="%{text}<extra></extra>"))
    lim = [float(min(j.log2fc.min(), j.log2fc_b.min())), float(max(j.log2fc.max(), j.log2fc_b.max()))]
    fig.add_trace(go.Scatter(x=lim, y=lim, mode="lines", line=dict(color=GRAY, dash="dot", width=1), name="same change",
                             hoverinfo="skip"))
    fig.update_xaxes(title=f"{COHORT[a]} log2FC (UC vs control)")
    fig.update_yaxes(title=f"{COHORT[b]} log2FC")
    return _html(fig, 430, title or f"{COHORT[a]} vs {COHORT[b]} (Spearman ρ = {rho:.2f} over {len(j):,} genes)"), rho, len(j)


def replication_heatmap(genes: list[str]) -> str:
    de = _p("de_results.parquet")
    d = de[de.contrast.isin(["inflamed_vs_control", "uc_vs_control"]) & de.gene.isin(genes)]
    m = d.pivot_table(index="gene", columns="cohort", values="log2fc").reindex(index=genes, columns=list(COHORT))
    lim = float(np.nanpercentile(np.abs(m.values), 95)) if m.size else 1
    fig = go.Figure(go.Heatmap(z=m.values, x=[COHORT[c] for c in m.columns], y=m.index, colorscale=DIV, zmid=0,
                               zmin=-lim, zmax=lim, colorbar=dict(title="log2FC", thickness=12),
                               hovertemplate="%{y} · %{x}<br>log2FC %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(380, 18 * len(genes) + 120), "Top candidates across the eight cohorts (UC vs control)")


def activity_bars(genes: list[str]) -> str:
    a = _p("activity.parquet")
    a = a[a.gene.isin(genes)]
    w = a.pivot_table(index="gene", columns="cohort", values="rho").reindex(genes).dropna(how="all")
    fig = go.Figure()
    for c, color in (("GSE193677", CAT[0]), ("GSE73661", CAT[1])):
        if c in w:
            fig.add_trace(go.Bar(y=w.index, x=w[c], name=COHORT[c], orientation="h", marker_color=color,
                                 hovertemplate="%{y} · " + COHORT[c] + " ρ %{x:.2f}<extra></extra>"))
    fig.update_layout(barmode="group")
    fig.add_vline(x=0, line_color=GRAY, line_width=1)
    fig.add_annotation(x=0.69, y=1.0, xref="x", yref="paper", text="higher in more active disease →", showarrow=False,
                       xanchor="right", yanchor="bottom", font=dict(size=11, color="#4d4c49"))
    fig.add_annotation(x=-0.69, y=1.0, xref="x", yref="paper", text="← lower in more active disease", showarrow=False,
                       xanchor="left", yanchor="bottom", font=dict(size=11, color="#4d4c49"))
    fig.update_xaxes(title="Spearman ρ between the gene and the endoscopic Mayo subscore", range=[-0.7, 0.7])
    fig.update_yaxes(autorange="reversed", title="gene")
    return _html(fig, max(380, 26 * len(w) + 130), "Does the gene follow disease activity? (MSCCR and Leuven)", top=70)


def activity_summary(genes: list[str], n: int = 3) -> dict:
    a = _p("activity.parquet")
    a = a[a.gene.isin(genes) & (a.cohort == "GSE193677")].dropna(subset=["rho"]).sort_values("rho")
    return {"up": [(r.gene, round(float(r.rho), 2)) for r in a.tail(n).iloc[::-1].itertuples() if r.rho > 0],
            "down": [(r.gene, round(float(r.rho), 2)) for r in a.head(n).itertuples() if r.rho < 0]}


# ------------------------------------------------------------------ cells
def dotplot(genes: list[str]) -> str:
    sh = _p("sc3_celltype_share.parquet")
    d = sh[(sh.health == "inflamed") & sh.gene.isin(genes)]
    d = d.groupby(["gene", "compartment"]).apply(lambda x: pd.Series({
        "pct": float((x.pct_expressing * x.n_cells).sum() / x.n_cells.sum()),
        "mean": float(np.log1p((x.mean_cp10k * x.n_cells).sum() / x.n_cells.sum())),
        "share": float(x.share.sum()), "cells": int(x.n_cells.sum())}), include_groups=False).reset_index()
    fig = go.Figure(go.Scatter(
        x=d.compartment, y=d.gene, mode="markers",
        marker=dict(size=np.sqrt(d.pct.clip(lower=0.5)) * 3.2, color=d["mean"], cmin=0, cmax=float(d["mean"].max()),
                    colorscale=SEQS, colorbar=dict(title="log(1 + CP10K)", thickness=12),
                    line=dict(width=0.6, color="rgba(137,135,129,.6)")),
        customdata=np.c_[d.pct, d.share, d.cells],
        hovertemplate="%{y} in %{x}<br>%{customdata[0]:.1f} % of %{customdata[2]:,.0f} cells · "
                      "%{customdata[1]:.1%} of the gene's transcripts<extra></extra>"))
    fig.update_xaxes(categoryorder="array", categoryarray=COMPARTMENT_ORDER, tickangle=-40)
    fig.update_yaxes(categoryorder="array", categoryarray=genes[::-1])
    return _html(fig, max(420, 20 * len(genes) + 160),
                 "Who expresses each candidate in inflamed mucosa? Dot size = % cells, colour = mean expression",
                 showlegend=False)


def share_stack(genes: list[str]) -> str:
    sh = _p("sc3_celltype_share.parquet")
    d = sh[(sh.health == "inflamed") & sh.gene.isin(genes)]
    piv = d.pivot_table(index="gene", columns="compartment", values="share", aggfunc="sum").fillna(0).reindex(genes).dropna(how="all")
    series = {k: piv[[c for c in comps if c in piv]].sum(axis=1) for k, _, _, comps in SHARE_GROUPS}
    order = series["neutrophil"].sort_values(ascending=False).index.tolist()
    fig = go.Figure()
    for k, label, color, _ in SHARE_GROUPS:
        s = series[k].reindex(order)
        fig.add_trace(go.Bar(y=order, x=s.values, name=label, orientation="h", marker_color=color,
                             marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y}: %{x:.1%} in " + label + "<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(tickformat=".0%", range=[0, 1])
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(420, 19 * len(order) + 120), "Share of each gene's transcripts by cell type (inflamed, SCP3755)")


def neutrophil_vs_retained(genes: list[str]) -> str:
    comp = _p("composition_effect.parquet").set_index("gene")
    sh = _p("sc3_celltype_share.parquet")
    sh = sh[(sh.health == "inflamed") & sh.gene.isin(genes)]
    neu = sh[sh.compartment == "Neutrophil"].groupby("gene")["share"].sum().reindex(genes).fillna(0)
    ret = comp["retained"].reindex(genes)
    fig = go.Figure()
    fig.add_trace(go.Bar(y=genes[::-1], x=neu[genes[::-1]], orientation="h", name="share of transcripts in neutrophils",
                         marker_color=CAT[1]))
    fig.add_trace(go.Bar(y=genes[::-1], x=ret[genes[::-1]].clip(0, 1), orientation="h",
                         name="share of the bulk change kept after composition adjustment", marker_color=CAT[0]))
    fig.update_layout(barmode="group")
    fig.update_xaxes(title="fraction")
    return _html(fig, max(420, 22 * len(genes) + 140), "Neutrophil share vs tissue change (top candidates)")


def sc_replication_figure() -> tuple[str, dict]:
    d = _p("sc_replication.parquet").dropna(subset=["share_2019", "share_2023"])
    lins = d.lineage_2019.value_counts().index[:3].tolist()
    fig = go.Figure()
    for lin, color in zip(lins, CAT):
        s = d[d.lineage_2019 == lin]
        fig.add_trace(go.Scatter(x=s.share_2019, y=s.share_2023, mode="markers", name=f"dominant in Smillie: {lin}",
                                 marker=dict(color=color, size=6, opacity=.7), text=s.gene,
                                 hovertemplate="%{text}: %{x:.2f} vs %{y:.2f}<extra></extra>"))
    o = d[~d.lineage_2019.isin(lins)]
    fig.add_trace(go.Scatter(x=o.share_2019, y=o.share_2023, mode="markers", name="other lineages",
                             marker=dict(color=GRAY, size=5, opacity=.5), text=o.gene + " (" + o.lineage_2019 + ")",
                             hovertemplate="%{text}: %{x:.2f} vs %{y:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=GRAY, dash="dot", width=1), showlegend=False,
                             hoverinfo="skip"))
    fig.update_xaxes(title="share in the dominant lineage, Smillie et al. 2019", range=[0, 1.02])
    fig.update_yaxes(title="share in the same lineage, second atlas", range=[0, 1.02])
    rho = float(d.share_2019.corr(d.share_2023, method="spearman"))
    stats_ = {"genes": int(len(d)), "same_lineage": int(d.same_lineage.sum()), "rho": round(rho, 2)}
    return _html(fig, 430, f"Does a second single-cell atlas agree? (ρ = {rho:.2f}, {len(d)} genes)"), stats_


def visium_map(section: str, title: str) -> str | None:
    from .spatial_deconv import CELL_GROUPS, VIEWER_GENES
    dec = _p("spatial_visium_deconvolution.parquet")
    d = dec[dec.section == section]
    if d.empty:
        return None
    from .spatial_visium import load
    m, genes = load(section)
    gi = {g: i for i, g in enumerate(genes)}
    lib = np.asarray(m.sum(axis=1)).ravel()
    variables = [(g, f"{g.replace('_', ' ')} share", False) for g in ("Epithelial", "Neutrophil", "Inflammatory_fibroblast",
                                                                     "Smooth_muscle", "Myeloid", "Plasma") if g in CELL_GROUPS]
    variables += [(g, f"{g} expression (CPM)", True) for g in VIEWER_GENES[:8] if g in gi]
    fig = go.Figure()
    buttons = []
    for i, (col, label, is_gene) in enumerate(variables):
        v = (m[:, gi[col]].toarray().ravel() / np.maximum(lib, 1) * 1e6) if is_gene else d[col].to_numpy(float)
        hi = np.quantile(v, 0.99) if is_gene else max(float(np.quantile(v, 0.99)), 0.05)
        fig.add_trace(go.Scattergl(x=d.x, y=d.y, mode="markers", visible=i == 0, name=label,
                                   marker=dict(size=5.5, color=v, cmin=0, cmax=max(hi, 1e-6), colorscale=SEQS,
                                               colorbar=dict(title=label.split(" (")[0], thickness=12)),
                                   hovertemplate=label + ": %{marker.color:.2f}<extra></extra>"))
        buttons.append(dict(label=label.replace(" share", "").replace(" expression (CPM)", " (gene)"), method="update",
                            args=[{"visible": [j == i for j in range(len(variables))]}]))
    fig.update_layout(updatemenus=[dict(buttons=buttons, direction="down", x=0, xanchor="left", y=1.0, yanchor="bottom",
                                        bgcolor="rgba(255,255,255,.9)", bordercolor="#c9c8c0", pad=dict(t=0, b=6))])
    fig.update_yaxes(autorange="reversed", scaleanchor="x", showgrid=False, showticklabels=False, zeroline=False)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    return _html(fig, 520, title, top=104, showlegend=False)


def neighbor_heatmap(genes: list[str]) -> str:
    nb = _p("spatial_visium_neighbors.parquet")
    nb = nb[nb.gene.isin(genes)]
    z = nb.pivot(index="gene", columns="cell_group", values="log2_enrichment").reindex(genes).dropna(how="all")
    p = nb.pivot(index="gene", columns="cell_group", values="padj").reindex(index=z.index, columns=z.columns)
    text = z.round(2).astype(str).where(p < 0.05, "")
    fig = go.Figure(go.Heatmap(z=z.values, x=[c.replace("_", " ") for c in z.columns], y=list(z.index), zmid=0,
                               colorscale=DIV, zmin=-1.5, zmax=1.5, text=text.values, texttemplate="%{text}",
                               textfont=dict(size=9), colorbar=dict(title="log2 obs/exp", thickness=12),
                               hovertemplate="%{y}-high spots · neighbours %{x}<br>log2 obs/exp %{z:.2f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(420, 20 * len(z) + 140), "Cell groups around gene-high spots (Visium; labelled: padj < 0.05)")


def expresser_stack(genes: list[str]) -> str:
    ex = _p("spatial_visium_expressers.parquet")
    ex = ex[ex.gene.isin(genes)]
    piv = ex.pivot(index="gene", columns="cell_group", values="in_situ_contribution").reindex(genes).dropna(how="all")
    order = piv.idxmax(axis=1).sort_values().index.tolist()
    fig = go.Figure()
    from .spatial_deconv import CELL_GROUPS
    colors = dict(zip(CELL_GROUPS, CAT + ["#5f5d58", "#b7b5ae", "#3b8ea5", "#a35d9b"]))
    for g in CELL_GROUPS:
        if g in piv and piv[g].max() > 0.02:
            fig.add_trace(go.Bar(y=order, x=piv.loc[order, g], name=g.replace("_", " "), orientation="h", marker_color=colors[g],
                                 marker_line=dict(width=1, color="rgba(255,255,255,.9)"),
                                 hovertemplate="%{y}: %{x:.0%} from " + g.replace("_", " ") + "<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(tickformat=".0%", range=[0, 1], title="share of the gene's in-situ signal")
    fig.update_yaxes(autorange="reversed")
    return _html(fig, max(420, 20 * len(order) + 130), "Which cell explains each gene in the tissue? (deconvolved Visium spots)")


# ------------------------------------------------------------------ structure
def _pockets() -> pd.DataFrame:
    d = _p("structure_pockets.parquet")
    return d[d["rank"] == 1].copy()


def pocket_bubble(report_genes: list[str], discovery_genes: list[str]) -> str:
    d = _pockets()

    def kind(r):
        if r.control:
            return "UC small-molecule target (control)"
        return "discovery shortlist" if r.gene in discovery_genes else "report candidate"
    d["kind"] = [kind(r) for r in d.itertuples()]
    fig = go.Figure()
    for k, c in {"report candidate": CAT[0], "discovery shortlist": CAT[1], "UC small-molecule target (control)": CAT[2]}.items():
        s = d[d.kind == k]
        fig.add_trace(go.Scatter(x=s.volume_a3, y=s.pocket_score, mode="markers+text", text=s.gene, name=k,
                                 textposition="top center", textfont=dict(size=9),
                                 marker=dict(size=6 + 1.2 * s.max_depth_a, color=c, opacity=0.85, line=dict(width=1.5, color="#fff")),
                                 customdata=np.c_[s.max_depth_a, s.mean_plddt, s.avg_hydropathy],
                                 hovertemplate="%{text}<br>%{x:.0f} Å³ · depth %{customdata[0]:.1f} Å · pLDDT %{customdata[1]:.0f}"
                                               "<br>hydropathy %{customdata[2]:.2f} · pocket score %{y:.2f}<extra></extra>"))
    fig.update_xaxes(title="pocket volume (Å³, log)", type="log")
    fig.update_yaxes(title="pocket score", range=[0, 1.05])
    return _html(fig, 470, "Best pocket per protein: size and score (bubble = depth)")


def pocket_components(weights: dict, n: int = 24) -> str:
    d = _pockets()
    d = d.sort_values("pocket_score", ascending=False).head(n)[::-1]
    fig = go.Figure()
    for i, k in enumerate(weights):
        fig.add_trace(go.Bar(y=d.gene, x=weights[k] * d[f"comp_{k}"], name=k, orientation="h", marker_color=CAT[i],
                             marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.3f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted component (sums to the pocket score)")
    return _html(fig, max(420, 22 * len(d) + 110), "Why a pocket scores high: model confidence, geometry, hydrophobicity")


def ligand_validation(st: dict) -> str | None:
    rows = [(r["gene"], r["ligand_site_validation"]) for r in st["results"]
            if r["modeled"] and r["ligand_site_validation"].get("checked")]
    if not rows:
        return None
    rows.sort(key=lambda x: x[1]["recall"])
    genes = [g for g, _ in rows]
    det, top = [v["recall"] for _, v in rows], [v["recall_in_top_pockets"] for _, v in rows]
    fig = go.Figure()
    for g, a, b in zip(genes, det, top):
        fig.add_trace(go.Scatter(x=[b, a], y=[g, g], mode="lines", line=dict(color="rgba(137,135,129,.55)", width=2),
                                 showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=det, y=genes, mode="markers", name="all detected cavities",
                             marker=dict(size=9, color=CAT[0], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{y}: cavities cover %{x:.0%} of ligand residues<extra></extra>"))
    fig.add_trace(go.Scatter(x=top, y=genes, mode="markers", name="top 3 pockets by score",
                             marker=dict(size=9, color=CAT[1], line=dict(width=1.5, color="#fff")),
                             hovertemplate="%{y}: top 3 cover %{x:.0%}<extra></extra>"))
    fig.update_xaxes(title="fraction of experimentally ligand-bound residues covered", tickformat=".0%", range=[-0.02, 1.02])
    return _html(fig, max(420, 20 * len(genes) + 120), "Check against PDBe-KB ligand sites")


def control_ranking() -> str | None:
    d = _pockets().sort_values("pocket_score", ascending=False).reset_index(drop=True)
    if not d.control.any():
        return None
    fig = go.Figure(go.Bar(x=d.gene, y=d.pocket_score, marker_color=[CAT[2] if c else "rgba(137,135,129,.45)" for c in d.control],
                           marker_cornerradius=3, hovertemplate="%{x}: %{y:.2f}<extra></extra>"))
    fig.update_xaxes(tickangle=-60, tickfont=dict(size=9))
    fig.update_yaxes(title="pocket score", range=[0, 1])
    return _html(fig, 360, "Best pocket score of every protein (green = UC small-molecule target)", showlegend=False)


# ------------------------------------------------------------------ discovery
def discovery_scatter() -> str:
    d = _p("discovery_scored.parquet")
    fig = go.Figure()
    for track, color in TRACK.items():
        s = d[d.track == track]
        fic = s.status == "first-in-class opportunity"
        for flag, sym, name in ((fic, "star", "first-in-class"), (~fic, "circle", "other")):
            t = s[flag]
            if t.empty:
                continue
            fig.add_trace(go.Scatter(
                x=t.uc_papers + 1, y=t.meta_log2fc, mode="markers", name=f"{track} · {name}",
                marker=dict(symbol=sym, size=5 + 14 * t.retained.fillna(0).clip(0, 1), color=color,
                            opacity=0.9 if sym == "star" else 0.4, line=dict(width=0.8, color="#fff")),
                text=t.gene, customdata=np.c_[t.novelty_score, t.lineage, t.status],
                hovertemplate="<b>%{text}</b> (%{customdata[1]})<br>%{customdata[2]}<br>novelty %{customdata[0]:.2f} · "
                              "UC papers %{x:.0f}<br>meta log2FC %{y:.2f}<extra></extra>"))
    top = d[d.status == "first-in-class opportunity"].sort_values("novelty_score", ascending=False).groupby("track").head(4)
    for r in top.itertuples():
        fig.add_annotation(x=np.log10(r.uc_papers + 1), y=r.meta_log2fc, text=r.gene, showarrow=False, yshift=11,
                           font=dict(size=9), xref="x", yref="y")
    fig.update_xaxes(title="UC papers + 1 (Europe PMC, log)", type="log")
    fig.update_yaxes(title="meta-analysis log2FC (UC vs control)")
    return _html(fig, 520, "Discovery map: genes up in UC with little UC literature (marker size = tissue change kept)")


def discovery_bars(weights: dict, n: int = 18) -> str:
    d = _p("discovery_scored.parquet")
    d = d[d.status == "first-in-class opportunity"].sort_values("novelty_score", ascending=False).head(n)[::-1]
    tot = sum(weights.values())
    fig = go.Figure()
    for i, (k, w) in enumerate(weights.items()):
        fig.add_trace(go.Bar(y=d.gene, x=w * d[f"c_{k}"] / tot, name=k.replace("_", " "), orientation="h", marker_color=CAT[i],
                             marker_line=dict(width=1.5, color="rgba(255,255,255,.9)"),
                             hovertemplate="%{y} · " + k + ": %{x:.3f}<extra></extra>"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="weighted component (weight × component ÷ total weight)")
    return _html(fig, max(420, 24 * len(d) + 120), "What makes a first-in-class candidate score high")


def safety_scatter() -> str:
    d = _p("discovery_scored.parquet")
    d = d[d.status.str.startswith("novel") | (d.status == "first-in-class opportunity")].dropna(subset=["vital_max_ntpm"])
    fig = go.Figure()
    for track, color in TRACK.items():
        s = d[d.track == track]
        fig.add_trace(go.Scatter(x=s.vital_max_ntpm + 0.1, y=s.novelty_score, mode="markers", name=track,
                                 marker=dict(size=7, color=color, opacity=0.75, line=dict(width=0.8, color="#fff")), text=s.gene,
                                 hovertemplate="<b>%{text}</b><br>%{x:.0f} nTPM in a vital organ · novelty %{y:.2f}<extra></extra>"))
    fig.add_vline(x=50, line=dict(color=GRAY, dash="dot", width=1))
    fig.update_xaxes(title="highest expression in a vital organ (HPA nTPM, log)", type="log")
    fig.update_yaxes(title="novelty score")
    return _html(fig, 400, "Off-target risk of novel candidates (right of the line = expressed in vital organs)")


# ------------------------------------------------------------------ treatment and subgroups
CHECK_GENES = ["IL11", "IL13RA2", "PTGS2", "STC1", "TNFRSF11B", "OSM", "TREM1"]


def treatment_controls() -> str:
    r = _p("response_genes.parquet")
    r = r[r.gene.isin(CHECK_GENES)]
    m = r.pivot_table(index="gene", columns="arm", values="t").reindex(CHECK_GENES)
    arms = [a for a in ["infliximab", "golimumab", "placebo (PURSUIT)", "ustekinumab", "placebo (UNIFI)", "vedolizumab"] if a in m]
    fig = go.Figure()
    for i, a in enumerate(arms):
        fig.add_trace(go.Bar(x=m.index, y=m[a], name=a, marker_color=CAT[i] if "placebo" not in a else ["#b7b5ae", "#898781"][i % 2],
                             hovertemplate="%{x} · " + a + ": t %{y:.1f}<extra></extra>"))
    fig.add_hline(y=0, line=dict(color=GRAY, width=1))
    fig.update_layout(barmode="group")
    fig.update_yaxes(title="t, non-responders minus responders (before treatment)")
    return _html(fig, 430, "Check: published infliximab non-response genes, every arm")


def treatment_arms(highlight: set[str]) -> str:
    r = _p("response_genes.parquet")
    w = r[r.arm.isin(["ustekinumab", "placebo (UNIFI)"])].pivot_table(index="gene", columns="arm", values="t").dropna()
    s = w.sample(min(len(w), 6000), random_state=1)
    fig = go.Figure(go.Scattergl(x=s["placebo (UNIFI)"], y=s["ustekinumab"], mode="markers", name="genes",
                                 marker=dict(color=GRAY, size=3, opacity=.35), text=s.index,
                                 hovertemplate="%{text}: placebo t %{x:.1f}, ustekinumab t %{y:.1f}<extra></extra>"))
    h = w[w.index.isin(highlight)]
    fig.add_trace(go.Scatter(x=h["placebo (UNIFI)"], y=h["ustekinumab"], mode="markers+text", name="report candidates",
                             text=h.index, textposition="top center", textfont=dict(size=9), marker=dict(color=CAT[0], size=8),
                             hovertemplate="%{text}<extra></extra>"))
    lim = [float(w.min().min()), float(w.max().max())]
    fig.add_trace(go.Scatter(x=lim, y=lim, mode="lines", line=dict(color=GRAY, dash="dot", width=1), name="same in both arms",
                             hoverinfo="skip"))
    rho = float(w["placebo (UNIFI)"].corr(w["ustekinumab"], method="spearman"))
    fig.update_xaxes(title="placebo arm: t, non-responders minus responders")
    fig.update_yaxes(title="ustekinumab arm: t")
    return _html(fig, 430, f"UNIFI: drug or severity? Every gene, ustekinumab vs placebo (ρ = {rho:.2f})")


def response_heatmap() -> str:
    d = _p("response_pathways.parquet")
    d = d[d.kind == "pathway"]
    m = d.pivot_table(index="name", columns="arm", values="t_nr_vs_r")
    fig = go.Figure(go.Heatmap(z=m.values, x=list(m.columns), y=list(m.index), colorscale=DIV, zmid=0,
                               colorbar=dict(title="t (NR − R)", thickness=12),
                               hovertemplate="%{y} · %{x}<br>t %{z:.1f}<extra></extra>"))
    return _html(fig, 460, "Pathway activity in non-responders vs responders, by arm")


def subgroup_counts(sub: dict) -> str:
    c = sub["counts"]
    x = [COHORT[k] + " (paediatric)" for k in c["paediatric"]] + [COHORT[k] for k in c["adult"]]
    y = list(c["paediatric"].values()) + list(c["adult"].values())
    col = [CAT[1]] * len(c["paediatric"]) + [CAT[0]] * len(c["adult"])
    s = sub["site_counts"]
    x += ["MSCCR rectum", "MSCCR proximal colon"]
    y += [s["rectum"], s["proximal"]]
    col += [CAT[2], CAT[2]]
    fig = go.Figure(go.Bar(x=x, y=y, marker_color=col, marker_cornerradius=4,
                           hovertemplate="%{x}: %{y} genes up (padj < 0.05, log2FC ≥ 1)<extra></extra>"))
    fig.update_yaxes(title="genes up in UC (padj < 0.05, log2FC ≥ 1)")
    fig.update_xaxes(tickangle=-30)
    return _html(fig, 400, "How many genes rise, by cohort and biopsy site", showlegend=False)


def activity_axis() -> str:
    d = _p("subgroup_activity_axis.parquet")
    order = ["control", "Mayo 0", "Mayo 1", "Mayo 2", "Mayo 3"]
    fig = go.Figure()
    for i, g in enumerate([o for o in order if o in set(d.mayo_group)]):
        s = d[d.mayo_group == g]
        fig.add_trace(go.Violin(x=[g] * len(s), y=s.neutrophil_score, name=g, box_visible=True, meanline_visible=False,
                                line_color=GRAY if g == "control" else CAT[0], fillcolor="rgba(42,120,214,.25)" if g != "control"
                                else "rgba(137,135,129,.25)", points=False,
                                hovertemplate=g + ": %{y:.2f}<extra></extra>"))
    fig.update_yaxes(title="neutrophil lineage score (z, per biopsy)")
    return _html(fig, 400, "Histological activity axis: neutrophil score by endoscopic Mayo (MSCCR biopsies)", showlegend=False)


# ------------------------------------------------------------------ biologics and immune
def biologics_figure() -> str:
    d = _p("biologics.parquet")
    d = d[d.in_universe | d.control]
    fig = go.Figure()
    s = d[~d.control]
    fig.add_trace(go.Scatter(x=s.meta_log2fc, y=s.patients_carrying_pct, mode="markers", name="genes up in UC",
                             marker=dict(size=4 + 10 * s.c_surface.fillna(0), color=s.vital_max_ntpm.fillna(0).clip(0, 100), colorscale=HEAT,
                                         cmin=0, cmax=100, colorbar=dict(title="vital-organ nTPM", thickness=12),
                                         line=dict(width=0.5, color="#fff")),
                             text=s.gene, customdata=np.c_[s.suitability, s.vital_max_ntpm.fillna(-1), s.c_surface],
                             hovertemplate="<b>%{text}</b><br>log2FC %{x:.2f} · %{y:.0f} % of patients<br>vital organs %{customdata[1]:.0f} nTPM"
                                           " · surface %{customdata[2]:.1f} · score %{customdata[0]:.2f}<extra></extra>"))
    c = d[d.control]
    fig.add_trace(go.Scatter(x=c.meta_log2fc, y=c.patients_carrying_pct, mode="markers+text", name="UC antibody targets (controls)",
                             text=c.gene, textposition="top center", marker=dict(symbol="diamond", size=10, color=CAT[0]),
                             hovertemplate="%{text}<extra></extra>"))
    top = s.head(10)
    fig.add_trace(go.Scatter(x=top.meta_log2fc, y=top.patients_carrying_pct, mode="text", text=top.gene, textposition="bottom center",
                             textfont=dict(size=9), showlegend=False, hoverinfo="skip"))
    fig.update_xaxes(title="meta-analysis log2FC (UC vs control)")
    fig.update_yaxes(title="inflamed patients above the control 95th percentile (%)")
    return _html(fig, 470, "Antibody and biologic suitability (marker size = surface or secreted)")


def immune_classes(imm: dict) -> str:
    m = imm["lineage_means_by_class"]
    order = [c for c in ("low", "intermediate", "high") if c in m]
    lins = list(next(iter(m.values())))
    fig = go.Figure()
    for i, lin in enumerate(lins):
        fig.add_trace(go.Bar(x=order, y=[m[c][lin] for c in order], name=lin, marker_color=CAT[i],
                             hovertemplate=lin + " · %{x}: %{y:.2f}<extra></extra>"))
    fig.update_layout(barmode="group")
    fig.update_yaxes(title="mean lineage score (z)")
    return _html(fig, 400, "Immune lineages in thirds of patients by immune score (MSCCR)")


# ------------------------------------------------------------------ validation
def composition_shift(genes: list[str]) -> str | None:
    c = _p("composition_effect.parquet").set_index("gene")
    d = c.reindex([g for g in genes if g in c.index]).dropna(subset=["adjusted_log2fc"])
    if d.empty:
        return None
    d = d.sort_values("adjusted_log2fc")
    fig = go.Figure()
    for g, r in d.iterrows():
        fig.add_trace(go.Scatter(x=[r.unadjusted_log2fc, r.adjusted_log2fc], y=[g, g], mode="lines",
                                 line=dict(color=GRAY, width=1.5), showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=d.unadjusted_log2fc, y=d.index, mode="markers", name="bulk biopsy, inflamed vs control",
                             marker=dict(color=GRAY, size=8), hovertemplate="%{y}: bulk %{x:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=d.adjusted_log2fc, y=d.index, mode="markers", name="after adjusting for the cell mix",
                             marker=dict(color=CAT[0], size=9), error_x=dict(type="data", array=1.96 * d.adjusted_se, color=CAT[0]),
                             hovertemplate="%{y}: adjusted %{x:.2f}<extra></extra>"))
    fig.update_xaxes(title="log2FC (MSCCR, patient means)")
    return _html(fig, max(380, 22 * len(d) + 120), "Cell influx or tissue change? Bulk vs composition-adjusted change")


def genetics_agreement() -> str | None:
    p = _p("mr_pqtl.parquet").set_index("gene")
    c = _p("mr_colon.parquet")
    c = c[c.layer == "colon_tissue"].set_index("gene")
    j = p[["beta", "q"]].join(c[["beta", "q"]], rsuffix="_colon", how="inner")
    if j.empty:
        return None
    sig = (j.q < 0.05) | (j.q_colon < 0.05)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=j.loc[~sig, "beta"], y=j.loc[~sig, "beta_colon"], mode="markers", name="neither significant",
                             marker=dict(color=GRAY, size=4, opacity=.4), text=j.index[~sig],
                             hovertemplate="%{text}: plasma %{x:.2f}, colon %{y:.2f}<extra></extra>"))
    s = j[sig]
    fig.add_trace(go.Scatter(x=s.beta, y=s.beta_colon, mode="markers+text", name="significant in at least one",
                             marker=dict(color=CAT[1], size=8), text=s.index, textposition="top center", textfont=dict(size=9),
                             hovertemplate="%{text}: plasma %{x:.2f}, colon %{y:.2f}<extra></extra>"))
    fig.add_hline(y=0, line=dict(color=GRAY, width=1))
    fig.add_vline(x=0, line=dict(color=GRAY, width=1))
    same = float((np.sign(j.beta) == np.sign(j.beta_colon)).mean())
    fig.update_xaxes(title="log OR per SD of plasma protein (UKB-PPP)")
    fig.update_yaxes(title="log OR per SD of colon expression (GTEx, CEDAR)")
    return _html(fig, 450, f"Blood and gut genetics: {len(j)} genes tested in both ({same:.0%} same sign)")


def mr_volcano(path: str, title: str, layer: str | None = None, xlab: str = "log odds ratio of UC per SD") -> str:
    r = _p(path)
    if layer:
        r = r[r.layer == layer]
    y = -np.log10(r.p.clip(lower=1e-300))
    sig = r.q < 0.05
    col = np.where(~sig, "rgba(137,135,129,.45)", np.where(r.beta > 0, "#c03130", "#184f95"))
    fig = go.Figure(go.Scattergl(x=r.beta, y=y, mode="markers", marker=dict(size=5, color=col), text=r.gene,
                                 hovertemplate="%{text}<br>log OR %{x:.2f}<br>-log10 p %{y:.1f}<extra></extra>"))
    lab = r[sig].nsmallest(15, "p")
    fig.add_trace(go.Scatter(x=lab.beta, y=-np.log10(lab.p.clip(lower=1e-300)), mode="markers+text", text=lab.gene,
                             textposition="top center", textfont=dict(size=9), marker=dict(size=7, color="#1c1c1b"),
                             hoverinfo="skip"))
    fig.update_xaxes(title=xlab)
    fig.update_yaxes(title="−log10 p")
    return _html(fig, 450, f"{title}: {len(r):,} genes", showlegend=False)


def mr_celltypes(genes: list[str]) -> str | None:
    c = _p("mr_colon_celltypes.parquet")
    c = c[c.gene.isin(genes)]
    if c.empty:
        return None
    c = c[c.groupby("gene")["p"].transform("min") < 0.05]
    if c.empty:
        return None
    top_cells = c.groupby("label")["p"].min().nsmallest(18).index
    c = c[c.label.isin(top_cells)]
    z = c.pivot_table(index="gene", columns="label", values="beta")
    p = c.pivot_table(index="gene", columns="label", values="p").reindex_like(z)
    text = np.where(p.values < 0.05, "•", "")
    lim = float(np.nanpercentile(np.abs(z.values), 95)) if z.size else 1
    fig = go.Figure(go.Heatmap(z=z.values, x=list(z.columns), y=list(z.index), colorscale=DIV, zmid=0, zmin=-lim, zmax=lim,
                               text=text, texttemplate="%{text}", colorbar=dict(title="log OR / SD", thickness=12),
                               hovertemplate="%{y} in %{x}<br>log OR %{z:.2f}<extra></extra>"))
    fig.update_xaxes(tickangle=-45, tickfont=dict(size=9))
    return _html(fig, max(380, 22 * len(z) + 220), "Gut cell types (IBDverse sc-eQTL): candidates with p < 0.05 in a cell type (• = p < 0.05)")


def tf_bars(n: int = 10) -> str:
    t = _p("tf_activity.parquet")
    t = t[(t.cohort == "GSE193677") & (t.kind == "tf")].sort_values("t")
    d = pd.concat([t.head(n), t.tail(n)])
    col = [CAT[7] if (x > 0 and q < 0.05) else CAT[0] if (x < 0 and q < 0.05) else GRAY for x, q in zip(d.t, d.padj)]
    fig = go.Figure(go.Bar(x=d.t, y=d.source, orientation="h", marker_color=col,
                           customdata=np.c_[d.cohorts_up, d.cohorts_down, d.padj],
                           hovertemplate="%{y}: t %{x:.1f}, up in %{customdata[0]} and down in %{customdata[1]} of 8 cohorts"
                                         ", q %{customdata[2]:.2g}<extra></extra>"))
    fig.update_xaxes(title="activity t, inflamed vs control (MSCCR; coloured: BH < 0.05)")
    return _html(fig, max(380, 20 * len(d) + 110), "Transcription-factor activity in inflamed mucosa (CollecTRI)", showlegend=False)


def pathway_heatmap() -> str:
    d = _p("pathway_activity.parquet")
    m = d.pivot_table(index="source", columns="cohort", values="t")
    m.columns = [COHORT.get(c, c) for c in m.columns]
    fig = go.Figure(go.Heatmap(z=m.values, x=list(m.columns), y=list(m.index), colorscale=DIV, zmid=0,
                               colorbar=dict(title="t", thickness=12), hovertemplate="%{y} · %{x}<br>t %{z:.1f}<extra></extra>"))
    return _html(fig, 460, "Pathway activity, UC vs control (PROGENy), eight cohorts")


def rank_intervals(name: str, n: int = 25) -> str | None:
    p = PROCESSED / f"sensitivity_{name}.parquet"
    if not p.exists():
        return None
    d = pd.read_parquet(p).head(n).iloc[::-1]
    fig = go.Figure(go.Scatter(
        x=d.rank_median, y=d.gene, mode="markers", marker=dict(color=CAT[0], size=8),
        error_x=dict(type="data", symmetric=False, array=d.rank_p95 - d.rank_median, arrayminus=d.rank_median - d.rank_p05,
                     color=CAT[0], thickness=1.4, width=3),
        customdata=np.c_[d.p_top10, d.p_top25, d.base_rank],
        hovertemplate="%{y}: median rank %{x}<br>P(top 10) %{customdata[0]:.2f} · P(top 25) %{customdata[1]:.2f}"
                      "<br>rank with configured weights %{customdata[2]}<extra></extra>"))
    fig.update_xaxes(title="rank under 1,000 random re-weightings (median, 5th–95th percentile)", type="log",
                     tickvals=[1, 2, 5, 10, 20, 50, 100, 200, 500], ticktext=["1", "2", "5", "10", "20", "50", "100", "200", "500"])
    return _html(fig, max(420, 19 * len(d) + 120), f"How stable is the {name} ranking?", showlegend=False)


def outcome_auroc(ev: dict) -> str:
    sets = [("clinical_targets", "UC clinical targets"), ("temporal_wide", "First UC trial ≥ 2016")]
    labels = ["convergence v1 (incl. literature + drugs)", "convergence without literature/tractability",
              "target profile (M14)", "predicted tractability alone", "expression alone", "biomarker profile (M14)"]
    labels = [l for l in labels if all(l in ev.get(s, {}) for s, _ in sets)]
    fig = go.Figure()
    for i, lab in enumerate(labels):
        y = [ev[s][lab]["auroc"] for s, _ in sets]
        lo = [ev[s][lab]["auroc"] - ev[s][lab]["auroc_ci95"][0] for s, _ in sets]
        hi = [ev[s][lab]["auroc_ci95"][1] - ev[s][lab]["auroc"] for s, _ in sets]
        fig.add_trace(go.Bar(x=[x for _, x in sets], y=y, name=lab.replace(" (M14)", ""), marker_color=CAT[i % len(CAT)],
                             error_y=dict(type="data", symmetric=False, array=hi, arrayminus=lo, color="#555")))
    fig.add_hline(y=0.5, line=dict(color=GRAY, dash="dot"))
    fig.update_yaxes(title="AUROC over the gene universe (95 % CI)", range=[0.3, 1])
    return _html(fig, 430, "Do the scores rank genes the clinic used? (0.5 = random)", barmode="group")


# ------------------------------------------------------------------ literature and drugs
def literature_charts() -> tuple[str, str]:
    rel = pd.read_parquet(LITERATURE / "relations.parquet")
    c = rel.relation.value_counts().sort_values()
    f1 = go.Figure(go.Bar(y=c.index, x=c.values, orientation="h", marker_color=CAT[0], marker_cornerradius=4,
                          hovertemplate="%{y}: %{x} relations<extra></extra>"))
    f1.update_xaxes(title="verbatim-validated relations")
    p = pd.read_parquet(LITERATURE / "publications.parquet")
    y = pd.to_numeric(p.year, errors="coerce").dropna().astype(int).value_counts().sort_index()
    f2 = go.Figure(go.Bar(x=y.index, y=y.values, marker_color=CAT[2], marker_cornerradius=3,
                          hovertemplate="%{x}: %{y} publications<extra></extra>"))
    f2.update_yaxes(title="publications in the local corpus")
    return (_html(f1, 330, "Relation types extracted and verbatim-validated", showlegend=False),
            _html(f2, 330, "Publication year of the local corpus", showlegend=False))


def drug_stage(ranking: list[dict]) -> str:
    stages = ["PHASE_1", "PHASE_2", "PHASE_3", "PHASE_4", "APPROVAL"]
    names = {"PHASE_1": "Phase 1", "PHASE_2": "Phase 2", "PHASE_3": "Phase 3", "PHASE_4": "Phase 4", "APPROVAL": "Approved"}
    rows = [(r["gene"], r["drugs"]) for r in ranking[:25] if r["drugs"]]
    fig = go.Figure()
    if rows:
        fig = go.Figure(go.Heatmap(z=[[sum(1 for d in drugs if d["max_stage"] == s) for s in stages] for _, drugs in rows],
                                   x=[names[s] for s in stages], y=[g for g, _ in rows], colorscale=SEQS, zmin=0,
                                   colorbar=dict(title="drugs", thickness=12),
                                   hovertemplate="%{y} · %{x}: %{z} drugs<extra></extra>"))
        fig.update_yaxes(autorange="reversed")
    return _html(fig, max(320, 22 * len(rows) + 110), "Known drugs by maximum clinical stage (top-25 candidates)")
