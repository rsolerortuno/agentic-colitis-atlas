"""Shortlist for experimental validation (added 2026-11-10): one row per first-in-class gene, with a traffic light per
criterion and a proposed first experiment. Everything comes from tables already in data/processed; nothing is
re-estimated, and the rules below were written before the table was looked at.

Criteria (green / amber / red; grey = not measured):
- reproducible: cohorts with a significant increase (green >= 6 of 8, amber 4-5, red < 4), and the covariate-adjusted
  MSCCR change holds (if it does not, at most amber);
- protein: the protein also rises in the independent proteomics study (green: replicates; red: measured, does not
  replicate; grey: not quantified);
- tissue change: the composition verdict the discovery screen already applied (green 'intrinsic': the change is
  of the tissue; amber 'partly_composition': part of it is cells arriving), with the share of the change kept;
- cell source: share of the gene's transcripts in its dominant lineage in inflamed single cells (Smillie 2019, the
  lineage the discovery screen uses; green >= 0.5, amber 0.3-0.5, red < 0.3); a lineage not replicated in the second
  atlas (Kong 2023) caps it at amber; the Visium in-situ expresser is shown when available;
- causal genetics: a significant Mendelian-randomisation estimate in plasma or colon in the direction that supports
  inhibition (genetically higher level raises risk) is green, the opposite direction red, tested and null amber, no
  instrument grey;
- druggable: an antibody route (surface or secreted, Open Targets buckets) or a pocket score >= 0.6 is green, a
  predicted tractability bucket only amber;
- safety: highest vital-organ RNA, loss-of-function constraint (LOEUF < 0.35) and known safety events; green with
  none, amber with one flag, red with constraint plus another flag;
- mouse: whether the gene has ever been perturbed in mouse colitis (information, not a score).
Proposed experiment: by dominant lineage (epithelial: patient-derived colon organoids; stromal: primary colonic
fibroblasts; immune: lamina propria immune cells) and by route (antibody, small molecule or genetic), followed by a
mouse colitis model when none exists.
Outputs: data/processed/validation_shortlist.parquet.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .paths import PROCESSED

ORDER = {"green": 2, "amber": 1, "grey": 0, "red": -1}


def _lineage_share() -> dict[str, tuple[str, float]]:
    # the same source as the lineage of the discovery screen (Smillie 2019, inflamed), so the two always agree
    from .sc_replication import shares_2019
    lin = shares_2019("inflamed").sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    return {g: (r.lineage, float(r.share)) for g, r in lin.iterrows()}


def _experiment(lineage: str, surface: bool, pocket: float | None, mouse_papers: int) -> str:
    if lineage == "epithelial":
        system = "knock-down or CRISPR knockout in patient-derived colon organoids stimulated with TNF and IFN-γ " \
                 "(read-outs: barrier integrity, chemokine release, epithelial death)"
    elif lineage in ("stromal", "endothelial"):
        system = "silencing in primary colonic fibroblasts from inflamed UC stimulated with OSM and IL-1β " \
                 "(read-outs: IL-6, CXCL1, MMP release; neutrophil attraction in co-culture)"
    else:
        system = "knock-down in lamina propria immune cells or blood-derived equivalents " \
                 "(read-outs: activation and cytokine profile)"
    route = ("a blocking antibody can be tested directly" if surface else
             "small-molecule screening is possible (pocket)" if (pocket or 0) >= 0.6 else "genetic tools first")
    mouse = " Then a knockout or blockade in DSS or T-cell transfer colitis (never tested in mouse colitis)." \
        if mouse_papers == 0 else f" Mouse colitis data exist ({mouse_papers} papers): check them first."
    return f"{system}; {route}.{mouse}"


def run() -> dict:
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    fic = d[d.status == "first-in-class opportunity"].copy()
    adj = pd.read_parquet(PROCESSED / "de_adjusted.parquet").set_index("gene")
    prot = pd.read_parquet(PROCESSED / "proteomics.parquet").set_index("gene")
    lin = _lineage_share()
    rep_p = PROCESSED / "sc_replication.parquet"
    screp = pd.read_parquet(rep_p).set_index("gene")["same_lineage"].to_dict() if rep_p.exists() else {}
    vis = {t["gene"]: t for t in json.loads((PROCESSED / "spatial_visium_summary.json").read_text())["targets"]}
    mr_rows = {}
    for f_, lab in (("mr_pqtl.parquet", "plasma"), ("mr_colon.parquet", "colon")):
        r = pd.read_parquet(PROCESSED / f_)
        if "layer" in r:
            r = r[r.layer == "colon_tissue"]
        for x in r.itertuples():
            mr_rows.setdefault(x.gene, []).append((lab, float(x.beta), float(x.q)))
    bio = pd.read_parquet(PROCESSED / "biologics.parquet").set_index("gene")
    pk = pd.read_parquet(PROCESSED / "structure_pockets.parquet")
    pocket = pk[pk["rank"] == 1].set_index("gene")["pocket_score"].to_dict()
    saf = pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene")
    ml = pd.read_parquet(PROCESSED / "mouse_literature.parquet").set_index("gene")["papers"].to_dict()
    rows = []
    for x in fic.itertuples():
        g = x.gene
        c = {}
        rep = "green" if x.cohorts_up >= 6 else "amber" if x.cohorts_up >= 4 else "red"
        holds = bool(adj.at[g, "holds"]) if g in adj.index else None
        if holds is False and rep == "green":
            rep = "amber"
        c["reproducible"] = (rep, f"up in {int(x.cohorts_up)}/8 cohorts" + ("" if holds is None else
                                                                             f"; adjusted: {'holds' if holds else 'does not hold'}"))
        if g in prot.index and pd.notna(prot.at[g, "padj"]):
            r_ = prot.at[g, "replicates"]
            c["protein"] = ("green" if r_ == "yes" else "red" if r_ == "no" else "amber",
                            f"log2 ratio {prot.at[g, 'log2_ratio']:.2f}, padj {prot.at[g, 'padj']:.1e}")
        else:
            c["protein"] = ("grey", "not quantified")
        ret = float(x.retained) if pd.notna(x.retained) else float("nan")
        c["tissue"] = ("green" if x.composition == "intrinsic" else "amber" if x.composition == "partly_composition"
                       else "red", f"{x.composition.replace('_', ' ')}; {ret:.0%} of the change kept after adjustment")
        lname, lshare = lin.get(g, (x.lineage, np.nan))
        vtxt = f"; in situ {vis[g]['top_in_situ_expresser'].replace('_', ' ')}" if g in vis and vis[g].get("top_in_situ_expresser") else ""
        col = "grey" if np.isnan(lshare) else "green" if lshare >= 0.5 else "amber" if lshare >= 0.3 else "red"
        rtxt = "" if g not in screp else "; same lineage in the second atlas" if screp[g] else "; not replicated in the second atlas"
        if col == "green" and screp.get(g) is False:
            col = "amber"
        c["cell"] = (col, f"{lname} {lshare:.0%} of transcripts{rtxt}{vtxt}" if not np.isnan(lshare) else "not in the single-cell atlas")
        mrs = mr_rows.get(g, [])
        sig = [m for m in mrs if m[2] < 0.05]
        if sig:
            up = any(m[1] > 0 for m in sig)
            c["genetics"] = ("green" if up else "red",
                             "; ".join(f"{m[0]} OR {np.exp(m[1]):.2f} (q {m[2]:.1e})" for m in sig))
        elif mrs:
            c["genetics"] = ("amber", "tested, no causal effect")
        else:
            c["genetics"] = ("grey", "no genetic instrument")
        surface = g in bio.index and bio.at[g, "c_surface"] >= 0.7
        pscore = pocket.get(g)
        mods = list(x.tractable)
        c["druggable"] = ("green" if surface or (pscore or 0) >= 0.6 else "amber" if mods else "red",
                          ", ".join(filter(None, ["surface or secreted" if surface else "",
                                                  f"pocket {pscore:.2f}" if pscore is not None else "",
                                                  "predicted " + "/".join(mods) if mods else ""])) or "no route")
        flags = []
        if pd.notna(x.vital_max_ntpm) and x.vital_max_ntpm > 50:
            flags.append(f"vital organs {x.vital_max_ntpm:.0f} nTPM")
        if g in saf.index:
            if saf.at[g, "constrained"] is True:
                flags.append(f"constrained (LOEUF {saf.at[g, 'loeuf']:.2f})")
            if saf.at[g, "n_safety_events"] > 0:
                flags.append("known safety events")
        c["safety"] = ("green" if not flags else "red" if len(flags) >= 2 and any("constrained" in f for f in flags)
                       else "amber", "; ".join(flags) or f"vital organs ≤ {x.vital_max_ntpm:.0f} nTPM, not constrained")
        papers = int(ml.get(g, 0))
        c["mouse"] = ("grey", "never tested in mouse colitis" if papers == 0 else f"{papers} mouse colitis papers")
        score = sum(ORDER[v[0]] for k, v in c.items() if k != "mouse")
        rows.append({"gene": g, "lineage": x.lineage, "track": x.track, "meta_log2fc": round(float(x.meta_log2fc), 2),
                     "uc_papers": int(x.uc_papers), "criteria": json.dumps(c, ensure_ascii=False),
                     "greens": sum(v[0] == "green" for k, v in c.items() if k != "mouse"),
                     "reds": sum(v[0] == "red" for k, v in c.items() if k != "mouse"), "score": score,
                     "experiment": _experiment(x.lineage, surface, pscore, papers)})
    out = pd.DataFrame(rows).sort_values(["score", "greens", "meta_log2fc"], ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out.to_parquet(PROCESSED / "validation_shortlist.parquet", index=False)
    return {"genes": len(out), "top5": out.gene.head(5).tolist(),
            "no_red": int((out.reds == 0).sum()), "with_genetics": int(out.criteria.str.contains('"genetics": \\["green"').sum())}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, ensure_ascii=False))
    d = pd.read_parquet(PROCESSED / "validation_shortlist.parquet")
    for r in d.head(12).itertuples():
        c = json.loads(r.criteria)
        print(r.rank, r.gene, r.lineage, r.score, {k: v[0] for k, v in c.items()})
