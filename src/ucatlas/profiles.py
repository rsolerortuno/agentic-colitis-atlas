"""M14: question-specific profiles without literature or drug history, and their weight sensitivity.

Components (all in [0, 1]):
- expression_up: c_expression when the meta-analysis change is positive (target profile only counts increases)
- tissue_change: composition-retained share of the change (0 when the verdict is 'composition')
- genetics: Open Targets germline genetic association
- tractability_predicted: 1 if Open Targets predicts small-molecule or antibody tractability, else 0
- activity: max |Spearman rho| with the Mayo endoscopic score over the cohorts where it is significant, scaled by 0.5
- expression_any, replication, inflammation: as in convergence v1
Weights in config.yaml (profiles). Sensitivity: 1,000 reweightings, log-uniform factor in [1/2, 2] per weight.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .config import load
from .paths import PROCESSED


def tractability_table() -> pd.DataFrame:
    """Open Targets 26.09 tractability for every target (same buckets as the kidney atlas): clinical precedence,
    predicted SM / AB tractability (strong buckets: SM high-quality ligand or pocket, AB high-confidence surface or
    secreted location) and the same without the ligand bucket (medicinal-chemistry history)."""
    from .paths import RAW
    parts = sorted((RAW / "extra").glob("ot_target_tractability_*.parquet"))
    t = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    t = t[t["value"]].copy()
    clinical = t["category"].isin(["Approved Drug", "Advanced Clinical", "Phase 1 Clinical"])
    pred = (((t["modality"] == "SM") & t["category"].isin(["High-Quality Ligand", "High-Quality Pocket"]))
            | ((t["modality"] == "AB") & (t["category"] == "UniProt loc high conf")))
    pred_nl = (((t["modality"] == "SM") & t["category"].eq("High-Quality Pocket"))
               | ((t["modality"] == "AB") & (t["category"] == "UniProt loc high conf")))
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    t["gene"] = t["targetId"].map(dict(zip(hg["ensembl_gene_id"], hg["symbol"])))
    t = t.dropna(subset=["gene"])
    genes = sorted(set(t["gene"]))
    out = pd.DataFrame(index=pd.Index(genes, name="gene"))
    out["clinical"] = out.index.isin(t.loc[clinical, "gene"])
    out["predicted"] = out.index.isin(t.loc[pred, "gene"])
    out["predicted_no_ligand"] = out.index.isin(t.loc[pred_nl, "gene"])
    return out


def components(t: pd.DataFrame) -> pd.DataFrame:
    t = t.copy()
    t["c_expression_up"] = np.where(t["meta_log2fc"] > 0, t["c_expression"], 0.0)
    t["c_expression_any"] = t["c_expression"]
    comp = pd.read_parquet(PROCESSED / "composition_effect.parquet").set_index("gene")
    ret = comp["retained"].where(comp["verdict"] != "composition", 0.0).clip(0, 1)
    t["c_tissue_change"] = t["gene"].map(ret).fillna(0.0)
    act = pd.read_parquet(PROCESSED / "activity.parquet")
    act = act[act["padj"] < 0.05].groupby("gene")["rho"].apply(lambda r: r.abs().max())
    t["c_activity"] = (t["gene"].map(act).fillna(0.0) / 0.5).clip(upper=1)
    # every gene from the Open Targets bulk table (v1 of the profiles used only the ~1,100 enriched genes, which
    # favoured genes already associated with UC; corrected in the review of 2026-11-10, evals/history.md)
    tr = tractability_table()
    t["c_tractability_predicted"] = t["gene"].map(tr["predicted"]).fillna(False).astype(float)
    t["c_tractability_predicted_no_ligand"] = t["gene"].map(tr["predicted_no_ligand"]).fillna(False).astype(float)
    return t


def score(t: pd.DataFrame, weights: dict) -> pd.Series:
    return sum(w * t[f"c_{k}"] for k, w in weights.items())


def sensitivity(t: pd.DataFrame, weights: dict, n: int = 1000, seed: int = 11) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keys = list(weights)
    C = t[[f"c_{k}" for k in keys]].to_numpy(float)
    w0 = np.array([weights[k] for k in keys])
    base = C @ w0
    top = np.argsort(-base)[:200]
    ranks = np.empty((n, len(top)), int)
    for i in range(n):
        w = w0 * np.exp(rng.uniform(np.log(0.5), np.log(2), len(keys)))
        s = C @ w
        order = np.argsort(-s)
        r = np.empty(len(s), int)
        r[order] = np.arange(1, len(s) + 1)
        ranks[i] = r[top]
    return pd.DataFrame({"gene": t["gene"].to_numpy()[top], "rank_p05": np.percentile(ranks, 5, axis=0).astype(int),
                         "rank_median": np.median(ranks, axis=0), "rank_p95": np.percentile(ranks, 95, axis=0).astype(int),
                         "p_top10": (ranks <= 10).mean(axis=0), "p_top25": (ranks <= 25).mean(axis=0)})


def add(t: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    prof = load().get("profiles", {})
    t = components(t)
    sens = {}
    for name, w in prof.items():
        t[f"score_{name}"] = score(t, w)
        sens[name] = sensitivity(t, w)
    return t, sens
