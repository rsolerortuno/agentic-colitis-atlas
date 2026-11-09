"""M8: evaluation against outcomes the scores do not see (deterministic, no LLM). Decided before measuring.

Gold set: UC clinical targets from Open Targets 26.09 clinical precedence: targets of drugs tested in ulcerative
colitis in phase 2 or later, with >= 2 trials or an approval. Drugs whose mechanism maps to more than 5 targets are
removed so one drug does not contribute a block of genes (same rule as the earlier oncology atlas).
Temporal subset: targets whose first UC trial (any phase) started on or after 2018-01-01. This is retrospective:
the features are not frozen at the cutoff.

Scores compared: the convergence score (includes literature and drug history, so partly circular), the same score
without literature and tractability ("data only"), expression alone, genetics alone and tractability alone.
Metrics: AUROC over the gene universe with a bootstrap 95 % CI over gold genes, recall@100 / @500, median rank,
fold enrichment at 100.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from .config import load
from .paths import PROCESSED, RAW

UC = re.compile(r"ulcerative colitis", re.I)
ADVANCED = {"PHASE_2", "PHASE_2_3", "PHASE_3", "PHASE_4", "APPROVAL"}
CUTOFF = "2018-01-01"
MAX_TARGETS_PER_DRUG = 5


def clinical_targets() -> pd.DataFrame:
    cp = pd.read_parquet(RAW / "extra" / "ot_clinical_precedence.parquet",
                         columns=["targetId", "drugId", "clinicalReportId", "clinicalStage", "studyStartDate",
                                  "diseaseFromSource"])
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    symbols = dict(zip(hg["ensembl_gene_id"], hg["symbol"]))
    n_targets = cp.groupby("drugId")["targetId"].nunique()
    r = cp[cp["diseaseFromSource"].fillna("").str.contains(UC)].copy()
    r["gene"] = r["targetId"].map(symbols)
    r["start"] = pd.to_datetime(r["studyStartDate"], errors="coerce")
    dates = r.dropna(subset=["gene"]).groupby("gene").agg(first_trial=("start", "min"),
                                                          date_incomplete=("start", lambda s: s.isna().any()))
    eligible = r[r["clinicalStage"].isin(ADVANCED) & (r["drugId"].map(n_targets) <= MAX_TARGETS_PER_DRUG)]
    g = eligible.dropna(subset=["gene"]).groupby("gene").agg(
        trials=("clinicalReportId", "nunique"), drugs=("drugId", "nunique"),
        approved=("clinicalStage", lambda s: "APPROVAL" in set(s)))
    g = g[(g["trials"] >= 2) | g["approved"]].join(dates)
    g["holdout"] = (g["first_trial"] >= pd.Timestamp(CUTOFF)) & ~g["date_incomplete"]
    return g.reset_index()


def auroc(score: np.ndarray, gold: np.ndarray) -> float:
    pos, neg = score[gold], score[~gold]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def metrics(t: pd.DataFrame, col: str, gold: set[str], exclude: set[str] | None = None, seed: int = 3) -> dict:
    d = t[~t["gene"].isin(exclude or set())]
    s = d[col].to_numpy(float)
    g = d["gene"].isin(gold).to_numpy()
    order = np.lexsort((np.random.default_rng(seed + 1).random(len(s)), -s))   # random order within ties
    rank = np.empty(len(s), int)
    rank[order] = np.arange(1, len(s) + 1)
    rng = np.random.default_rng(seed)
    pos_idx = np.nonzero(g)[0]
    boots = [auroc(np.concatenate([s[smp], s[~g]]), np.concatenate([np.ones(len(smp), bool), np.zeros((~g).sum(), bool)]))
             for smp in (rng.choice(pos_idx, size=len(pos_idx), replace=True) for _ in range(500))] if len(pos_idx) else [np.nan]
    n = int(g.sum())
    rec100 = float((rank[g] <= 100).mean()) if n else float("nan")
    return {"n_gold": n, "universe": int(len(s)), "auroc": round(auroc(s, g), 3),
            "auroc_ci95": [round(float(np.nanpercentile(boots, 2.5)), 3), round(float(np.nanpercentile(boots, 97.5)), 3)],
            "recall_at_100": round(rec100, 3), "recall_at_500": round(float((rank[g] <= 500).mean()), 3) if n else None,
            "median_rank": int(np.median(rank[g])) if n else None,
            "fold_enrichment_at_100": round(rec100 / (100 / len(s)), 1) if n else None,
            "top_hits": d.iloc[order[:100]].loc[lambda x: x["gene"].isin(gold), "gene"].tolist()[:15]}


def run() -> dict:
    from .integrate import compute
    t = compute()
    sw = load()["convergence"]["stream_weights"]
    t["score_data_only"] = t["score"] - sw["literature"] * t["c_literature"] - sw["tractability"] * t["c_tractability"]
    noise = np.random.default_rng(0).random(len(t)) * 1e-9
    for k in ("expression", "genetics", "tractability"):
        t[f"only_{k}"] = t[f"c_{k}"] + noise
    if load().get("profiles") and (PROCESSED / "activity.parquet").exists():
        from .profiles import components
        t = components(t)
        prof = load()["profiles"]
        for name, w in prof.items():
            t[f"score_{name}"] = sum(wt * t[f"c_{k}"] for k, wt in w.items())
        tw = dict(prof["target"])
        t["score_target_no_tractability"] = t["score_target"] - tw.get("tractability_predicted", 0) * t["c_tractability_predicted"]
        t["only_tractability_predicted"] = t["c_tractability_predicted"] + noise
    clin = clinical_targets()
    allc, hold = set(clin["gene"]), set(clin.loc[clin["holdout"], "gene"])
    pre = allc - hold
    scores = {"convergence v1 (incl. literature + drugs)": "score", "convergence without literature/tractability": "score_data_only",
              "expression alone": "only_expression", "genetics alone": "only_genetics", "tractability alone": "only_tractability"}
    if "score_target" in t:
        scores.update({"target profile (M14)": "score_target", "target profile without tractability": "score_target_no_tractability",
                       "predicted tractability alone": "only_tractability_predicted", "biomarker profile (M14)": "score_biomarker"})
    out = {"gold_sets": {"uc_clinical_targets": sorted(allc), "temporal_holdout_targets": sorted(hold), "cutoff": CUTOFF},
           "clinical_targets": {k: metrics(t, c, allc) for k, c in scores.items()},
           "temporal_holdout": {k: metrics(t, c, hold, exclude=pre) for k, c in scores.items()},
           "clinical_table": clin.assign(first_trial=clin["first_trial"].dt.date.astype(str)).to_dict("records")}
    (PROCESSED.parents[1] / "evals").mkdir(exist_ok=True)
    (PROCESSED.parents[1] / "evals" / "results.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    return out


if __name__ == "__main__":
    r = run()
    print(len(r["gold_sets"]["uc_clinical_targets"]), "clinical targets;", len(r["gold_sets"]["temporal_holdout_targets"]), "after", CUTOFF)
    for key in ("clinical_targets", "temporal_holdout"):
        print(key)
        for k, v in r[key].items():
            print(f"  {k}: AUROC {v['auroc']} {v['auroc_ci95']} recall@100 {v['recall_at_100']} median rank {v['median_rank']}")
