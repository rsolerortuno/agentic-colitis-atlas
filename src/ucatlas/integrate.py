"""Deterministic convergence scoring (no LLM). Writes runs/<run>/convergence.json.

Each stream is a component in [0, 1]; score = sum(stream_weight * component), weights in config.yaml. Every component
carries the evidence ids that justify it; only ids present in the graph are kept. Streams (v1):
- expression: meta-analysis change, min(1, |log2FC| / cap) when significant in the meta-analysis
- replication: share of the five cohorts with a significant change in the same direction as the meta-analysis
- inflammation: paired inflamed vs non-inflamed change within the same patient (local, not between people)
- genetics: Open Targets germline genetic association with the disease
- literature: max(min(1, validated relations / 3), Open Targets literature)
- tractability: 1 approved or clinical drug (any indication), 0.5 predicted tractable, 0 otherwise
"""
from __future__ import annotations

import json
from datetime import date

import numpy as np
import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .literature_entities import resolve_gene
from .paths import LITERATURE, PROCESSED
from .runs import run_dir

COHORTS = ["GSE193677", "GSE206285", "GSE92415", "GSE16879", "GSE73661"]


def compute() -> pd.DataFrame:
    cfg = load()["convergence"]
    meta = pd.read_parquet(PROCESSED / "meta.parquet").set_index("gene")
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    details = json.loads((PROCESSED / "external_targets.json").read_text())
    rel = pd.read_parquet(LITERATURE / "relations.parquet")

    genes = sorted(set(meta.index) | set(ot.index))
    t = pd.DataFrame(index=pd.Index(genes, name="gene"))
    t["meta_log2fc"] = meta["meta_log2fc"]
    t["meta_padj"] = meta["meta_padj"]
    sig = (t["meta_padj"] < 0.05).fillna(False)
    t["c_expression"] = np.where(sig, (t["meta_log2fc"].abs() / cfg["expression_cap"]).clip(upper=1), 0.0)
    t["direction"] = np.sign(t["meta_log2fc"]).fillna(0)
    t["n_cohorts_same_direction"] = np.where(t["direction"] > 0, meta["n_sig_up"].reindex(t.index),
                                             meta["n_sig_down"].reindex(t.index))
    t["n_cohorts_same_direction"] = t["n_cohorts_same_direction"].fillna(0)
    t["c_replication"] = np.where(sig, t["n_cohorts_same_direction"] / len(COHORTS), 0.0)

    pair = de[(de.cohort == "GSE193677") & (de.contrast == "inflamed_vs_noninflamed_paired")].set_index("gene")
    t["paired_log2fc"] = pair["log2fc"].reindex(t.index)
    t["paired_padj"] = pair["padj"].reindex(t.index)
    psig = (t["paired_padj"] < 0.05).fillna(False)
    t["c_inflammation"] = np.where(psig, (t["paired_log2fc"].abs() / cfg["paired_cap"]).clip(upper=1), 0.0)

    t["ot_score"] = ot["score"].reindex(t.index).fillna(0)
    t["ot_genetic"] = ot["dt_genetic_association"].reindex(t.index).fillna(0)
    t["ot_literature"] = ot["dt_literature"].reindex(t.index).fillna(0)
    t["c_genetics"] = t["ot_genetic"]

    lit_ids: dict[str, set[str]] = {}
    for r in rel.itertuples():
        for txt in (r.subject, r.object):
            g = resolve_gene(txt)
            if g:
                lit_ids.setdefault(g, set()).add(evidence_id("literature", r.relation_id))
    t["n_literature"] = [len(lit_ids.get(g, ())) for g in t.index]
    t["c_literature"] = np.maximum((t["n_literature"] / 3).clip(upper=1), t["ot_literature"])

    def tract(g):
        d = (details.get(g) or {}).get("open_targets")
        if not d:
            return 0.0
        return 1.0 if d["has_approved_or_clinical_drug"] else (0.5 if d["tractable_modalities"] else 0.0)

    t["c_tractability"] = [tract(g) for g in t.index]
    sw = cfg["stream_weights"]
    t["score"] = sum(sw[k] * t[f"c_{k}"] for k in sw)
    t["n_streams"] = sum((t[f"c_{k}"] > 0).astype(int) for k in sw)
    t["_lit_ids"] = [sorted(lit_ids.get(g, ())) for g in t.index]
    return t.reset_index().sort_values(["score", "n_streams"], ascending=False).reset_index(drop=True)


def _evidence(row, ot: pd.DataFrame, details: dict, store: GraphStore) -> dict:
    g = row["gene"]
    ot_ev = [evidence_id("open_targets_assoc", ot.loc[g, "ensembl"])] if g in ot.index else []
    ev = {"expression": [evidence_id("expression_meta", g)],
          "replication": [evidence_id("omics_de", c, "inflamed_vs_control" if c == "GSE193677" else "uc_vs_control", g)
                          for c in COHORTS],
          "inflammation": [evidence_id("omics_de", "GSE193677", "inflamed_vs_noninflamed_paired", g)],
          "genetics": ot_ev if row["ot_genetic"] > 0 else [],
          "literature": row["_lit_ids"] + (ot_ev if row["ot_literature"] > 0 else []),
          "tractability": []}
    d = (details.get(g) or {}).get("open_targets")
    if d:
        ev["tractability"] = [evidence_id("open_targets_tract", g)] + [
            evidence_id("open_targets_drug", x["chembl"], g) for x in d["drugs"]][:5]
    valid = set(store.check_evidence_ids([i for v in ev.values() for i in v])["valid"])
    return {k: [i for i in v if i in valid] for k, v in ev.items()}


def _r(x, nd=3):
    return None if x is None or pd.isna(x) else round(float(x), nd)


def _row(i, row, cfg, details, ot, store) -> dict:
    d = (details.get(row["gene"]) or {}).get("open_targets") or {}
    return {
        "rank": i + 1, "gene": row["gene"], "score": round(float(row["score"]), 4), "n_streams": int(row["n_streams"]),
        "components": {k: round(float(row[f"c_{k}"]), 3) for k in cfg["convergence"]["stream_weights"]},
        "meta_log2fc": _r(row["meta_log2fc"]), "meta_padj": None if pd.isna(row["meta_padj"]) else float(f"{row['meta_padj']:.3g}"),
        "cohorts_same_direction": int(row["n_cohorts_same_direction"]),
        "paired_log2fc": _r(row["paired_log2fc"]),
        "open_targets_score": _r(row["ot_score"], 4), "ot_genetic_association": _r(row["ot_genetic"], 4),
        "ot_literature": _r(row["ot_literature"], 4), "n_literature_relations": int(row["n_literature"]),
        "drugs": [{"name": x["name"], "max_stage": x["max_stage"]} for x in d.get("drugs", [])][:5],
        "evidence_ids": _evidence(row, ot, details, store),
    }


def write(n: int | None = None, day: str | None = None) -> dict:
    cfg = load()
    n = n or cfg["convergence"]["top_n_report"]
    t = compute()
    store = GraphStore()
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").set_index("symbol")
    details = json.loads((PROCESSED / "external_targets.json").read_text())
    ranking = [_row(i, row, cfg, details, ot, store) for i, row in t.head(n).iterrows()]
    profiles, weights = {}, cfg.get("profiles", {})
    if weights and (PROCESSED / "activity.parquet").exists():
        from .profiles import add
        t, sens = add(t)
        for name, w in weights.items():
            tp = t.sort_values(f"score_{name}", ascending=False).reset_index(drop=True)
            si = sens[name].set_index("gene")
            rows = []
            for i, row in tp.head(cfg["convergence"]["top_n_report"]).iterrows():
                r = _row(i, row, cfg, details, ot, store)
                r["score"] = round(float(row[f"score_{name}"]), 4)
                r["components"] = {k: round(float(row[f"c_{k}"]), 3) for k in w}
                if row["gene"] in si.index:
                    x = si.loc[row["gene"]]
                    r["rank_interval_90"] = {"rank_p05": int(x.rank_p05), "rank_median": float(x.rank_median),
                                             "rank_p95": int(x.rank_p95), "p_top10": round(float(x.p_top10), 3),
                                             "p_top25": round(float(x.p_top25), 3)}
                rows.append(r)
            profiles[name] = rows
    out = {"generated": date.today().isoformat(), "config": cfg["convergence"], "profile_weights": weights,
           "universe_size": int(len(t)), "ranking": ranking, "profiles": profiles}
    path = run_dir(day) / "convergence.json"
    path.write_text(json.dumps(out, indent=1))
    from .ranking_evidence import load as load_ranking_evidence
    ids = load_ranking_evidence(out, store)
    t.drop(columns="_lit_ids").to_parquet(PROCESSED / "convergence_table.parquet", index=False)
    return {"path": str(path), "top10": [r["gene"] for r in ranking[:10]], "universe": len(t),
            "ranking_evidence_ids": ids}
