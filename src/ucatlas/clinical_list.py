"""Targets in clinical validation that the omics support, with the approved ones for comparison (added 2026-11-27).

The four other lists answer one interest each (first-in-class, well known without a drug, repurposing, genetics
first). This one answers: which genes already have a drug tried in UC or IBD, and how strongly do the data of this
atlas back them? A gene enters when a drug acting on it reached clinical trials in UC or IBD (drug_indications.parquet:
UC or Crohn's/IBD at any stage, or a trial recorded outside Open Targets) and, unless the drug is approved there:

- its change is significant in the eight-cohort meta-analysis (BH < 0.05) and goes the same way in at least
  MIN_COHORTS cohorts;
- it is well studied: more than MIN_PAPERS Europe PMC records name it with ulcerative colitis.

Genes with a drug approved in UC or IBD are all listed, whatever their data, and marked approved: they show where
validated targets fall on the same scores. Rows are ordered by the convergence score; the data-only rank (expression,
replication, inflammation and genetics, without the literature and drug-history streams, which favour approved
targets) is shown next to it. A row whose only UC or IBD drug is doxycycline is flagged: that trial tested an
antibiotic, not a drug made against the gene.
"""
from __future__ import annotations

import pandas as pd

from .paths import PROCESSED

MIN_COHORTS = 5
MIN_PAPERS = 100
OUT = PROCESSED / "clinical_list.parquet"
STAGE_ORDER = ["APPROVAL", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1", "EARLY_PHASE_1"]
NOT_DIRECTED = {"doxycycline"}


def _data_only(c: pd.DataFrame) -> pd.Series:
    from .config import load
    w = load()["convergence"]["stream_weights"]
    return (w["expression"] * c["c_expression"] + w["replication"] * c["c_replication"]
            + w["inflammation"] * c["c_inflammation"] + w["genetics"] * c["c_genetics"])


def build() -> pd.DataFrame:
    from .external import DISEASE_TERMS, europepmc_count
    from .graph import evidence_id
    di = pd.read_parquet(PROCESSED / "drug_indications.parquet")
    di = di[di.category.isin(["UC", "Crohn's / IBD", "tried in UC/IBD (outside Open Targets)"])].copy()
    c = pd.read_parquet(PROCESSED / "convergence_table.parquet").set_index("gene")
    c["data_only"] = _data_only(c)
    c["rank"] = c["score"].rank(ascending=False, method="min")
    c["rank_data_only"] = c["data_only"].rank(ascending=False, method="min")
    rows = []
    for r in di.itertuples():
        if r.category == "UC":
            setting, stage, drugs = "UC", r.uc_stage, list(r.uc_drugs)
        elif r.category.startswith("Crohn"):
            setting, stage, drugs = "Crohn's / IBD", r.ibd_stage, list(r.ibd_drugs)
        else:
            setting, stage, drugs = "UC / IBD (outside Open Targets)", None, []
        approved = stage == "APPROVAL"
        x = c.loc[r.gene] if r.gene in c.index else None
        sig = x is not None and pd.notna(x["meta_padj"]) and x["meta_padj"] < 0.05
        coh = int(x["n_cohorts_same_direction"]) if x is not None and pd.notna(x["n_cohorts_same_direction"]) else 0
        if not approved and not (sig and coh >= MIN_COHORTS):
            continue
        papers = europepmc_count(f'"{r.gene}" AND {DISEASE_TERMS}')
        if not approved and papers <= MIN_PAPERS:
            continue
        names = [str(d) for d in drugs]
        rows.append({
            "gene": r.gene, "approved": approved, "setting": setting, "stage": stage or "tried",
            "drugs": names[:4], "n_drugs": len(names),
            "not_directed": bool(names) and {n.lower() for n in names} <= NOT_DIRECTED,
            "meta_log2fc": None if x is None or pd.isna(x["meta_log2fc"]) else round(float(x["meta_log2fc"]), 2),
            "meta_padj": None if x is None or pd.isna(x["meta_padj"]) else float(x["meta_padj"]),
            "cohorts_same_direction": coh, "uc_papers": int(papers),
            "score": None if x is None else round(float(x["score"]), 3),
            "rank": None if x is None else int(x["rank"]),
            "rank_data_only": None if x is None else int(x["rank_data_only"]),
            "universe": len(c),
            "evidence_ids": [evidence_id("clinical_indications", r.gene), evidence_id("expression_meta", r.gene)]})
    d = pd.DataFrame(rows).sort_values("score", ascending=False, na_position="last").reset_index(drop=True)
    d.to_parquet(OUT, index=False)
    return d


def view() -> dict:
    if not OUT.exists():
        return {}
    d = pd.read_parquet(OUT)
    rows = []
    for r in d.to_dict("records"):
        r["drugs"] = [x.upper() if any(ch.isdigit() for ch in x) else x for x in r["drugs"]]
        r["evidence_ids"] = list(r["evidence_ids"])
        rows.append(r)
    appr = [r for r in rows if r["approved"]]
    return {"rows": rows, "n": len(rows), "approved": len(appr),
            "approved_genes": [r["gene"] for r in appr],
            "trial": len(rows) - len(appr), "min_cohorts": MIN_COHORTS, "min_papers": MIN_PAPERS,
            "universe": int(d["universe"].iloc[0]) if len(d) else None,
            "approved_median_rank": int(pd.Series([r["rank"] for r in appr if r["rank"]]).median()) if appr else None,
            "approved_median_rank_data_only": int(pd.Series([r["rank_data_only"] for r in appr
                                                             if r["rank_data_only"]]).median()) if appr else None}


if __name__ == "__main__":
    d = build()
    print(len(d), "genes;", int(d.approved.sum()), "approved")
    print(d[["gene", "approved", "setting", "stage", "meta_log2fc", "cohorts_same_direction", "uc_papers", "rank",
             "rank_data_only"]].head(40).to_string())
