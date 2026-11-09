"""M14: disease activity (replaces survival in the kidney atlas). Per gene, Spearman correlation of expression with
the Mayo endoscopic score across UC colon biopsies of the discovery cohort (MSCCR, one value per patient: the mean
of the patient's biopsies and the patient's maximum score), replicated in GSE73661 (baseline Mayo endoscopic
subscore). Output: data/processed/activity.parquet."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED


def _spearman(x: pd.DataFrame, score: pd.Series) -> pd.DataFrame:
    xr = x.rank(axis=1)
    sr = score.rank()
    xr = xr.sub(xr.mean(axis=1), axis=0)
    sr = sr - sr.mean()
    rho = (xr @ sr) / np.sqrt((xr ** 2).sum(axis=1) * (sr ** 2).sum())
    n = len(score)
    t = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return pd.DataFrame({"rho": rho, "p": p, "n": n})


def run() -> pd.DataFrame:
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    out = []
    # MSCCR: UC colon biopsies with a Mayo endoscopic score
    m = s[(s.cohort == "GSE193677") & s.group.isin(["UC_inflamed", "UC_noninflamed"]) & s.mayo.notna()]
    x = pd.read_parquet(PROCESSED / "expr_GSE193677.parquet")[m["sample"].tolist()]
    xp = x.T.groupby(m.set_index("sample")["patient"].reindex(x.columns).values).mean().T
    score = m.groupby("patient")["mayo"].max().reindex(xp.columns)
    out.append(_spearman(xp, score).assign(cohort="GSE193677"))
    # GSE73661: UC at baseline with the Mayo endoscopic subscore
    v = s[(s.cohort == "GSE73661") & (s.group == "UC") & (s.timepoint == "baseline") & s.mayo.notna()]
    xv = pd.read_parquet(PROCESSED / "expr_GSE73661.parquet")
    xv = xv[[c for c in v["geo"] if c in xv.columns]]
    out.append(_spearman(xv, v.set_index("geo")["mayo"].reindex(xv.columns)).assign(cohort="GSE73661"))
    d = pd.concat(out).rename_axis("gene").reset_index()
    d["padj"] = d.groupby("cohort")["p"].transform(lambda p: multipletests(p.fillna(1), method="fdr_bh")[1])
    d.to_parquet(PROCESSED / "activity.parquet", index=False)
    return d


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "activity.parquet").dropna(subset=["rho"])
    n = 0
    with store.batch():
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            ev = store.add_evidence(evidence_id("disease_activity", r.cohort, r.gene), "disease_activity",
                                    f"DATASET:{r.cohort}", date.today().isoformat(), gene=r.gene, cohort=r.cohort,
                                    score="Mayo endoscopic score", spearman_rho=round(float(r.rho), 3),
                                    padj=float(f"{r.padj:.3g}"), n=int(r.n))
            store.add_claim(node["id"], "tracks_activity_in", f"DATASET:{r.cohort}", ev)
            n += 1
    return n


if __name__ == "__main__":
    d = run()
    print(d.groupby("cohort").apply(lambda x: (x.padj < 0.05).sum()).to_dict(), "records:", load())
    print(d[d.gene.isin(["LCN2", "S100A8", "CXCL8", "AQP8", "HMGCS2", "DUOX2"])].round(3).to_string())
