"""M15: the pathway layer (new compared with the kidney atlas).

1. Pathway activity (PROGENy, 14 pathways, top 500 target genes) and transcription-factor activity (CollecTRI) per
   biopsy with a univariate linear model (decoupler ulm), on the patient-level discovery data (MSCCR) and the four
   microarray cohorts. Per cohort: inflamed / UC vs control (Welch), BH over sources; then the number of cohorts
   with a significant change in the same direction.
2. Cell-cell communication (OmniPath ligand-receptor pairs reported by >= 2 resources, receptor side annotated as
   a receptor in CellPhoneDB and ligand side not; complexes use their least-expressed subunit) on the
   single-cell data with neutrophils (Gubatan 2026, inflamed vs uninflamed UC): for every sender cell state,
   receiver cell state and pair, score = mean CP10K of the ligand in the sender x mean CP10K of the receptor in the
   receiver (CellPhoneDB-like product). Reported: the pairs whose score rises most in inflamed mucosa, with both
   genes expressed in >= 10 % of the cells involved.
Outputs: data/processed/pathway_activity.parquet, tf_activity.parquet, communication.parquet.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED, RAW

from .expression import COHORTS  # noqa: E402  (eight bulk cohorts)


def _cohort_matrix(cohort: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    s = s[(s.cohort == cohort) & (s.timepoint == "baseline") & ~s.in_gse16879]   # same samples as expression.py
    x = pd.read_parquet(PROCESSED / f"expr_{cohort}.parquet")
    if cohort == "GSE193677":
        from .expression import patient_means
        s = s[s.group.isin(["UC_inflamed", "healthy"])]
        pm = patient_means(x, s)
        return pm, [c for c in pm.columns if c.endswith("|UC_inflamed")], [c for c in pm.columns if c.endswith("|healthy")]
    s = s.set_index("geo")
    x = x[[c for c in x.columns if c in s.index]]
    return x, [c for c in x.columns if s.loc[c, "group"] == "UC"], [c for c in x.columns if s.loc[c, "group"] == "healthy"]


def _activity(net: pd.DataFrame, kind: str) -> pd.DataFrame:
    import decoupler as dc
    rows = []
    for cohort in COHORTS:
        x, a, b = _cohort_matrix(cohort)
        x = x.sub(x.mean(axis=1), axis=0)                      # gene-centred so activities compare samples
        est, pv = dc.mt.ulm(data=x.T, net=net, tmin=5, verbose=False)
        t, p = stats.ttest_ind(est.loc[a], est.loc[b], equal_var=False)
        d = pd.DataFrame({"source": est.columns, "diff": est.loc[a].mean().values - est.loc[b].mean().values,
                          "t": t, "p": p, "cohort": cohort, "n_disease": len(a), "n_control": len(b)})
        d["padj"] = multipletests(d["p"].fillna(1), method="fdr_bh")[1]
        rows.append(d)
    d = pd.concat(rows)
    sig = d[d.padj < 0.05]
    summary = d.groupby("source").agg(mean_t=("t", "mean")).join(
        sig.assign(up=sig.t > 0).groupby("source")["up"].agg(cohorts_up="sum", cohorts_sig="count"), how="left").fillna(0)
    summary["cohorts_down"] = summary["cohorts_sig"] - summary["cohorts_up"]
    d = d.merge(summary[["cohorts_up", "cohorts_down"]], left_on="source", right_index=True)
    d["kind"] = kind
    return d


def pathway_activity() -> pd.DataFrame:
    import decoupler as dc
    d = _activity(dc.op.progeny(organism="human", top=500), "pathway")
    d.to_parquet(PROCESSED / "pathway_activity.parquet", index=False)
    return d


def tf_activity() -> pd.DataFrame:
    import decoupler as dc
    d = _activity(dc.op.collectri(organism="human"), "tf")
    d.to_parquet(PROCESSED / "tf_activity.parquet", index=False)
    return d


def communication(min_pct: float = 10.0, min_cells: int = 50) -> pd.DataFrame:
    lr = pd.read_csv(RAW / "extra" / "omnipath_ligrec.tsv", sep="\t")
    # curated pairs only: >= 2 resources, the receptor side annotated as a receptor and the ligand side not
    import decoupler as dc
    ann = dc.op.resource("CellPhoneDB").drop_duplicates("genesymbol").set_index("genesymbol")
    is_rec = lambda g: all(bool(ann["receptor"].get(x, False)) for x in g.split("_"))
    lr["n_sources"] = lr["sources"].fillna("").str.split(";").map(len)
    lr = lr[lr["is_directed"] & (lr["n_sources"] >= 2)]
    lr = lr[lr["target_genesymbol"].map(is_rec) & ~lr["source_genesymbol"].map(is_rec)]
    lr = lr[["source_genesymbol", "target_genesymbol"]].drop_duplicates()
    sh = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet")
    sh = sh[sh["n_cells"] >= min_cells]
    mean = sh.pivot_table(index=["health", "cluster"], columns="gene", values="mean_cp10k", fill_value=0.0)
    pct = sh.pivot_table(index=["health", "cluster"], columns="gene", values="pct_expressing", fill_value=0.0)

    def level(genes: str, tab):  # complexes "A_B": least-expressed subunit
        parts = [g for g in genes.split("_")]
        if not all(g in tab.columns for g in parts):
            return None
        return tab[parts].min(axis=1)

    rows = []
    for lig, rec in lr.itertuples(index=False):
        ml, mr, pl, pr = level(lig, mean), level(rec, mean), level(lig, pct), level(rec, pct)
        if ml is None or mr is None:
            continue
        for h in ("inflamed", "noninflamed"):
            if h not in ml.index.get_level_values(0):
                continue
            L, R, PL, PR = ml.loc[h], mr.loc[h], pl.loc[h], pr.loc[h]
            snd = L[(PL >= min_pct) & (L > 0)]
            rcv = R[(PR >= min_pct) & (R > 0)]
            if snd.empty or rcv.empty:
                continue
            for s_, lv in snd.items():
                for r_, rv in rcv.items():
                    rows.append((lig, rec, s_, r_, h, float(lv * rv)))
    d = pd.DataFrame(rows, columns=["ligand", "receptor", "sender", "receiver", "health", "score"])
    w = d.pivot_table(index=["ligand", "receptor", "sender", "receiver"], columns="health", values="score", fill_value=0.0).reset_index()
    w["log2_change"] = np.log2((w.get("inflamed", 0) + 1e-3) / (w.get("noninflamed", 0) + 1e-3))
    w = w.sort_values("inflamed", ascending=False)
    w.to_parquet(PROCESSED / "communication.parquet", index=False)
    return w


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    today = date.today().isoformat()
    n = 0
    with store.batch():
        store.add_node("ANALYSIS:pathways", "Dataset", "Pathway (PROGENy) and TF (CollecTRI) activity, eight cohorts")
        store.add_node("ANALYSIS:communication", "Dataset", "Ligand-receptor communication, Gubatan 2026 single cell")
        for kind, f in (("pathway", "pathway_activity"), ("tf", "tf_activity")):
            d = pd.read_parquet(PROCESSED / f"{f}.parquet")
            for src, g in d.groupby("source"):
                per = {r.cohort: {"t": round(float(r.t), 2), "padj": float(f"{r.padj:.3g}")} for r in g.itertuples()}
                rec = {"name": src, "kind": kind, "cohorts_up": int(g["cohorts_up"].iloc[0]),
                       "cohorts_down": int(g["cohorts_down"].iloc[0]), "per_cohort": per,
                       "method": "decoupler ulm, inflamed/UC vs control, Welch t, BH"}
                nid = f"{'PATHWAY' if kind == 'pathway' else 'TF'}ACT:{src}"
                store.add_node(nid, "PathwayActivity", src, kind=kind)
                ev = store.add_evidence(evidence_id(f"{kind}_activity", src), f"{kind}_activity", "ANALYSIS:pathways", today, **rec)
                store.add_claim(nid, "activity_assessed_in", "ANALYSIS:pathways", ev)
                tf_node = store.resolve(src) if kind == "tf" else None
                if tf_node and tf_node["type"] == "Gene":
                    store.add_claim(tf_node["id"], "tf_activity_assessed_in", "ANALYSIS:pathways", ev)
                n += 1
        c = pd.read_parquet(PROCESSED / "communication.parquet")
        # deterministic selection, independent of what was loaded before: the 400 highest inflamed scores overall plus,
        # per ligand-receptor pair, the three sender-receiver combinations with the highest inflamed score
        c = c[(c.get("inflamed", 0) > 0)]
        c = pd.concat([c.head(400), c.groupby(["ligand", "receptor"], group_keys=False).head(3)]).drop_duplicates(
            ["ligand", "receptor", "sender", "receiver"])
        for r in c.itertuples():
            rec = {"ligand": r.ligand, "receptor": r.receptor, "sender": r.sender, "receiver": r.receiver,
                   "score_inflamed": round(float(r.inflamed), 3), "score_noninflamed": round(float(getattr(r, "noninflamed", 0.0)), 3),
                   "log2_change": round(float(r.log2_change), 3), "dataset": "Gubatan 2026 (SCP3755)"}
            eid = evidence_id("communication", r.ligand, r.receptor, r.sender, r.receiver)
            nid = f"COMM:{r.ligand}->{r.receptor}:{r.sender}->{r.receiver}"
            store.add_node(nid, "CellCommunication", f"{r.ligand} -> {r.receptor} ({r.sender} -> {r.receiver})")
            store.add_evidence(eid, "cell_communication", "ANALYSIS:communication", today, **rec)
            store.add_claim(nid, "observed_in", "ANALYSIS:communication", eid)
            for gname in (r.ligand.split("_")[0], r.receptor.split("_")[0]):
                g = store.resolve(gname)
                if g and g["type"] == "Gene":
                    store.add_claim(g["id"], "participates_in", nid, eid)
            n += 1
    return n


if __name__ == "__main__":
    p = pathway_activity()
    print(p.drop_duplicates("source")[["source", "cohorts_up", "cohorts_down"]].sort_values("cohorts_up", ascending=False).to_string())
    t = tf_activity()
    print(t.drop_duplicates("source").sort_values(["cohorts_up"], ascending=False)[["source", "cohorts_up"]].head(15).to_string())
    c = communication()
    print(c.head(20)[["ligand", "receptor", "sender", "receiver", "inflamed", "log2_change"]].to_string())
    print("records:", load())
