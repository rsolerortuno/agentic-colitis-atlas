"""M16: response to four biologics, from baseline biopsies (deterministic, no LLM).

Cohorts and endpoints (data/README.md): UNIFI ustekinumab vs placebo (mucosal healing at week 8), PURSUIT golimumab vs
placebo (clinical response at week 6), GSE16879 infliximab (endoscopic and histological response), GSE73661
vedolizumab (earliest post-treatment response label). GSE73661's infliximab patients probably overlap GSE16879 and
are not used.

1. Per drug arm and gene: non-responders vs responders at baseline (Welch t on log2 expression, BH). Positive
   difference = higher in non-responders.
2. Drug-specific vs general: logistic regression response ~ expression * arm in UNIFI and PURSUIT (drug vs placebo),
   for a candidate gene set; the interaction term asks whether a gene predicts response to the drug more than to
   placebo. No gene is expected to survive correction; reported as association.
3. Pathway level: PROGENy and CollecTRI activities, non-responders vs responders per drug arm.
4. Positive control: the infliximab non-response genes of Arijs et al. 2009 (TNFRSF11B, STC1, PTGS2, IL13RA2, IL11)
   should be higher in non-responders of GSE16879.
Outputs: data/processed/response_genes.parquet, response_interaction.parquet, response_pathways.parquet.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .paths import PROCESSED

ARMS = {"ustekinumab": ("GSE206285", "ustekinumab"), "golimumab": ("GSE92415", "golimumab"),
        "infliximab": ("GSE16879", "infliximab"), "vedolizumab": ("GSE73661", "vedolizumab"),
        "placebo (UNIFI)": ("GSE206285", "placebo"), "placebo (PURSUIT)": ("GSE92415", "placebo")}
ARIJS_2009 = ["TNFRSF11B", "STC1", "PTGS2", "IL13RA2", "IL11"]


def _arm(cohort: str, treatment: str) -> tuple[pd.DataFrame, pd.Series]:
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    s = s[(s.cohort == cohort) & (s.timepoint == "baseline") & (s.group == "UC") & (s.treatment == treatment)
          & s.response.notna()].set_index("geo")
    x = pd.read_parquet(PROCESSED / f"expr_{cohort}.parquet")
    keep = [c for c in x.columns if c in s.index]
    return x[keep], s.loc[keep, "response"].astype(bool)


def response_genes() -> pd.DataFrame:
    rows = []
    for arm, (cohort, tr) in ARMS.items():
        x, r = _arm(cohort, tr)
        nr, rs = x.loc[:, ~r.values], x.loc[:, r.values]
        t, p = stats.ttest_ind(nr, rs, axis=1, equal_var=False)
        d = pd.DataFrame({"gene": x.index, "arm": arm, "cohort": cohort, "diff_nr_minus_r": nr.mean(1) - rs.mean(1),
                          "t": t, "p": p, "n_nonresponders": nr.shape[1], "n_responders": rs.shape[1]})
        d["padj"] = multipletests(d["p"].fillna(1), method="fdr_bh")[1]
        rows.append(d)
    out = pd.concat(rows, ignore_index=True)
    out.to_parquet(PROCESSED / "response_genes.parquet", index=False)
    return out


def interaction(genes: list[str]) -> pd.DataFrame:
    import statsmodels.api as sm
    rows = []
    for drug, (cohort, tr) in (("ustekinumab", ARMS["ustekinumab"]), ("golimumab", ARMS["golimumab"])):
        xd, rd = _arm(cohort, tr)
        xp, rp = _arm(cohort, "placebo")
        for g in genes:
            if g not in xd.index:
                continue
            e = np.concatenate([xd.loc[g].to_numpy(float), xp.loc[g].to_numpy(float)])
            e = (e - e.mean()) / (e.std() or 1)
            arm = np.concatenate([np.ones(xd.shape[1]), np.zeros(xp.shape[1])])
            y = np.concatenate([rd.to_numpy(float), rp.to_numpy(float)])
            X = sm.add_constant(np.column_stack([e, arm, e * arm]))
            try:
                f = sm.Logit(y, X).fit(disp=0)
                rows.append({"gene": g, "drug": drug, "cohort": cohort, "log_or_interaction": float(f.params[3]),
                             "p_interaction": float(f.pvalues[3]), "log_or_expression_placebo": float(f.params[1])})
            except Exception:
                continue
    d = pd.DataFrame(rows)
    d["padj_interaction"] = d.groupby("drug")["p_interaction"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    d.to_parquet(PROCESSED / "response_interaction.parquet", index=False)
    return d


def response_pathways() -> pd.DataFrame:
    import decoupler as dc
    nets = {"pathway": dc.op.progeny(organism="human", top=500), "tf": dc.op.collectri(organism="human")}
    rows = []
    for arm, (cohort, tr) in ARMS.items():
        x, r = _arm(cohort, tr)
        x = x.sub(x.mean(axis=1), axis=0)
        for kind, net in nets.items():
            est, _ = dc.mt.ulm(data=x.T, net=net, tmin=5, verbose=False)
            t, p = stats.ttest_ind(est.loc[~r.values], est.loc[r.values], equal_var=False)
            d = pd.DataFrame({"name": est.columns, "kind": kind, "arm": arm, "t_nr_vs_r": t, "p": p})
            d["padj"] = multipletests(d["p"].fillna(1), method="fdr_bh")[1]
            rows.append(d)
    out = pd.concat(rows, ignore_index=True)
    out.to_parquet(PROCESSED / "response_pathways.parquet", index=False)
    return out


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    today = date.today().isoformat()
    g = pd.read_parquet(PROCESSED / "response_genes.parquet")
    inter = pd.read_parquet(PROCESSED / "response_interaction.parquet")
    pw = pd.read_parquet(PROCESSED / "response_pathways.parquet")
    keep = set(pd.read_parquet(PROCESSED / "meta.parquet").query("meta_padj < 0.05")["gene"])
    n = 0
    with store.batch():
        store.add_node("ANALYSIS:treatment", "Dataset", "Baseline biopsy expression vs response to ustekinumab, "
                                                        "golimumab, infliximab, vedolizumab (and placebo)")
        for gene, d in g[g.gene.isin(keep)].groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            arms = {r.arm: {"diff_nr_minus_r": round(float(r.diff_nr_minus_r), 3), "padj": float(f"{r.padj:.3g}"),
                            "n_nr": int(r.n_nonresponders), "n_r": int(r.n_responders)} for r in d.itertuples()}
            it = inter[inter.gene == gene]
            rec = {"gene": gene, "arms": arms,
                   "interaction": {r.drug: {"log_or_interaction": round(r.log_or_interaction, 3),
                                            "padj": float(f"{r.padj_interaction:.3g}")} for r in it.itertuples()},
                   "reading": "diff_nr_minus_r > 0 = higher at baseline in non-responders; association, not prediction"}
            ev = store.add_evidence(evidence_id("treatment_response", gene), "treatment_response", "ANALYSIS:treatment",
                                    today, **rec)
            store.add_claim(node["id"], "response_assessed_in", "ANALYSIS:treatment", ev)
            n += 1
        for (name, kind), d in pw.groupby(["name", "kind"]):
            rec = {"name": name, "kind": kind, "arms": {r.arm: {"t_nr_vs_r": round(float(r.t_nr_vs_r), 2),
                                                               "padj": float(f"{r.padj:.3g}")} for r in d.itertuples()}}
            ev = store.add_evidence(evidence_id("response_pathway", kind, name), "response_pathway", "ANALYSIS:treatment", today, **rec)
            store.add_claim("ANALYSIS:treatment", "assessed", "ANALYSIS:treatment", ev) if False else None
            nid = f"{'PATHWAY' if kind == 'pathway' else 'TF'}ACT:{name}"
            if store.get_node(nid):
                store.add_claim(nid, "response_assessed_in", "ANALYSIS:treatment", ev)
                n += 1
    return n


if __name__ == "__main__":
    from .cli import expression_top
    g = response_genes()
    print(g.groupby("arm").apply(lambda d: int((d.padj < 0.05).sum())).to_dict())
    a = g[(g.arm == "infliximab") & g.gene.isin(ARIJS_2009)][["gene", "diff_nr_minus_r", "p", "padj"]]
    print("Arijs 2009 control:\n", a.round(4).to_string())
    genes = list(dict.fromkeys(expression_top(300) + ["TNF", "IL23A", "IL12B", "ITGB7", "ITGA4", "MADCAM1", "OSM", "OSMR", "TREM1", "JAK1", "S1PR1", "TNFSF15", "IL11"]))
    it = interaction(genes)
    print("interaction padj<0.05:", it.groupby("drug").apply(lambda d: int((d.padj_interaction < 0.05).sum())).to_dict())
    pw = response_pathways()
    print(pw[(pw.kind == "pathway") & (pw.padj < 0.05)].sort_values("t_nr_vs_r", ascending=False)[["name", "arm", "t_nr_vs_r", "padj"]].head(20).to_string())
    print("records:", load())
