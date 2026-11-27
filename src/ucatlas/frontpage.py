"""Which genes may be on the front page, by a fixed rule (added 2026-11-27).

Until now the orchestrator chose the 10-14 prioritised candidates by judgment, and a reader could not tell why JAK2
(rank 36 of the convergence ranking) was there while genes ranked above it were not. The rule below makes the choice
checkable. A gene may be on the front page when it meets at least one criterion, and the card says which:

- control: a gene of the convergence ranking (top 100) on which, or on whose class, a drug approved in UC acts
  (drug_indications.parquet, UC at approval), with an Open Targets genetic association with UC of at least
  CONTROL_GENETICS. These are mandatory: each must be on the
  page, as a positive control, whatever its convergence rank, so the method is checked against what is known.
- convergence: among the first TOP_RANK genes of the convergence ranking.
- exploratory: the first first-in-class candidate of the discovery screen by pooled log2FC, or by the share of the
  change retained after composition adjustment.
- genetics: a significant Mendelian randomisation estimate (q < 0.05) in plasma (pQTL) or gut (eQTL) that
  colocalises or replicates, the genetic argument for or against blocking the gene.

The orchestrator still orders the genes and writes the rationale; verify.py fails a run whose front page has a gene
that meets no criterion or misses a mandatory control.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .paths import PROCESSED

CONTROL_GENETICS = 0.8
TOP_RANK = 30

LABEL = {
    "control": "positive control: drug approved in UC on this gene or its class, genetic association {g:.2f}",
    "convergence": "convergence rank {r} (first {n})",
    "exploratory_lfc": "first-in-class candidate with the largest pooled increase",
    "exploratory_retained": "first-in-class candidate with the highest retained share after composition adjustment",
    "genetics": "Mendelian randomisation q < 0.05 ({layer})",
}

RULE_TEXT = (f"A gene is on this page when it meets at least one fixed criterion: a drug approved in UC on the gene or "
             f"its class with Open Targets genetic association of at least {CONTROL_GENETICS} among the ranked genes (a positive control, "
             f"always shown); one of the first {TOP_RANK} of the convergence ranking; the leading first-in-class "
             f"candidate by pooled increase or by retained share; or a significant Mendelian randomisation estimate. "
             f"The orchestrator chooses and orders among them.")


def _mr_genes() -> dict[str, str]:
    out = {}
    p = PROCESSED / "mr_pqtl.parquet"
    if p.exists():
        m = pd.read_parquet(p)
        ok = (m["q"] < 0.05) & (m["colocalises"].fillna(False) | m.get("replicates_finngen", False).fillna(False))
        out.update({g: "plasma protein" for g in m.loc[ok, "gene"]})
    p = PROCESSED / "mr_colon.parquet"
    if p.exists():
        m = pd.read_parquet(p)
        ok = (m["q"] < 0.05) & (m["colocalises"].fillna(False) | m["replicates_finngen"].fillna(False))
        for g in m.loc[ok, "gene"]:
            out.setdefault(g, "gut expression")
    return out


def criteria(findings: dict) -> dict[str, list[str]]:
    """gene -> the criteria it meets (labels), for every gene that meets one."""
    ranking = findings["convergence"]["ranking"]
    out: dict[str, list[str]] = {}
    ot = {r["gene"]: r.get("ot_genetic_association") or 0.0 for r in ranking}
    di = pd.read_parquet(PROCESSED / "drug_indications.parquet")
    approved = di[(di.category == "UC") & (di.uc_stage == "APPROVAL")].gene
    for g in approved:
        if ot.get(g, 0.0) >= CONTROL_GENETICS:
            out.setdefault(g, []).append(LABEL["control"].format(g=ot[g]))
    for r in ranking[:TOP_RANK]:
        out.setdefault(r["gene"], []).append(LABEL["convergence"].format(r=r["rank"], n=TOP_RANK))
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    f = d[d.status == "first-in-class opportunity"]
    if len(f):
        out.setdefault(f.sort_values("meta_log2fc", ascending=False).gene.iloc[0], []).append(LABEL["exploratory_lfc"])
        out.setdefault(f.sort_values("retained", ascending=False).gene.iloc[0], []).append(LABEL["exploratory_retained"])
    for g, layer in _mr_genes().items():
        out.setdefault(g, []).append(LABEL["genetics"].format(layer=layer))
    return out


def mandatory(findings: dict) -> list[str]:
    return sorted(g for g, c in criteria(findings).items() if any(x.startswith("positive control") for x in c))


def check(findings: dict) -> dict:
    """The front-page genes against the rule: genes meeting no criterion, and mandatory controls left out."""
    crit = criteria(findings)
    genes = [t["gene"] for t in findings.get("targets", [])]
    return {"rule": RULE_TEXT, "genes": len(genes),
            "unjustified": [g for g in genes if g not in crit],
            "missing_controls": [g for g in mandatory(findings) if g not in genes],
            "reasons": {g: crit.get(g, []) for g in genes}}


def write(run_dir: Path) -> Path:
    """Eligible genes for the orchestrator (runs/<run>/frontpage_rule.json), from the agent-step convergence."""
    conv = json.loads((run_dir / "convergence.json").read_text(encoding="utf-8"))
    f = {"convergence": conv}
    out = run_dir / "frontpage_rule.json"
    out.write_text(json.dumps({"rule": RULE_TEXT, "mandatory_controls": mandatory(f), "eligible": criteria(f)},
                              indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    import sys
    from .runs import current_run, run_dir
    d = run_dir(sys.argv[1]) if len(sys.argv) > 1 else current_run()
    p = write(d)
    r = json.loads(p.read_text(encoding="utf-8"))
    print(f"{p}: {len(r['eligible'])} eligible genes; mandatory controls: {', '.join(r['mandatory_controls']) or 'none'}")
