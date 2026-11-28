"""The orchestrator's decision: four targets to test first from each list of targets with no therapy in UC (2026-11-28).

The atlas answers "which new targets should be tested first?". Four lists hold targets that no drug treats in UC: the
first-in-class shortlist (no drug for any disease), well known in UC with no drug, repurposing (an inhibitor exists
for another disease, never tried in UC) and genetics first (inherited risk variants). Each list sorts its genes into
groups A, B and C by how stable they are under several ways of scoring. The orchestrator chooses four genes from
group A of each list and explains each choice with fixed criteria (direction, cell, drug route, safety, first
experiment); a reviewer subagent then gives a cited verdict per candidate (expert_review.py).

The rule, checked by verify.py final:
- exactly PER_LIST genes per list, each from the list's pool (group A; group B only when group A has fewer than
  PER_LIST genes, marked);
- no gene that already has a drug tried or approved in UC or IBD (clinical_list.py);
- no gene twice; each choice cites evidence.

The pool and a dossier per gene (lights, values, cell, drugs, papers, proposed experiment, evidence ids) are written
to RUN/decision_pool.json for the specialists and the orchestrator. This replaces the front-page selection rule
(frontpage.py, 2026-11-27): positive controls (IL12B, JAK2) left the decision and stay as reference in the clinical
list, which marks every approved target.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .graph import evidence_id
from .paths import PROCESSED

PER_LIST = 4
PREVIOUS_RUN = "2026-11-27"   # the last run with the old selection, for "genes considered and not chosen"

# key, title, file, literature file, evidence types cited by the list's own table
LISTS = [
    ("first_in_class", "First-in-class shortlist", "validation_shortlist", "shortlist_literature",
     ["discovery_profile", "shortlist_literature"]),
    ("known", "Well known in UC, no drug yet", "validation_known", "known_literature",
     ["discovery_profile", "known_literature"]),
    ("repurposing", "Repurposing: an inhibitor for another disease", "validation_repurposing", "repurposing_literature",
     ["discovery_profile", "clinical_indications", "repurposing_literature"]),
    ("genetic", "Genetics first", "genetic_targets", None, ["genetic_target"]),
]
TITLE = {k: t for k, t, *_ in LISTS}
CRITERIA = ["direction", "cell", "drug route", "safety", "first experiment"]

RULE_TEXT = (f"The orchestrator chose {PER_LIST} genes from group A of each list of targets with no therapy in UC "
             f"(first-in-class, well known with no drug, repurposing, genetics first): group A holds the genes that "
             f"stay among the list's first ten under every way of scoring. No chosen gene has a drug tried or approved "
             f"in UC or IBD. Each choice is explained by its direction, cell, drug route, safety and first experiment, "
             f"and a reviewer gives a verdict on each, with PubMed citations for every fact it adds.")


def _clinical_genes() -> set[str]:
    p = PROCESSED / "clinical_list.parquet"
    return set(pd.read_parquet(p).gene) if p.exists() else set()


def _papers(lit_file: str | None) -> dict[str, list[dict]]:
    if not lit_file or not (PROCESSED / f"{lit_file}.parquet").exists():
        return {}
    lit = pd.read_parquet(PROCESSED / f"{lit_file}.parquet").set_index("gene")
    return {g: json.loads(r.key_papers) for g, r in lit.iterrows()}


def _biomarker_papers() -> dict[str, list[str]]:
    """Genes with a biomarker result in an IBD setting in any of the literature reviews: gene -> PMIDs."""
    out: dict[str, list[str]] = {}
    for f in ("shortlist_literature", "known_literature", "repurposing_literature"):
        p = PROCESSED / f"{f}.parquet"
        if not p.exists():
            continue
        for g, r in pd.read_parquet(p).set_index("gene").iterrows():
            for x in json.loads(r.findings):
                if x.get("kind") == "biomarker" and x.get("ibd_setting") and x.get("pmid"):
                    out.setdefault(g, []).append(str(x["pmid"]))
    return out


def pool() -> dict[str, list[dict]]:
    """list key -> its pool, in the list's own order, each gene with its dossier."""
    clin = _clinical_genes()
    bio = _biomarker_papers()
    out: dict[str, list[dict]] = {}
    for key, title, f, lit_f, ev_types in LISTS:
        p = PROCESSED / f"{f}.parquet"
        if not p.exists():
            out[key] = []
            continue
        d = pd.read_parquet(p).sort_values("rank")
        d = d[~d.gene.isin(clin)]
        a = d[d.group == "A"]
        take = a if len(a) >= PER_LIST else pd.concat([a, d[d.group == "B"].head(PER_LIST - len(a))])
        papers = _papers(lit_f)
        rows = []
        for r in take.itertuples():
            g = r.gene
            dossier = {"gene": g, "list": key, "group": r.group, "from_group_b": r.group != "A",
                       "evidence_ids": [evidence_id(t, g) for t in ev_types],
                       "biomarker_pmids": sorted(set(bio.get(g, [])))}
            if key == "genetic":
                dossier.update({
                    "direction": r.hypothesis, "risk_allele": r.risk_allele, "tier": r.tier,
                    "predicted_effect": r.contexts, "mendelian_randomisation": r.mr,
                    "biopsy_change": ("no change" if r.tissue_change == "no change" or pd.isna(r.meta_log2fc)
                                      else f"{r.meta_log2fc:+.2f} log2FC, {r.tissue_change}"),
                    "cell_in_single_cells": r.atlas_lineage if isinstance(r.atlas_lineage, str) else None,
                    "cell_matches_prediction": bool(r.cell_match) if pd.notna(r.cell_match) else None,
                    "drug": r.drug, "points": int(r.points)})
            else:
                crit = json.loads(r.criteria)
                dossier.update({
                    "direction": "inhibit (rises in UC)" if (r.meta_log2fc or 0) > 0 else "restore (falls in UC)",
                    "stability": f"in the top 10 in {int(r.stability)} of 6 scorings",
                    "cell_family": r.cell_family, "meta_log2fc": round(float(r.meta_log2fc), 2),
                    "uc_papers": int(r.uc_papers),
                    "lights": {k: {"light": v[0], "detail": v[1]} for k, v in crit.items()},
                    "proposed_experiment": r.experiment,
                    "key_papers": [{k: x.get(k) for k in ("pmid", "year", "kind", "direction", "summary")}
                                   for x in papers.get(g, [])]})
                if key == "repurposing":
                    dossier["inhibitor"] = {"stage": r.inhibitor_stage, "drugs": list(r.inhibitor_drugs),
                                            "made_for_cancer": bool(r.cancer_drug)}
            rows.append(dossier)
        out[key] = rows
    return out


def previous_selection() -> list[dict]:
    """Genes of the last run's selection that are not in the decision pool, and where each one is now."""
    from .paths import ROOT
    p = ROOT / "runs" / PREVIOUS_RUN / "findings.json"
    if not p.exists():
        return []
    prev = [t["gene"] for t in json.loads(p.read_text(encoding="utf-8")).get("targets", [])]
    pl = {g["gene"]: k for k, rows in pool().items() for g in rows}
    clin = pd.read_parquet(PROCESSED / "clinical_list.parquet").set_index("gene") \
        if (PROCESSED / "clinical_list.parquet").exists() else pd.DataFrame()
    where = {}
    for key, title, f, *_ in LISTS:
        fp = PROCESSED / f"{f}.parquet"
        if fp.exists():
            for r in pd.read_parquet(fp).itertuples():
                where.setdefault(r.gene, (key, r.group))
    out = []
    for g in prev:
        if g in pl:
            continue
        if g in clin.index:
            st = "approved" if bool(clin.at[g, "approved"]) else "tried"
            reason = (f"a drug acting on it is {st} in UC or IBD, so it is not a target without therapy; it stays in "
                      f"the clinical list as reference")
            place = "In clinical validation"
        elif g in where:
            k, grp = where[g]
            reason = f"group {grp} of its list: its lights do not keep it among the first ten under every scoring"
            place = f"{TITLE[k]}, group {grp}"
        else:
            reason = "in none of the target lists (the screens set it aside)"
            place = "none"
        out.append({"gene": g, "now_in": place, "reason": reason})
    return out


def check(findings: dict) -> dict:
    """The decision against the rule."""
    P = pool()
    allowed = {k: {g["gene"] for g in rows} for k, rows in P.items()}
    clin = _clinical_genes()
    targets = findings.get("targets", [])
    problems = []
    seen = set()
    per = {k: 0 for k in allowed}
    for t in targets:
        g, k = t.get("gene"), t.get("list")
        if k not in allowed:
            problems.append(f"{g}: list {k!r} is not one of {sorted(allowed)}")
            continue
        per[k] += 1
        if g not in allowed[k]:
            problems.append(f"{g}: not in the pool of {k}")
        if g in clin:
            problems.append(f"{g}: has a drug tried or approved in UC or IBD")
        if g in seen:
            problems.append(f"{g}: chosen twice")
        seen.add(g)
        if not t.get("evidence_ids"):
            problems.append(f"{g}: cites no evidence")
    for k, n in per.items():
        if allowed[k] and n != min(PER_LIST, len(allowed[k])):
            problems.append(f"{k}: {n} genes chosen, {min(PER_LIST, len(allowed[k]))} required")
    return {"rule": RULE_TEXT, "chosen": len(targets), "per_list": per, "problems": problems}


def write(run_dir: Path) -> Path:
    P = pool()
    out = run_dir / "decision_pool.json"
    out.write_text(json.dumps({
        "rule": RULE_TEXT, "per_list": PER_LIST, "criteria": CRITERIA,
        "lists": [{"key": k, "title": t, "pool": P.get(k, [])} for k, t, *_ in LISTS],
        "previous_selection_not_in_pool": previous_selection()}, indent=1, ensure_ascii=False), encoding="utf-8")
    return out


def view(findings: dict, verdicts: dict | None = None) -> dict:
    """The decision for the report: per list, the chosen genes with their dossier and the reviewer's verdict, and the
    pool genes not chosen."""
    P = pool()
    verdicts = verdicts or {}
    chosen = {t["gene"]: t for t in findings.get("targets", []) if t.get("list")}
    blocks = []
    for k, title, *_ in LISTS:
        rows = P.get(k, [])
        picks = [t for t in findings.get("targets", []) if t.get("list") == k]
        dos = {r["gene"]: r for r in rows}
        blocks.append({"key": k, "title": title,
                       "cards": [{**t, "dossier": dos.get(t["gene"], {}), "verdict": verdicts.get(t["gene"])}
                                 for t in picks],
                       "not_chosen": [r for r in rows if r["gene"] not in chosen]})
    return {"rule": RULE_TEXT, "blocks": blocks, "previous": previous_selection(),
            "n": len(chosen), "pool": sum(len(v) for v in P.values())}


if __name__ == "__main__":
    import sys
    from .runs import current_run, run_dir
    d = run_dir(sys.argv[1]) if len(sys.argv) > 1 else current_run()
    p = write(d)
    r = json.loads(p.read_text(encoding="utf-8"))
    print(f"{p}: " + ", ".join(f"{L['key']} {len(L['pool'])}" for L in r["lists"]) +
          f"; previous selection not in the pool: {len(r['previous_selection_not_in_pool'])}")
