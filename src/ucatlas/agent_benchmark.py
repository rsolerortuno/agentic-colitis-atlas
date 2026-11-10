"""M17: one agent vs the team on the full pipeline (same question, data, tools, verifier and narrative rules).

Team (/informe): Opus orchestrator + six Sonnet specialists in parallel; three runs (main + two repeats).
One agent (/informe_single): one Opus agent doing the six specialists' work itself; one run.
Measures: verification (statements, unsupported, stale), coverage of the agents' evidence in the narrative, a
checklist of 12 facts the data establish (fixed before reading the single-agent report; regexes on the narrative),
overlap of the target lists, cost, wall time and turns from Claude Code's JSON result.
Output: evals/agent_benchmark.json.
"""
from __future__ import annotations

import itertools
import json
import re

from .paths import ROOT

RUNS = {"team": ["2026-11-09_full", "2026-11-09_rep1", "2026-11-09_rep2"], "single": ["2026-11-09_single"]}


def _para(text: str, *pats: str) -> bool:
    return any(all(re.search(p, para, re.I) for p in pats) for para in text.split("\n\n"))


CHECKS = {
    "neutrophil_influx_explains_cxcr2": lambda t, tg: _para(t, r"\bCXCR2\b|\bCXCL8\b", r"neutrophil"),
    "composition_adjustment_used": lambda t, tg: bool(re.search(r"composition|adjust", t, re.I)),
    "epithelial_loss_biomarkers": lambda t, tg: bool({"HMGCS2", "AQP8", "SLC26A2", "CLDN8"} & set(tg)),
    "osm_osmr_circuit": lambda t, tg: _para(t, r"\bOSM\b", r"\bOSMR\b"),
    "tnf_nfkb_in_nonresponders": lambda t, tg: _para(t, r"TNF|NFkB|NF-κB|NF-kB", r"non-?respond"),
    "placebo_severity_caveat": lambda t, tg: bool(re.search(r"placebo", t, re.I)),
    "two_exploratory_targets": lambda t, tg: True,   # filled from the targets' rationales below
    "known_targets_as_controls": lambda t, tg: bool(re.search(r"IL-?23|\bJAK|integrin|ITGB7|TNF\b", t)),
    "il11_inflammatory_fibroblasts": lambda t, tg: _para(t, r"\bIL11\b|inflammatory fibroblast", r"fibroblast"),
    "infliximab_signature": lambda t, tg: bool(re.search(r"IL13RA2|Arijs|TNFRSF11B", t)),
    "replication_across_cohorts": lambda t, tg: bool(re.search(r"cohorts", t, re.I)),
    "genetic_support_named": lambda t, tg: _para(t, r"genetic", r"\bJAK2\b|\bFCGR2A\b|\bIL23R\b|\bOSMR\b"),
}


def _session(run: str) -> dict:
    b = (ROOT / "runs" / run / "session.json").read_bytes()
    return json.loads(b.decode("utf-16") if b[:2] == b"\xff\xfe" else b.decode("utf-8"))


def _text(n: dict) -> str:
    parts = [p.get("text", "") for p in n.get("summary", [])]
    for v in (n.get("sections") or {}).values():
        parts += [p.get("text", "") for p in v]
    parts += [t.get("rationale", "") for t in n.get("targets", [])]
    return "\n\n".join(parts)


def score(run: str) -> dict:
    d = ROOT / "runs" / run
    n = json.loads((d / "narrative.json").read_text(encoding="utf-8"))
    v = json.loads((d / "verification_final.json").read_text(encoding="utf-8"))
    cov = json.loads((d / "coverage.json").read_text(encoding="utf-8"))
    s = _session(run)
    targets = [t["gene"] for t in n.get("targets", [])]
    text = _text(n)
    checks = {k: bool(f(text, targets)) for k, f in CHECKS.items()}
    checks["two_exploratory_targets"] = sum(t.get("rationale", "").startswith("Exploratory") for t in n.get("targets", [])) >= 2
    return {"run": run, "statements": v["statements"], "unsupported": v["unsupported"], "stale": v.get("stale"),
            "key_findings": cov["key_findings"], "key_findings_uncited": len(cov["key_findings_not_cited"]),
            "evidence_cited_share": cov["share_of_finding_evidence_cited"], "targets": targets,
            "checklist": checks, "checklist_score": sum(checks.values()), "checklist_max": len(checks),
            "cost_usd": round(s.get("total_cost_usd", 0), 2), "minutes": round(s["duration_ms"] / 60000, 1),
            "turns": s["num_turns"], "models": sorted((s.get("modelUsage") or {}).keys()),
            "words": len(text.split())}


def run() -> dict:
    out = {arm: [score(r) for r in runs] for arm, runs in RUNS.items()}
    allr = [r for v in out.values() for r in v]
    jac = {}
    for a, b in itertools.combinations(allr, 2):
        A, B = set(a["targets"]), set(b["targets"])
        jac[f"{a['run']} vs {b['run']}"] = round(len(A & B) / len(A | B), 2)
    out["target_jaccard"] = jac
    (ROOT / "evals" / "agent_benchmark.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    r = run()
    for arm in RUNS:
        for x in r[arm]:
            print(arm, x["run"], x["statements"], x["unsupported"], x["checklist_score"], "/", x["checklist_max"],
                  [k for k, v in x["checklist"].items() if not v], x["evidence_cited_share"], x["cost_usd"], x["minutes"], x["turns"])
    print(r["target_jaccard"])
