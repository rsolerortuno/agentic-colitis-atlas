"""Does each statement say only what its cited evidence says? A reader subagent checks every claim (added 2026-11-23).

The deterministic verifier (verify.py, semantic.py) proves that a statement cites real evidence about the right gene and
that its numbers and direction words match signed values. It cannot tell whether a claim is contained in what is
cited: a review of the run of 2026-11-14 found about ten literature claims taken from an abstract or a paper that the
statement does not cite, two sources described wrongly, and direction flips that no rule can see (a comparison of AUCs,
a gene the statement cites nothing about). This module adds a reading check:

1. prepare(run): every statement of the run that passes verification, with the full records of the evidence it cites
   and, for literature evidence with a PubMed id, the title and abstract (PubTator3, data/literature/pubtator). The same
   is written for the direction-flipped copies of the verifier sabotage test (sabotage.corrupt_direction), in separate
   batches with opaque ids (key.json says which is which). JSONL batches in data/claim_check/<run>/batches (not in
   the repository: they contain abstracts).
2. The reader subagent (.claude/agents/lector-afirmaciones.md, Claude Code subscription) splits each statement into
   claims about data or papers and gives each one a verdict: supported, contradicted (the evidence says otherwise) or
   not_in_evidence (nothing cited says it), with a quote copied verbatim from the evidence for the first two.
   Interpretation, recommendations and placement in the ranking are not claims about the evidence and are skipped.
3. validate(): a supported or contradicted claim keeps its verdict only when its quote is verbatim in the evidence
   given to the reader (whitespace and case ignored); otherwise it is counted as not_in_evidence with the reason.
4. summarize(): per statement, the claims contradicted and not in evidence; for the sabotage copies, the share whose
   flipped claim is caught (a contradicted or not_in_evidence claim), compared with the deterministic verifier on the
   same copies. Outputs: data/claim_check/<run>/read/ (reader output, kept), evals/claim_check.json.

The reader is a language model, so its verdicts are measured (the sabotage copies) and shown as a review aid next to the
deterministic checks; it does not remove evidence or fail the run.

Revision loop (2026-11-23, steps 3b and 4c of .claude/commands/informe.md): the same reader runs on the specialists'
verified findings (prepare_agents -> data/claim_check/<run>/agents/) and on the narrative before assembly
(prepare_narrative -> .../narrative/); revisions() and revisions_narrative() write RUN/revisions_*.json with the claims
to fix, which go back to the specialist or the orchestrator once, before verification and the final check above.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from . import sabotage, semantic
from .paths import LITERATURE, ROOT

BASE = ROOT / "data" / "claim_check"
PUBTATOR = LITERATURE / "pubtator"
BATCH_SIZE = 10
MAX_RECORD = 3000
MAX_ABSTRACT = 4000
VERDICTS = {"supported", "contradicted", "not_in_evidence"}


def _abstract(pmid) -> str | None:
    p = PUBTATOR / f"{pmid}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    return " ".join(x.get("text", "") for x in d.get("passages", []))[:MAX_ABSTRACT] or None


def _evidence_block(e: dict) -> dict:
    rec = {k: v for k, v in e.items() if k not in ("date",)}
    out = {"evidence_id": e.get("evidence_id"), "record": json.dumps(rec, ensure_ascii=False, default=str)[:MAX_RECORD]}
    if e.get("pmid"):
        a = _abstract(e["pmid"])
        if a:
            out["abstract"] = a
    return out


def _verifier(run: str):
    from .provenance import load_basis
    from .verify import Verifier
    return Verifier(basis=load_basis(ROOT / "runs" / run))


def _passing(run: str):
    from .verify import _statements
    findings = json.loads((ROOT / "runs" / run / "findings.json").read_text(encoding="utf-8"))
    v = _verifier(run)
    for loc, s in _statements(findings):
        ok, _ = v.check_statement(s)
        if ok is None or len(ok["evidence_ids"]) != len(s.get("evidence_ids", [])):
            continue
        yield loc, s, v.store.get_evidence([i for i in s.get("evidence_ids", []) if isinstance(i, str)])


def prepare(run: str) -> dict:
    """Write the reader's input batches: real statements and their direction-flipped copies, in separate batches."""
    real, flipped = [], []
    for loc, s, ev in _passing(run):
        key = sabotage._text_key(s)
        item = {"gene": s.get("gene"), "evidence": [_evidence_block(e) for e in ev]}
        real.append({"id": f"real|{loc}", "statement": s[key], **item})
        m = sabotage.corrupt_direction(s)
        if m:
            flipped.append({"id": f"flip|{loc}", "statement": m[key], **item})
    n = _write_batches([real, flipped], BASE / run)
    return {"real": len(real), "flipped": len(flipped), "batches": n}


def _write_batches(groups: list[list[dict]], base: Path) -> int:
    """Batches with opaque ids; which statement is which (real or flipped, and where) stays in key.json, outside the
    batches the reader sees."""
    d = base / "batches"
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.jsonl"):
        f.unlink()
    key, n = {}, 0
    for rows in groups:
        for i in range(0, len(rows), BATCH_SIZE):
            n += 1
            chunk = []
            for r in rows[i:i + BATCH_SIZE]:
                oid = f"S{len(key) + 1:03d}"
                key[oid] = r["id"]
                chunk.append({**r, "id": oid})
            (d / f"batch_{n:02d}.jsonl").write_text(
                "\n".join(json.dumps(r, ensure_ascii=False) for r in chunk) + "\n", encoding="utf-8")
    (base / "key.json").write_text(json.dumps(key, indent=0), encoding="utf-8")
    return n


# ------------------------------------------------------------------ revision loop inside /informe (2026-11-23)
def _rows_for(statements: list[tuple[str, dict]], store) -> list[dict]:
    out = []
    for loc, s in statements:
        ids = [i for i in s.get("evidence_ids", []) if isinstance(i, str)]
        if not ids:
            continue
        out.append({"id": loc, "statement": s[sabotage._text_key(s)], "gene": s.get("gene"),
                    "evidence": [_evidence_block(e) for e in store.get_evidence(ids)]})
    return out


def prepare_agents(run: str) -> dict:
    """Batches of the specialists' verified findings, for the reading check before the narrative is written."""
    from .graph import GraphStore
    from .verify import AGENT_FILES
    stmts = []
    for name in AGENT_FILES:
        p = ROOT / "runs" / run / name.replace(".json", ".verified.json")
        if p.exists():
            for i, f in enumerate(json.loads(p.read_text(encoding="utf-8")).get("findings", [])):
                stmts.append((f"{name}|{i}", f))
    rows = _rows_for(stmts, GraphStore())
    n = _write_batches([rows], BASE / run / "agents")
    return {"findings": len(rows), "batches": n, "dir": str(BASE / run / "agents" / "batches")}


def _narrative_statements(nar: dict) -> list[tuple[str, dict]]:
    out = [(f"narrative.summary[{i}]", p) for i, p in enumerate(nar.get("summary", []))]
    for sec, ps in (nar.get("sections") or {}).items():
        out += [(f"narrative.sections.{sec}[{i}]", p) for i, p in enumerate(ps)]
    out += [(f"targets.{t.get('gene')}", t) for t in nar.get("targets", [])]
    return out


def prepare_narrative(run: str) -> dict:
    from .graph import GraphStore
    from .provenance import extend_basis
    nar = json.loads((ROOT / "runs" / run / "narrative.json").read_text(encoding="utf-8"))
    store = GraphStore()
    stmts = _narrative_statements(nar)
    # the orchestrator may cite records the specialists did not (a rank_context record for a comparison, a drug record):
    # they enter the run's evidence basis now, with their current fingerprint, as the specialists' citations did at
    # verify-agents time
    extend_basis(ROOT / "runs" / run, store, [i for _, s in stmts for i in s.get("evidence_ids", []) if isinstance(i, str)])
    rows = _rows_for(stmts, store)
    n = _write_batches([rows], BASE / run / "narrative")
    return {"paragraphs": len(rows), "batches": n, "dir": str(BASE / run / "narrative" / "batches")}


def _flagged(rows: list[dict], texts: dict[str, str]) -> dict[str, dict]:
    return {r["id"]: {"statement": texts.get(r["id"], ""),
                      "claims": [c for c in r["claims"] if c["verdict"] != "supported"]}
            for r in rows if any(c["verdict"] != "supported" for c in r["claims"])}


def _texts(base: Path) -> dict[str, str]:
    key = json.loads((base / "key.json").read_text(encoding="utf-8"))
    out = {}
    for f in sorted((base / "batches").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                out[key[r["id"]]] = r["statement"]
    return out


def revisions(run: str) -> dict:
    """RUN/revisions_<agent>.json: per finding, the claims the reader found contradicted or not in the cited records."""
    base = BASE / run / "agents"
    flagged = _flagged(validate(run, "agents"), _texts(base))
    by_agent: dict[str, list] = {}
    for loc, v in flagged.items():
        name, i = loc.split("|")
        by_agent.setdefault(name, []).append({"finding": int(i), **v})
    # findings the deterministic verifier dropped, with the reason, go back in the same file
    dropped: dict[str, list] = {}
    vp = ROOT / "runs" / run / "verification_agents.json"
    if vp.exists():
        keys = ("unmatched_numbers", "direction_contradiction", "numeric_binding_errors", "uncited_drugs",
                "stale_evidence", "nonexistent", "not_linked_to_gene")
        for name, r in json.loads(vp.read_text(encoding="utf-8"))["agents"].items():
            for d in r.get("dropped_detail", []):
                dropped.setdefault(name, []).append({"gene": d.get("gene"), "text_start": d.get("text"),
                                                     "why": {k: d[k] for k in keys if d.get(k)}})
    out = {}
    for name in sorted(set(by_agent) | set(dropped)):
        items = by_agent.get(name, [])
        p = ROOT / "runs" / run / f"revisions_{name.replace('agent_', '').replace('.json', '')}.json"
        p.write_text(json.dumps({"agent_file": name, "findings_verified_by_index": "index into the verified file",
                                 "claims_to_fix": items, "dropped_by_verifier": dropped.get(name, [])},
                                indent=1, ensure_ascii=False), encoding="utf-8")
        out[name] = {"findings_to_revise": len(items), "claims": sum(len(x["claims"]) for x in items),
                     "dropped": len(dropped.get(name, [])), "file": p.name}
    return out


def revisions_narrative(run: str) -> dict:
    base = BASE / run / "narrative"
    flagged = _flagged(validate(run, "narrative"), _texts(base))
    p = ROOT / "runs" / run / "revisions_narrative.json"
    p.write_text(json.dumps(flagged, indent=1, ensure_ascii=False), encoding="utf-8")
    return {"paragraphs_to_revise": len(flagged), "claims": sum(len(v["claims"]) for v in flagged.values())}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip().lower()


def validate(run: str, sub: str = "") -> list[dict]:
    """Reader lines checked against the batches: one line per statement, quotes verbatim in the evidence given."""
    base = BASE / run / sub if sub else BASE / run
    key = json.loads((base / "key.json").read_text(encoding="utf-8"))
    src = {}
    for f in sorted((base / "batches").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                src[r["id"]] = r
    rows, seen = [], set()
    for f in sorted((base / "read").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("id") not in src or r["id"] in seen:
                continue
            seen.add(r["id"])
            text = _norm(" ".join(e["record"] + " " + e.get("abstract", "") for e in src[r["id"]]["evidence"]))
            claims = []
            for c in r.get("claims") or []:
                v = c.get("verdict")
                if v not in VERDICTS:
                    continue
                c = {"claim": str(c.get("claim", "")), "verdict": v, "quote": str(c.get("quote") or ""),
                     "evidence_id": c.get("evidence_id"), "reason": str(c.get("reason") or "")}
                if v != "not_in_evidence" and (not c["quote"] or _norm(c["quote"]) not in text):
                    c.update(verdict="not_in_evidence", reason="quote not verbatim in the cited evidence; " + c["reason"])
                claims.append(c)
            rows.append({"id": key[r["id"]], "claims": claims})
    missing = sorted(set(src) - seen)
    if missing:
        raise ValueError(f"{len(missing)} statements have no reader line, e.g. {missing[:3]}")
    return rows


def summarize(run: str, write: bool = True) -> dict:
    rows = {r["id"]: r for r in validate(run)}
    key = json.loads((BASE / run / "key.json").read_text(encoding="utf-8"))
    real = {k.split("|", 1)[1]: r for k, r in rows.items() if k.startswith("real|")}
    flip = {k.split("|", 1)[1]: r for k, r in rows.items() if k.startswith("flip|")}

    def bad(r, v):
        return [c for c in r["claims"] if c["verdict"] == v]

    src = {}
    for f in sorted((BASE / run / "batches").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                src[key[r["id"]]] = r["statement"]

    # a flipped copy is caught when a claim the reader flags (contradicted or not in evidence) contains the flipped
    # word: the reader must have read that very claim as wrong, not something else in the statement
    def caught(loc):
        a, b = src[f"real|{loc}"], src[f"flip|{loc}"]
        i = next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)
        while i > 0 and (b[i - 1].isalnum() or b[i - 1] == "-"):
            i -= 1
        word = re.match(r"[\w-]+", b[i:]).group(0)
        return any(re.search(rf"\b{re.escape(word)}\b", c["claim"], re.I)
                   for c in bad(flip[loc], "contradicted") + bad(flip[loc], "not_in_evidence"))

    det = json.loads((ROOT / "evals" / "verifier_sabotage.json").read_text(encoding="utf-8"))["corruptions"]["direction"]
    caught_locs = sorted(loc for loc in flip if loc in real and caught(loc))
    out = {
        "run": run,
        "reader": ".claude/agents/lector-afirmaciones.md (Claude Code subscription)",
        "statements": len(real),
        "claims": sum(len(r["claims"]) for r in real.values()),
        "claims_contradicted": sum(len(bad(r, "contradicted")) for r in real.values()),
        "claims_not_in_evidence": sum(len(bad(r, "not_in_evidence")) for r in real.values()),
        "statements_with_a_contradicted_claim": sorted(k for k, r in real.items() if bad(r, "contradicted")),
        "statements_with_an_unsupported_claim": sorted(k for k, r in real.items() if bad(r, "not_in_evidence")),
        "flagged": {k: [c for c in r["claims"] if c["verdict"] != "supported"] for k, r in sorted(real.items())
                    if any(c["verdict"] != "supported" for c in r["claims"])},
        "sabotage": {"flipped": len(flip), "caught_by_reader": len(caught_locs),
                     "caught_by_deterministic_checks": det["caught"], "deterministic_applied": det["applied"],
                     "caught_locations": caught_locs},
        "scope": "A language-model reader; verdicts are a review aid next to the deterministic checks, measured on "
                 "the direction-flipped copies of the sabotage test.",
    }
    if write:
        (ROOT / "evals" / "claim_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                         encoding="utf-8")
    return out


if __name__ == "__main__":
    import sys
    cmd, run = sys.argv[1], sys.argv[2]
    fn = {"prepare": prepare, "summarize": summarize, "prepare-agents": prepare_agents, "revisions": revisions,
          "prepare-narrative": prepare_narrative, "revisions-narrative": revisions_narrative}[cmd]
    print(json.dumps(fn(run), indent=1, ensure_ascii=False)[:4000])
