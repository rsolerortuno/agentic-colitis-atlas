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
    d = BASE / run / "batches"
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.jsonl"):
        f.unlink()
    real, flipped = [], []
    for loc, s, ev in _passing(run):
        key = sabotage._text_key(s)
        item = {"gene": s.get("gene"), "evidence": [_evidence_block(e) for e in ev]}
        real.append({"id": f"real|{loc}", "statement": s[key], **item})
        m = sabotage.corrupt_direction(s)
        if m:
            flipped.append({"id": f"flip|{loc}", "statement": m[key], **item})
    # opaque ids: which copy is real and which is flipped stays in key.json, outside the batches the reader sees
    key, n = {}, 0
    for rows in (real, flipped):
        for i in range(0, len(rows), BATCH_SIZE):
            n += 1
            chunk = []
            for r in rows[i:i + BATCH_SIZE]:
                oid = f"S{len(key) + 1:03d}"
                key[oid] = r["id"]
                chunk.append({**r, "id": oid})
            (d / f"batch_{n:02d}.jsonl").write_text(
                "\n".join(json.dumps(r, ensure_ascii=False) for r in chunk) + "\n", encoding="utf-8")
    (BASE / run / "key.json").write_text(json.dumps(key, indent=0), encoding="utf-8")
    return {"real": len(real), "flipped": len(flipped), "batches": n}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip().lower()


def validate(run: str) -> list[dict]:
    """Reader lines checked against the batches: one line per statement, quotes verbatim in the evidence given."""
    key = json.loads((BASE / run / "key.json").read_text(encoding="utf-8"))
    src = {}
    for f in sorted((BASE / run / "batches").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                src[r["id"]] = r
    rows, seen = [], set()
    for f in sorted((BASE / run / "read").glob("*.jsonl")):
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
    print(json.dumps(prepare(run) if cmd == "prepare" else summarize(run), indent=1, ensure_ascii=False)[:4000])
