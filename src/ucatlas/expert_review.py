"""Expert review of the front page, with every factual objection cited (added 2026-11-27).

The verifier and the reading check (verify.py, claim_check.py) prove that each claim is in the evidence it cites. They
cannot see a claim that is true of its record but wrong or misleading to a UC clinician: an anti-p19 antibody called
an IL12B drug, a JAK2 card read as a JAK2-drug proposal, a failed drug class presented as open. Four independent
expert-style reviews of run 2026-11-23 found such problems by hand. This module makes that review a step of /informe:

1. prepare(run): the front page as a reader sees it, from RUN/narrative.json: the summary paragraphs and the
   prioritised candidates with their block, title and drug tag, numbers written out, in data/expert_review/<run>/
   input.json, with opaque ids (key.json maps them back).
2. The reviewer subagent (.claude/agents/revisor-experto.md) reads it as a UC clinician would, with PubMed search, and
   writes review.json. Two kinds of objection:
   - factual: the text says something the literature contradicts, or leaves out a fact a clinician would expect.
     It must give a PMID and a sentence copied verbatim from that paper's abstract.
   - framing: emphasis, title, order, a missing caveat. It may only ask to reword, reorder or remove what the text
     already says; it adds no fact, so it needs no citation.
3. validate(run): a factual objection is kept only when its quote is verbatim in the PubMed abstract and the paper is a
   peer-reviewed journal article, not retracted (pubmed_quotes.problems). Kept quotes are registered as pubmed_quote
   evidence, so the orchestrator can cite them when it corrects the text. Writes RUN/revisions_expert.json (kept and
   discarded objections, with reasons). The orchestrator applies them once; the corrected text then goes through the
   reading check and verification like any other, which is what stops a framing change from adding an uncited fact.

The reviewer is a language model, so it is measured (seed() and score()): framing and factual errors of the kind the
manual reviews found are written into copies of real cards, mixed blind with unaltered cards, and the share caught is
reported with the objections raised on unaltered cards (evals/expert_review.json).
"""
from __future__ import annotations

import json
import random
import re

from .paths import ROOT

BASE = ROOT / "data" / "expert_review"
KINDS = {"factual", "framing"}


# ------------------------------------------------------------------ the front page as text
def plain(stmt: dict, store) -> str:
    """The statement as the report prints it: numeric annotations replaced by their values, other brackets removed."""
    from .semantic import numeric_claim_text
    text = stmt.get("statement") or stmt.get("text") or stmt.get("rationale") or ""
    recs = {e["evidence_id"]: e for e in store.get_evidence(
        [c["evidence_id"] for c in stmt.get("numeric_claims", []) if c.get("evidence_id")])} if stmt.get("numeric_claims") else {}
    for c in stmt.get("numeric_claims", []):
        r = recs.get(c.get("evidence_id"))
        if r is None:
            continue
        try:
            ann = numeric_claim_text(c, r)
        except (ValueError, TypeError):
            continue
        value = c["display"] + ("%" if c.get("unit") in ("percent", "native_percent") else "")
        text = text.replace(ann, f"({'.'.join(map(str, c['field']))} = {value})", 1)
    return re.sub(r"\s+", " ", re.sub(r"\s*\[[^\]]*\]", "", text)).strip()


def front_page(nar: dict) -> list[dict]:
    """Summary paragraphs and candidate cards, with block and tags as rendered."""
    import pandas as pd
    from .graph import GraphStore
    from .paths import PROCESSED
    from .render import TARGET_GROUPS, drug_indication_view, target_groups, drug_tags
    from .drug_literature import DATA as DRUG_LIT
    store = GraphStore()
    out = [{"location": f"narrative.summary[{i}]", "kind": "summary paragraph", "text": plain(p, store)}
           for i, p in enumerate(nar.get("summary", []))]
    conv = json.loads((_run_dir(nar) / "convergence.json").read_text(encoding="utf-8"))
    f = {"convergence": conv, "targets": nar.get("targets", [])}
    disc = pd.read_parquet(PROCESSED / "discovery_scored.parquet").set_index("gene")
    tags = drug_tags()
    if DRUG_LIT.exists():
        tags.update(json.loads(DRUG_LIT.read_text(encoding="utf-8")).get("tags", {}))
    rank = {r["gene"]: r["rank"] for r in conv["ranking"]}
    titles = {k: t for k, t, _ in TARGET_GROUPS}
    for g in target_groups(f, drug_indication_view(), disc):
        for t in g["targets"]:
            out.append({"location": f"targets.{t['gene']}", "kind": "candidate card", "block": titles[g["key"]],
                        "gene": t["gene"], "title": t["gene"] + (f" · {t['subtitle']}" if t.get("subtitle") else ""),
                        "role": t.get("category"), "convergence_rank": rank.get(t["gene"]),
                        "drug_tag": tags.get(t["gene"]), "selected_because": t.get("selected_because"),
                        "text": plain(t, store)})
    return out


_RUN: dict = {}


def _run_dir(nar: dict):
    return _RUN["dir"]


def _write(items: list[dict], base, header: dict) -> dict:
    base.mkdir(parents=True, exist_ok=True)
    key = {}
    rows = []
    for i, it in enumerate(items, 1):
        oid = f"C{i:02d}"
        key[oid] = it["location"]
        rows.append({"id": oid, **{k: v for k, v in it.items() if k != "location"}})
    (base / "input.json").write_text(json.dumps({**header, "items": rows}, indent=1, ensure_ascii=False), encoding="utf-8")
    (base / "key.json").write_text(json.dumps(key, indent=0), encoding="utf-8")
    return key


HEADER = {"question": "What are the candidate therapeutic targets, immune pathways and biomarkers in ulcerative colitis "
                      "with the strongest convergent evidence?",
          "what_you_read": "The front page of an ulcerative colitis evidence atlas: the summary of the findings and the "
                           "prioritised candidates, in the blocks the page shows them in."}


def prepare(run: str) -> dict:
    from .frontpage import RULE_TEXT
    from .render import TARGET_GROUPS
    rd = ROOT / "runs" / run
    _RUN["dir"] = rd
    nar = json.loads((rd / "narrative.json").read_text(encoding="utf-8"))
    items = front_page(nar)
    header = {**HEADER, "selection_rule": RULE_TEXT,
              "blocks": [{"title": t, "note": n} for _, t, n in TARGET_GROUPS]}
    _write(items, BASE / run, header)
    return {"items": len(items), "input": str(BASE / run / "input.json"),
            "output": str(BASE / run / "review.json")}


# ------------------------------------------------------------------ validation of the reviewer's objections
def check_objections(objs: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split objections into kept and discarded (with the reason), without writing anything."""
    from .pubmed_quotes import fetch_many, problems
    kept, dropped = [], []
    pmids = sorted({int(o["pmid"]) for o in objs if o.get("kind") == "factual" and str(o.get("pmid", "")).isdigit()})
    recs = fetch_many(pmids) if pmids else {}
    for o in objs:
        kind = o.get("kind")
        if kind not in KINDS:
            dropped.append({**o, "reason": f"unknown kind {kind!r}"})
        elif kind == "framing":
            kept.append(o)
        elif not str(o.get("pmid", "")).isdigit() or not o.get("quote"):
            dropped.append({**o, "reason": "factual objection without a PMID and a verbatim quote"})
        elif int(o["pmid"]) not in recs:
            dropped.append({**o, "reason": f"PMID {o['pmid']} not found in PubMed"})
        else:
            p = problems(recs[int(o["pmid"])], [o["quote"]])
            (dropped if p else kept).append({**o, **({"reason": "; ".join(p)} if p else {})})
    return kept, dropped


def validate(run: str) -> dict:
    """Keep the objections that may be acted on, register their quotes as citable evidence, and write
    RUN/revisions_expert.json for the orchestrator."""
    from .pubmed_quotes import add
    base = BASE / run
    key = json.loads((base / "key.json").read_text(encoding="utf-8"))
    review = json.loads((base / "review.json").read_text(encoding="utf-8"))
    objs = [{**o, "location": key.get(o.get("id"), o.get("id"))} for o in review.get("objections", [])]
    kept, dropped = check_objections(objs)
    for o in kept:
        if o["kind"] == "factual":
            genes = o.get("genes") or ([o["location"].split(".", 1)[1]] if o["location"].startswith("targets.") else [])
            try:
                o["evidence_ids"] = add(run, "reviewer", int(o["pmid"]), genes, [o["quote"]], o.get("drug"))
            except SystemExit as exc:   # no gene node to attach the quote to
                o["evidence_ids"] = []
                o["note"] = str(exc)
        else:
            o["rule"] = "framing: reword, reorder or remove what the text already says; add no fact"
    out = {"run": run, "kept": kept, "discarded": dropped,
           "instruction": "For each kept objection, correct the text at `location` in narrative.json once: a factual "
                          "objection is fixed citing its evidence_ids (or by removing the claim it disputes); a framing "
                          "objection only rewords, reorders or removes. Then continue with the reading check."}
    (ROOT / "runs" / run / "revisions_expert.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                                 encoding="utf-8")
    return {"kept": len(kept), "factual_kept": sum(o["kind"] == "factual" for o in kept),
            "discarded": len(dropped), "file": f"runs/{run}/revisions_expert.json"}


# ------------------------------------------------------------------ measurement with seeded errors
# Errors of the kinds the manual reviews of 2026-11-26 found, written into copies of real cards. Each applies to the
# first front-page item whose text matches; patterns that match nothing in a run are skipped.
SEEDS = [
    ("subunit", "factual", r"IL12B is the p40 subunit", "IL12B is the p19 subunit"),
    ("anti-p19 target", "factual", r"(anti-p19 antibodies mirikizumab, risankizumab and guselkumab, which bind) IL23A",
     r"\1 IL12B"),
    ("selectivity", "factual", r"selective JAK1 inhibitor upadacitinib", "selective JAK2 inhibitor upadacitinib"),
    ("trial outcome", "factual", r"did not show significant benefit in UC phase 2", "showed significant benefit in UC phase 3"),
    ("class as drug proposal", "framing", r"Positive control for the JAK class, not a proposal for a JAK2-selective drug\.",
     "A strong new candidate for a JAK2-selective drug in UC."),
    ("failed class as open", "factual", r"went no further than phase 2 in UC", "were approved in UC"),
    ("drug effect direction", "factual", r"lowered blood neutrophil counts", "raised blood neutrophil counts"),
    ("genetic direction", "framing", r"so inhibition runs against genetics", "so genetics supports inhibition"),
    ("overclaimed novelty", "framing", r"^Exploratory: the largest", "Validated target: the largest"),
    ("stage", "factual", r"nangibotide is recorded in phase 2", "nangibotide is approved for UC"),
    ("mechanism", "factual", r"which shows that blocking interleukin-23 alone is enough",
     "which shows that blocking interleukin-12 is required"),
    ("influx as tissue target", "framing", r"a tissue marker of neutrophil influx", "a T-cell-intrinsic drug target"),
]


def seed(run: str, seed_: int = 11) -> dict:
    """Two blind batches: each seeded card is in one batch, its unaltered original in the other."""
    from .frontpage import RULE_TEXT
    from .render import TARGET_GROUPS
    rd = ROOT / "runs" / run
    _RUN["dir"] = rd
    items = front_page(json.loads((rd / "narrative.json").read_text(encoding="utf-8")))
    seeded = []
    for name, kind, pat, rep in SEEDS:
        for it in items:
            if re.search(pat, it["text"]):
                new = re.sub(pat, rep, it["text"], count=1)
                inserted = re.sub(pat, rep, re.search(pat, it["text"]).group(0))
                seeded.append({**it, "text": new, "seed": name, "seed_kind": kind, "inserted": inserted,
                               "original_location": it["location"]})
                break
    rng = random.Random(seed_)
    rng.shuffle(seeded)
    halves = [seeded[0::2], seeded[1::2]]
    base = BASE / run / "seeded"
    header = {**HEADER, "selection_rule": RULE_TEXT, "blocks": [{"title": t, "note": n} for _, t, n in TARGET_GROUPS]}
    keys = {}
    for name, mine, other in (("A", halves[0], halves[1]), ("B", halves[1], halves[0])):
        taken = {s["original_location"] for s in mine}
        batch = [{**s, "location": f"seed:{s['seed']}"} for s in mine] + \
                [it for it in items if it["location"] not in taken]
        rng.shuffle(batch)
        k = _write([{kk: vv for kk, vv in b.items() if kk not in ("seed", "seed_kind", "inserted", "original_location")}
                    | {"location": b["location"]} for b in batch], base / name, header)
        keys[name] = {"key": k, "seeds": {s["seed"]: {"kind": s["seed_kind"], "inserted": s["inserted"],
                                                      "text": s["text"]} for s in mine}}
    (base / "seeds.json").write_text(json.dumps(keys, indent=1, ensure_ascii=False), encoding="utf-8")
    return {"seeded": len(seeded), "of": len(SEEDS), "batches": [str(base / n / "input.json") for n in ("A", "B")]}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s or "").strip().lower()


def _hits(o: dict, s: dict) -> bool:
    """The objection points at the seeded words: its quoted claim overlaps the inserted span of the seeded text."""
    text, ins, claim = _norm(s["text"]), _norm(s["inserted"]), _norm(o.get("claim"))
    if not claim:
        return False
    if ins in claim or claim in ins:
        return True
    i, j = text.find(ins), text.find(claim)
    return i >= 0 and j >= 0 and j < i + len(ins) and i < j + len(claim)


def score(run: str) -> dict:
    base = BASE / run / "seeded"
    meta = json.loads((base / "seeds.json").read_text(encoding="utf-8"))
    per_seed, clean = {}, {"items": 0, "objections": 0, "factual": 0, "factual_cited_validly": 0, "framing": 0,
                           "detail": []}
    for name in ("A", "B"):
        rp = base / name / "review.json"
        if not rp.exists():
            continue
        key, seeds = meta[name]["key"], meta[name]["seeds"]
        objs = json.loads(rp.read_text(encoding="utf-8")).get("objections", [])
        kept, dropped = check_objections(objs)
        valid = {id(o) for o in kept}
        by_id: dict[str, list] = {}
        for o in kept + dropped:
            by_id.setdefault(o.get("id"), []).append((o, o in kept))
        for oid, loc in key.items():
            if loc.startswith("seed:"):
                s = seeds[loc[5:]]
                os_ = by_id.get(oid, [])
                caught = [o for o, ok in os_ if _hits(o, s)]
                per_seed[loc[5:]] = {"kind": s["kind"], "caught": bool(caught),
                                     "caught_with_valid_objection": any(ok and _hits(o, s) for o, ok in os_),
                                     "objections_on_card": len(os_)}
            else:
                clean["items"] += 1
                for o, ok in by_id.get(oid, []):
                    clean["objections"] += 1
                    clean[o.get("kind") if o.get("kind") in KINDS else "framing"] += 1
                    if o.get("kind") == "factual" and ok:
                        clean["factual_cited_validly"] += 1
                    clean["detail"].append({"location": loc, "kind": o.get("kind"), "claim": o.get("claim"),
                                            "objection": o.get("objection"), "kept": ok,
                                            **({"reason": o["reason"]} if o.get("reason") else {})})
        del valid
    n = len(per_seed)
    out = {"run": run, "seeded": n,
           "caught": sum(v["caught"] for v in per_seed.values()),
           "caught_with_valid_objection": sum(v["caught_with_valid_objection"] for v in per_seed.values()),
           "by_kind": {k: {"seeded": sum(v["kind"] == k for v in per_seed.values()),
                           "caught": sum(v["caught"] for v in per_seed.values() if v["kind"] == k)} for k in KINDS},
           "per_seed": per_seed, "unaltered": clean}
    (ROOT / "evals" / "expert_review.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    return {k: v for k, v in out.items() if k not in ("per_seed", "unaltered")} | {
        "unaltered_items": clean["items"], "objections_on_unaltered": clean["objections"]}


def main(argv=None) -> None:
    import sys
    argv = argv or sys.argv[1:]
    cmd, run = argv[0], argv[1]
    print(json.dumps({"prepare": prepare, "validate": validate, "seed": seed, "score": score}[cmd](run), indent=1))


if __name__ == "__main__":
    main()
