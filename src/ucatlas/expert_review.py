"""Expert review of the summary and the orchestrator's decision, with every factual statement cited (added 2026-11-27;
per-candidate verdicts since 2026-11-28).

The verifier and the reading check (verify.py, claim_check.py) prove that each claim is in the evidence it cites. They
cannot see a claim that is true of its record but wrong or misleading to a UC clinician: an anti-p19 antibody called
an IL12B drug, a JAK2 card read as a JAK2-drug proposal, a failed drug class presented as open. Four independent
expert-style reviews of run 2026-11-23 found such problems by hand. This module makes that review a step of /informe:

1. prepare(run): what a reader sees, from RUN/narrative.json: the summary paragraphs (front page) and the
   orchestrator's decision (four targets per list of targets with no therapy in UC, decision.py), each card with its
   list, the gene's dossier and the orchestrator's rationale, numbers written out, in data/expert_review/<run>/
   input.json, with opaque ids (key.json maps them back).
2. The reviewer subagent (.claude/agents/revisor-experto.md) reads it as a UC clinician would, with PubMed search, and
   writes review.json: a verdict per decision card (support, concerns or against, with reasons; a reason that states a
   fact needs a PMID and a verbatim abstract sentence) and objections to the text. Two kinds of objection:
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


# ------------------------------------------------------------------ the summary and the selection as text
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


def _brief(dossier: dict) -> dict:
    """The dossier as the card shows it, without evidence ids."""
    keep = ("group", "direction", "stability", "cell_family", "meta_log2fc", "uc_papers", "lights",
            "proposed_experiment", "inhibitor", "risk_allele", "tier", "predicted_effect", "mendelian_randomisation",
            "biopsy_change", "cell_in_single_cells", "cell_matches_prediction", "drug", "key_papers")
    return {k: dossier[k] for k in keep if k in dossier}


def front_page(nar: dict, genes: list[str] | None = None) -> list[dict]:
    """Summary paragraphs and the decision cards, as rendered (genes: only these cards, for a second pass)."""
    from .decision import TITLE, pool
    from .graph import GraphStore
    store = GraphStore()
    out = [] if genes else [{"location": f"narrative.summary[{i}]", "kind": "summary paragraph", "text": plain(p, store)}
                            for i, p in enumerate(nar.get("summary", []))]
    dos = {r["gene"]: r for rows in pool().values() for r in rows}
    for tg in nar.get("targets", []):
        if genes and tg["gene"] not in genes:
            continue
        out.append({"location": f"targets.{tg['gene']}", "kind": "decision card", "list": TITLE.get(tg.get("list"), tg.get("list")),
                    "gene": tg["gene"], "role": tg.get("category"), "dossier": _brief(dos.get(tg["gene"], {})),
                    "text": plain(tg, store)})
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


HEADER = {"question": "What are the candidate therapeutic targets and biomarkers in ulcerative colitis, and how strong "
                      "is the evidence for each?",
          "decision_question": "Which new therapeutic targets for ulcerative colitis, not yet treated by any drug in UC, "
                               "have the strongest evidence to be tested first?",
          "what_you_read": "An ulcerative colitis evidence atlas: the summary of its findings (front page) and the "
                           "orchestrator's decision: four targets to test first from each list of targets with no therapy "
                           "in UC. Each decision card shows the list, the gene's dossier from the atlas data and the "
                           "orchestrator's rationale."}


def _header() -> dict:
    from .decision import LISTS, RULE_TEXT
    return {**HEADER, "selection_rule": RULE_TEXT,
            "blocks": [{"title": title, "note": {
                "first_in_class": "Genes that rise in UC as a change of the tissue, with few UC papers and no clinical drug for any disease.",
                "known": "Established UC biology (more than 15 UC papers) that no drug reaches.",
                "repurposing": "An inhibitor exists in the clinic for another disease; never tried in UC, Crohn's disease or IBD.",
                "genetic": "Inherited UC risk variants: the gene, cell and direction predicted from DNA sequence (AlphaGenome), checked against the rest of the atlas."}[k]}
                for k, title, *_ in LISTS]}


def prepare(run: str, genes: list[str] | None = None) -> dict:
    """genes: a second pass on these decision cards only (after the orchestrator replaces a gene), in <run>/pass2/."""
    rd = ROOT / "runs" / run
    _RUN["dir"] = rd
    nar = json.loads((rd / "narrative.json").read_text(encoding="utf-8"))
    items = front_page(nar, genes)
    base = BASE / run / "pass2" if genes else BASE / run
    _write(items, base, _header())
    return {"items": len(items), "input": str(base / "input.json"), "output": str(base / "review.json")}


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


VERDICTS = ("support", "concerns", "against")
REASON_KINDS = {"factual", "judgement"}


def check_verdicts(verdicts: list[dict], key: dict) -> tuple[dict, list[dict]]:
    """gene -> verdict with its kept reasons; and the reasons discarded (with why). A factual reason needs a PMID and
    a verbatim abstract sentence of a peer-reviewed, non-retracted paper; a judgement reason weighs facts the card
    already shows and adds none."""
    from .pubmed_quotes import fetch_many, problems
    pmids = sorted({int(r["pmid"]) for v in verdicts for r in v.get("reasons", [])
                    if r.get("kind") == "factual" and str(r.get("pmid", "")).isdigit()})
    recs = fetch_many(pmids) if pmids else {}
    out, dropped = {}, []
    for v in verdicts:
        loc = key.get(v.get("id"), "")
        if not loc.startswith("targets.") or v.get("verdict") not in VERDICTS:
            dropped.append({**v, "reason_dropped": "not a decision card or unknown verdict"})
            continue
        gene = loc.split(".", 1)[1]
        kept = []
        for r in v.get("reasons", []):
            kind = r.get("kind")
            if kind not in REASON_KINDS:
                dropped.append({"gene": gene, **r, "reason_dropped": f"unknown kind {kind!r}"})
            elif kind == "judgement":
                kept.append(r)
            elif not str(r.get("pmid", "")).isdigit() or not r.get("quote"):
                dropped.append({"gene": gene, **r, "reason_dropped": "factual reason without a PMID and a verbatim quote"})
            elif int(r["pmid"]) not in recs:
                dropped.append({"gene": gene, **r, "reason_dropped": f"PMID {r['pmid']} not found in PubMed"})
            else:
                pr = problems(recs[int(r["pmid"])], [r["quote"]])
                if pr:
                    dropped.append({"gene": gene, **r, "reason_dropped": "; ".join(pr)})
                else:
                    kept.append({**r, "title": recs[int(r["pmid"])].get("title"), "journal": recs[int(r["pmid"])].get("journal"),
                                 "year": recs[int(r["pmid"])].get("year")})
        out[gene] = {"verdict": v["verdict"], "summary": v.get("summary", ""), "reasons": kept}
    return out, dropped


def validate(run: str, second_pass: bool = False) -> dict:
    """Keep the objections that may be acted on and the verdicts' reasons that are cited, register their quotes as
    citable evidence, and write RUN/revisions_expert.json (objections, for the orchestrator) and
    RUN/expert_verdicts.json (verdicts, shown on the decision cards)."""
    from .pubmed_quotes import add
    base = BASE / run / "pass2" if second_pass else BASE / run
    key = json.loads((base / "key.json").read_text(encoding="utf-8"))
    review = json.loads((base / "review.json").read_text(encoding="utf-8"))
    vfile = ROOT / "runs" / run / "expert_verdicts.json"
    verdicts = json.loads(vfile.read_text(encoding="utf-8")) if second_pass and vfile.exists() else {"verdicts": {}, "discarded": []}
    new, vdrop = check_verdicts(review.get("verdicts", []), key)
    for gene, v in new.items():
        for r in v["reasons"]:
            if r["kind"] == "factual":
                try:
                    r["evidence_ids"] = add(run, "reviewer", int(r["pmid"]), r.get("genes") or [gene], [r["quote"]], r.get("drug"))
                except SystemExit as exc:
                    r["evidence_ids"], r["note"] = [], str(exc)
    verdicts["verdicts"].update(new)
    verdicts["discarded"] += vdrop
    vfile.write_text(json.dumps(verdicts, indent=1, ensure_ascii=False), encoding="utf-8")
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
    rf = ROOT / "runs" / run / ("revisions_expert_pass2.json" if second_pass else "revisions_expert.json")
    rf.write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    return {"kept": len(kept), "factual_kept": sum(o["kind"] == "factual" for o in kept),
            "discarded": len(dropped), "file": rf.relative_to(ROOT).as_posix(),
            "verdicts": {g: v["verdict"] for g, v in new.items()},
            "verdict_reasons_kept": sum(len(v["reasons"]) for v in new.values()), "verdict_reasons_discarded": len(vdrop)}


# ------------------------------------------------------------------ measurement with seeded errors
# Errors of the kinds the manual reviews of 2026-11-26 found, written into copies of real cards. Each applies to the
# first item whose text matches; patterns that match nothing in a run are skipped. Rewritten for the decision cards of
# 2026-11-28 (the earlier seeds targeted the IL12B and JAK2 cards, which left the decision).
SEEDS = [
    ("approval status", "factual", r"an inhibitor is approved for pulmonary arterial hypertension",
     "an inhibitor is approved for ulcerative colitis"),
    ("mouse direction", "factual", r"knockdown or ablation protected against colitis", "knockdown or ablation worsened colitis"),
    ("variant direction", "factual", r"implies its function protects", "implies its function drives the disease"),
    ("genetic direction", "factual", r"higher expression in sigmoid colon raises UC risk",
     "higher expression in sigmoid colon lowers UC risk"),
    ("knockout phenotype", "factual", r"knockout mice develop more severe colitis", "knockout mice are protected from colitis"),
    ("trial stage", "factual", r"an inhibitor reached phase 3 in sepsis and covid-19", "an inhibitor was approved for UC"),
    ("protein function", "factual", r"UniProt describes it as a tumour suppressor", "UniProt describes it as an oncogene"),
    ("influx as tissue target", "framing", r"so blocking it would act on neutrophil recruitment rather than on a change of the tissue",
     "so blocking it would correct a change of the tissue itself"),
    ("overclaimed validation", "framing", r"Direction: inhibit, as a hypothesis to test\.", "Direction: inhibit; a validated target."),
    ("genetics overclaimed", "framing", r"plasma Mendelian randomisation is null", "plasma Mendelian randomisation supports inhibition"),
    ("severity as drug-specific", "framing", r"so it marks severity\.", "so it predicts ustekinumab failure specifically."),
    ("druggability overclaimed", "framing", r"a transcription factor is hard to drug", "a transcription factor is easy to drug"),
]


def seed(run: str, seed_: int = 11) -> dict:
    """Two blind batches: each seeded card is in one batch, its unaltered original in the other."""
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
    header = _header()
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
        rv = json.loads(rp.read_text(encoding="utf-8"))
        objs = rv.get("objections", [])
        # a verdict reason that points at the card's words counts like an objection (2026-11-28)
        objs += [{"id": v.get("id"), "kind": "factual" if r.get("kind") == "factual" else "framing",
                  "claim": r.get("claim"), "objection": r.get("text"), "pmid": r.get("pmid"), "quote": r.get("quote")}
                 for v in rv.get("verdicts", []) for r in v.get("reasons", []) if r.get("claim")]
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
    if cmd == "prepare" and len(argv) > 2:      # prepare <run> GENE1,GENE2: second pass on replaced genes
        print(json.dumps(prepare(run, argv[2].split(",")), indent=1))
    elif cmd == "validate" and "--second-pass" in argv:
        print(json.dumps(validate(run, second_pass=True), indent=1))
    else:
        print(json.dumps({"prepare": prepare, "validate": validate, "seed": seed, "score": score}[cmd](run), indent=1))


if __name__ == "__main__":
    main()
