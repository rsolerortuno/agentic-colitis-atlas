"""Verbatim PubMed quotes the orchestrator and the expert reviewer may add during a run (added 2026-11-27).

The literature layers are fixed before the agents run (the PubTator3 corpus, Europe PMC per candidate, the curated
drug trial papers). When the orchestrator or the reviewer needs a fact none of them holds (what a drug binds, how a
trial ended), it may search PubMed, but it cannot cite what it read. It registers the paper here instead:

    uv run python -m ucatlas.pubmed_quotes add --run 2026-11-27 --by reviewer --pmid 28845751 --genes IL1B \
        [--drug canakinumab] --quote "verbatim sentence from the abstract" [--quote ...]

add() fetches the PubMed record (NCBI E-utilities) and refuses the paper unless every quote is verbatim in the
abstract (whitespace normalised) and the record is a peer-reviewed journal article: not a preprint, comment,
editorial, letter, news item or erratum, and not retracted. Accepted quotes are stored in data/pubmed_quotes.json
(kept in the repository, so a rebuilt graph loads them again), written to the graph as pubmed_quote evidence
(gene -> quoted_in -> PUBMED:<pmid>) and added to the run's evidence basis, and the new evidence ids are printed for
citing. Statements citing them pass the same checks as any other.
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import date

from .paths import ROOT

DATA = ROOT / "data" / "pubmed_quotes.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NOT_PEER_REVIEWED = {"Preprint", "Comment", "Editorial", "Letter", "News", "Newspaper Article", "Published Erratum",
                     "Retracted Publication", "Retraction of Publication", "Expression of Concern"}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def _get(url: str, params: dict):
    """E-utilities allow three requests a second without a key: retry briefly when refused."""
    import time
    import requests
    for i in range(4):
        r = requests.get(url, params=params, timeout=60)
        if r.status_code == 200 and not r.text.lstrip().startswith('{"error'):
            return r
        time.sleep(1 + i)
    r.raise_for_status()
    raise RuntimeError(r.text[:200])


def _parse(art) -> dict:
    return {"pmid": int(art.findtext(".//PMID")),
            "title": "".join(art.find(".//ArticleTitle").itertext()),
            "journal": art.findtext(".//Journal/ISOAbbreviation") or art.findtext(".//Journal/Title"),
            "year": art.findtext(".//JournalIssue/PubDate/Year") or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4],
            "doi": next((i.text for i in art.iter("ArticleId") if i.get("IdType") == "doi"), None),
            "publication_types": [t.text for t in art.iter("PublicationType")],
            "retracted": any(c.get("RefType") == "RetractionIn" for c in art.iter("CommentsCorrections")),
            "abstract": _norm(" ".join("".join(t.itertext()) for t in art.iter("AbstractText")))}


def fetch_many(pmids: list[int]) -> dict[int, dict]:
    import xml.etree.ElementTree as ET
    if not pmids:
        return {}
    xml = _get(EUTILS, {"db": "pubmed", "id": ",".join(map(str, pmids)), "rettype": "abstract", "retmode": "xml"}).text
    return {r["pmid"]: r for r in map(_parse, ET.fromstring(xml).iter("PubmedArticle"))}


def fetch(pmid: int) -> dict:
    r = fetch_many([int(pmid)]).get(int(pmid))
    if r is None:
        raise ValueError(f"PMID {pmid} not found in PubMed")
    return r


def search(query: str, n: int = 10) -> list[dict]:
    """PubMed search: the records of the first n hits, with their abstracts and whether they may be cited."""
    ids = _get(EUTILS.replace("efetch", "esearch"), {"db": "pubmed", "term": query, "retmax": n,
                                                     "retmode": "json"}).json()["esearchresult"]["idlist"]
    recs = fetch_many([int(i) for i in ids])
    out = []
    for i in ids:
        r = recs.get(int(i))
        if r:
            pr = problems(r, [])
            out.append({**r, "citable": not pr, "problems": pr})
    return out


def problems(rec: dict, quotes: list[str]) -> list[str]:
    """Why a paper and its quotes may not be cited (empty when they may)."""
    out = []
    bad = sorted(set(rec["publication_types"]) & NOT_PEER_REVIEWED)
    if bad:
        out.append(f"not a peer-reviewed article ({', '.join(bad)})")
    if rec["retracted"]:
        out.append("retracted")
    if "Journal Article" not in rec["publication_types"] and not set(rec["publication_types"]) & {
            "Review", "Randomized Controlled Trial", "Clinical Trial", "Meta-Analysis", "Systematic Review"}:
        out.append(f"not a journal article ({', '.join(rec['publication_types'])})")
    if not rec["abstract"]:
        out.append("no abstract to check the quote against")
    out += [f"not verbatim in the abstract: {q[:80]}" for q in quotes if _norm(q) not in rec["abstract"]]
    return out


def entries() -> list[dict]:
    return json.loads(DATA.read_text(encoding="utf-8"))["entries"] if DATA.exists() else []


def _key(e: dict) -> str:
    import hashlib
    return hashlib.sha1(json.dumps([e["pmid"], e["quotes"]], ensure_ascii=False).encode()).hexdigest()[:10]


def _write_entry(store, e: dict) -> list[str]:
    from .graph import evidence_id
    pub = f"PUBMED:{e['pmid']}"
    store.add_node(pub, "Publication", e["title"], pmid=e["pmid"], doi=e["doi"])
    ids = []
    for g in e["genes"]:
        node = store.resolve(g)
        if not node or node["type"] != "Gene":
            continue
        ev = store.add_evidence(evidence_id("pubmed_quote", e["key"], g), "pubmed_quote", f"PMID:{e['pmid']}",
                                e["added"], gene=g, pmid=e["pmid"], doi=e["doi"], journal=e["journal"],
                                year=e["year"], title=e["title"], drug=e.get("drug"), added_by=e["added_by"],
                                run=e["run"], sentence=" ".join(e["quotes"]))
        store.add_claim(node["id"], "quoted_in", pub, ev)
        ids.append(ev)
    return ids


def add(run: str, by: str, pmid: int, genes: list[str], quotes: list[str], drug: str | None = None) -> list[str]:
    from .graph import GraphStore
    from .provenance import extend_basis
    rec = fetch(pmid)
    bad = problems(rec, quotes)
    if bad:
        raise SystemExit(f"PMID {pmid} refused: " + "; ".join(bad))
    e = {"pmid": int(pmid), "doi": rec["doi"], "journal": rec["journal"], "year": rec["year"], "title": rec["title"],
         "publication_types": rec["publication_types"], "genes": genes, "drug": drug, "quotes": quotes,
         "added_by": by, "run": run, "added": date.today().isoformat()}
    e["key"] = _key(e)
    store = GraphStore()
    with store.batch():
        ids = _write_entry(store, e)
    if not ids:
        raise SystemExit(f"none of {genes} is a gene node of the graph")
    all_e = [x for x in entries() if x["key"] != e["key"]] + [e]
    DATA.write_text(json.dumps({"about": __doc__.split("\n\n")[0], "entries": all_e}, indent=1, ensure_ascii=False),
                    encoding="utf-8")
    extend_basis(ROOT / "runs" / run, store, ids)
    return ids


def load(store=None) -> int:
    """Rebuild step: write every stored entry to the graph again."""
    from .graph import GraphStore
    store = store or GraphStore()
    n = 0
    with store.batch():
        for e in entries():
            n += len(_write_entry(store, e))
    return n


def check() -> dict:
    """Re-fetch every stored paper: quotes still verbatim, still peer-reviewed, not retracted since."""
    bad = {}
    recs = fetch_many(sorted({e["pmid"] for e in entries()}))
    for e in entries():
        p = problems(recs[e["pmid"]], e["quotes"]) if e["pmid"] in recs else ["not found in PubMed"]
        if p:
            bad[e["key"]] = p
    return {"entries": len(entries()), "problems": bad}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="python -m ucatlas.pubmed_quotes")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("add")
    a.add_argument("--run", required=True)
    a.add_argument("--by", required=True, choices=["orchestrator", "reviewer"])
    a.add_argument("--pmid", required=True, type=int)
    a.add_argument("--genes", required=True)
    a.add_argument("--drug")
    a.add_argument("--quote", action="append", required=True)
    sub.add_parser("check")
    sub.add_parser("load")
    s = sub.add_parser("search")
    s.add_argument("query")
    s.add_argument("-n", type=int, default=8)
    s = sub.add_parser("abstract")
    s.add_argument("pmid", type=int)
    args = ap.parse_args(argv)
    if args.cmd == "search":
        for r in search(args.query, args.n):
            print(json.dumps({k: r[k] for k in ("pmid", "title", "journal", "year", "citable")}, ensure_ascii=False))
            print("   ", r["abstract"][:600])
    elif args.cmd == "abstract":
        r = fetch(args.pmid)
        print(json.dumps({**r, "citable": not problems(r, []), "problems": problems(r, [])}, indent=1, ensure_ascii=False))
    elif args.cmd == "add":
        ids = add(args.run, args.by, args.pmid, [g.strip() for g in args.genes.split(",")], args.quote, args.drug)
        print(json.dumps({"evidence_ids": ids}))
    elif args.cmd == "check":
        print(json.dumps(check(), indent=1))
    else:
        print("loaded:", load())


if __name__ == "__main__":
    main()
