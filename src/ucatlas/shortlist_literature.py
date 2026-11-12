"""What the literature says about each shortlist candidate in relation to UC (added 2026-11-11).

The discovery screen counted UC papers per gene (Europe PMC full text) but nobody read them. This layer reads them.
1. retrieve(): for each first-in-class gene, Europe PMC records that mention the gene (symbol or an unambiguous HGNC
   alias) with ulcerative colitis anywhere in the text (the same search as the paper count), plus records with the gene
   in the title or abstract and colitis, inflammatory bowel disease or Crohn's disease. For each record the title and
   abstract, and for open-access papers whose abstract does not name the gene, up to four full-text sentences that do
   (Europe PMC full-text XML). Written as JSONL batches for the classifier.
2. The classifier subagent (.claude/agents/clasificador-literatura.md, Claude Code subscription) reads each record and
   writes one JSON line: whether the paper is about this gene, the kind of evidence (list mention, expression,
   biomarker, functional, genetic), the species or system, the direction (pathogenic: higher activity worsens
   inflammation, so inhibition should help; protective: the gene limits inflammation or helps healing, so inhibition
   may harm; neutral), and one sentence copied verbatim from the text given.
3. validate(): keeps lines whose sentence is verbatim in the record text (whitespace and case ignored), with allowed
   values; everything else is rejected with a reason.
4. summarize(): a paper is relevant when colitis, IBD, Crohn's disease or UC appears in its title, abstract or the
   full-text sentences given (deterministic; most retrieved papers only mention UC somewhere in the full text). Only
   relevant papers set the light; functional studies in other inflammatory settings are shown as context. Per gene,
   counts by kind and direction and a literature light for the shortlist:
   green: at least one functional study (human, mouse or cells) with a pathogenic direction and none protective;
   red: functional evidence of a protective role and none pathogenic (inhibiting it may harm);
   amber: a specific expression, biomarker or genetic result, or functional evidence in both directions;
   grey: no paper about the gene, or the gene is only named in lists of many genes or in reviews.
Two gene sets (2026-11-12): "shortlist" (the first-in-class genes, as above) and "repurposing" (the genes with an
inhibitor for another disease, never tried in UC; validation_shortlist.repurposing_genes). These genes have hundreds
or thousands of UC papers, so only the second search is used for them (gene in the title or abstract with colitis,
IBD or Crohn's disease), the 25 most cited per gene: functional studies are found there. Same classifier, same
validation and the same light rules; separate folders (data/literature/repurposing/), outputs
(repurposing_literature*.parquet) and evidence type (repurposing_literature).
Outputs: data/literature/shortlist/records.parquet and batches/ (abstracts and full-text sentences: not in the
repository, rebuilt by retrieve() from Europe PMC), classified/ (agent output, kept in the repository),
data/processed/shortlist_literature.parquet, shortlist_literature_rejected.json.
"""
from __future__ import annotations

import json
import re
import time
import xml.etree.ElementTree as ET

import pandas as pd

from .paths import LITERATURE, PROCESSED, RAW

DIR = LITERATURE / "shortlist"
BATCHES = DIR / "batches"
CLASSIFIED = DIR / "classified"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
MAX_RECORDS = 60
BATCH_SIZE = 25
# aliases that name something else (another gene, a drug, a common word) are not searched
AMBIGUOUS = {"MIG2", "LAX", "UNG2", "UDG2", "CCNU", "PP20", "DOC-1", "CPX1", "CPX-1", "CLEC1", "PRO0082", "CHGN",
             "INAM", "BNF1"}
CLONE = re.compile(r"^(FLJ|MGC|KIAA|DKFZ|dJ|RP\d|OTTHUMP|PRO\d|HsT|C\d+orf)", re.I)
IBD = re.compile(r"colitis|inflammatory bowel|crohn|\bIBD\b|\bUC\b", re.I)
KINDS = {"not_about_gene", "list_mention", "expression", "biomarker", "functional", "genetic", "review"}
DIRECTIONS = {"pathogenic", "protective", "neutral", "unclear"}
SYSTEMS = {"human_tissue", "human_blood", "mouse", "cells", "organoids", "in_silico", "other"}


SETS = {
    "shortlist": {"dir": LITERATURE / "shortlist", "out": "shortlist_literature", "max_records": 60,
                  "queries": ("uc_any_text", "ibd_title_abstract"),
                  "dataset": ("DATASET:SHORTLIST_LIT", "Literature review of the shortlist candidates (Europe PMC, "
                                                       "classified by a subagent, verbatim-checked)")},
    "repurposing": {"dir": LITERATURE / "repurposing", "out": "repurposing_literature", "max_records": 25,
                    "queries": ("ibd_title_abstract",),
                    "dataset": ("DATASET:REPURPOSING_LIT", "Literature review of the repurposing candidates (Europe "
                                                           "PMC, classified by a subagent, verbatim-checked)")},
}


def _genes(which: str) -> list[str]:
    if which == "repurposing":
        from .validation_shortlist import repurposing_genes
        return repurposing_genes()
    return first_in_class()


def first_in_class() -> list[str]:
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    return d.gene[d.status == "first-in-class opportunity"].tolist()


def names(gene: str) -> list[str]:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", dtype=str,
                     usecols=["symbol", "alias_symbol", "prev_symbol"]).set_index("symbol")
    out = [gene]
    if gene in hg.index:
        for col in ("alias_symbol", "prev_symbol"):
            v = hg.at[gene, col]
            if isinstance(v, str):
                out += [a for a in v.split("|") if len(a) >= 3 and not CLONE.match(a) and a.upper() not in AMBIGUOUS]
    return list(dict.fromkeys(out))


def _get(url: str, params: dict | None = None, kind: str = "json"):
    """Cached GET; a missing full text (404) is cached as None so it is not asked again."""
    import requests
    from .external import _cached

    def fetch():
        for attempt in range(4):
            try:
                r = requests.get(url, params=params, timeout=40)
            except requests.RequestException:
                time.sleep(2 * (attempt + 1))
                continue
            if r.status_code == 404:
                return None
            if r.status_code == 200:
                return r.json() if kind == "json" else r.text
            time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"Europe PMC failed: {url}")
    return _cached("europepmc_shortlist", url + json.dumps(params or {}, sort_keys=True), fetch)["data"]


def _search(query: str, max_records: int = MAX_RECORDS) -> list[dict]:
    out, cursor = [], "*"
    while len(out) < max_records:
        d = _get(f"{EPMC}/search", {"query": query, "format": "json", "pageSize": 100, "resultType": "core",
                                    "cursorMark": cursor, "sort": "CITED desc"})
        res = d.get("resultList", {}).get("result", [])
        out += res
        if not res or d.get("nextCursorMark") in (None, cursor):
            break
        cursor = d["nextCursorMark"]
    return out[:max_records]


def _token(nm: list[str]) -> re.Pattern:
    return re.compile(r"(?<![A-Za-z0-9])(" + "|".join(re.escape(n) for n in nm) + r")(?![A-Za-z0-9])", re.I)


def _fulltext_sentences(pmcid: str, tok: re.Pattern, k: int = 4) -> list[str]:
    try:
        xml = _get(f"{EPMC}/{pmcid}/fullTextXML", kind="text")
        if not xml:
            return []
        root = ET.fromstring(xml)
    except (ET.ParseError, RuntimeError):   # malformed XML or the service kept failing
        return []
    text = " ".join(" ".join(p.itertext()) for p in root.iter("p"))
    sents = re.split(r"(?<=[.!?])\s+(?=[A-Z])", re.sub(r"\s+", " ", text))
    hits = [s.strip() for s in sents if tok.search(s) and 30 <= len(s) <= 600]
    return hits[:k]


def retrieve(genes: list[str] | None = None, which: str = "shortlist") -> dict:
    cfg = SETS[which]
    if genes is None:
        genes = _genes(which)
    from concurrent.futures import ThreadPoolExecutor
    from functools import partial
    with ThreadPoolExecutor(4) as ex:   # Europe PMC answers slowly; four genes at a time
        per_gene = list(ex.map(partial(_gene_records, queries=cfg["queries"], max_records=cfg["max_records"]), genes))
    rows = [r for rs in per_gene for r in rs]
    d = pd.DataFrame(rows)
    return _write(d, genes, cfg["dir"])


def _gene_records(g: str, queries: tuple = ("uc_any_text", "ibd_title_abstract"), max_records: int = MAX_RECORDS) -> list[dict]:
    from .external import DISEASE_TERMS
    rows = []
    if True:
        nm = names(g)
        tok = _token(nm)
        q_names = " OR ".join(f'"{n}"' for n in nm)
        q_abs = " OR ".join(f'ABSTRACT:"{n}" OR TITLE:"{n}"' for n in nm)
        found = {}
        for q, why in ((f"({q_names}) AND {DISEASE_TERMS}", "uc_any_text"),
                       (f'({q_abs}) AND (colitis OR "inflammatory bowel" OR crohn*)', "ibd_title_abstract")):
            if why not in queries:
                continue
            for r in _search(q, max_records):
                key = r.get("pmid") or r.get("id")
                if key not in found:
                    found[key] = (r, why)
        for key, (r, why) in found.items():
            title, abstract = r.get("title") or "", re.sub(r"<[^>]+>", "", r.get("abstractText") or "")
            in_abs = bool(tok.search(title + " " + abstract))
            snippets = [] if in_abs or not r.get("pmcid") else _fulltext_sentences(r["pmcid"], tok)
            rows.append({"gene": g, "id": str(key), "pmid": r.get("pmid"), "pmcid": r.get("pmcid"),
                         "year": int(r["pubYear"]) if str(r.get("pubYear", "")).isdigit() else None,
                         "title": title, "abstract": abstract, "fulltext_sentences": snippets,
                         "gene_in_abstract": in_abs, "found_by": why, "names_searched": nm,
                         "pub_type": "; ".join((r.get("pubTypeList") or {}).get("pubType", []))[:120]})
    return rows


def _write(d: pd.DataFrame, genes: list[str], base=DIR) -> dict:
    base.mkdir(parents=True, exist_ok=True)
    batches = base / "batches"
    d.to_parquet(base / "records.parquet", index=False)
    # only records where the gene is visible in the text given to the classifier
    usable = d[d.gene_in_abstract | (d.fulltext_sentences.map(len) > 0)].reset_index(drop=True)
    batches.mkdir(parents=True, exist_ok=True)
    for f in batches.glob("batch_*.jsonl"):
        f.unlink()
    for i in range(0, len(usable), BATCH_SIZE):
        with open(batches / f"batch_{i // BATCH_SIZE:03d}.jsonl", "w", encoding="utf-8") as fh:
            for r in usable.iloc[i:i + BATCH_SIZE].itertuples():
                fh.write(json.dumps({"gene": r.gene, "id": r.id, "names": list(r.names_searched), "title": r.title,
                                     "abstract": r.abstract, "fulltext_sentences": list(r.fulltext_sentences)},
                                    ensure_ascii=False) + "\n")
    return {"genes": len(genes), "records": len(d), "usable": len(usable),
            "gene_only_in_inaccessible_full_text": int(len(d) - len(usable)),
            "batches": len(list(batches.glob("batch_*.jsonl")))}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip().lower()


def validate(which: str = "shortlist") -> dict:
    cfg = SETS[which]
    rec = pd.read_parquet(cfg["dir"] / "records.parquet")
    text = {(r.gene, r.id): _norm(" ".join([r.title, r.abstract] + list(r.fulltext_sentences))) for r in rec.itertuples()}
    kept, rejected = [], []
    for f in sorted((cfg["dir"] / "classified").glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                rejected.append({"reason": "invalid_json", "raw": line[:200]})
                continue
            k = (r.get("gene"), str(r.get("id")))
            reason = None
            if k not in text:
                reason = "record_not_retrieved"
            elif r.get("kind") not in KINDS:
                reason = "kind_not_allowed"
            elif r.get("kind") != "not_about_gene" and r.get("direction") not in DIRECTIONS:
                reason = "direction_not_allowed"
            elif r.get("kind") != "not_about_gene" and r.get("system") not in SYSTEMS:
                reason = "system_not_allowed"
            elif r.get("kind") != "not_about_gene" and (not r.get("quote") or _norm(r["quote"]) not in text[k]):
                reason = "quote_not_in_text"
            if reason:
                rejected.append({"reason": reason, **{x: r.get(x) for x in ("gene", "id", "kind", "quote")}})
            else:
                kept.append({"gene": k[0], "id": k[1], "kind": r["kind"], "direction": r.get("direction", "unclear"),
                             "system": r.get("system", "other"), "quote": r.get("quote", ""),
                             "summary": (r.get("summary") or "")[:300]})
    k = pd.DataFrame(kept).drop_duplicates(["gene", "id"], keep="last")
    k = k.merge(rec[["gene", "id", "pmid", "pmcid", "year", "title"]], on=["gene", "id"], how="left")
    k.to_parquet(PROCESSED / f"{cfg['out']}_papers.parquet", index=False)
    (PROCESSED / f"{cfg['out']}_rejected.json").write_text(json.dumps(rejected, indent=1), encoding="utf-8")
    return {"kept": len(k), "rejected": len(rejected),
            "rejected_by_reason": pd.Series([x["reason"] for x in rejected]).value_counts().to_dict() if rejected else {}}


def summarize(which: str = "shortlist") -> dict:
    cfg = SETS[which]
    p = pd.read_parquet(PROCESSED / f"{cfg['out']}_papers.parquet")
    rec = pd.read_parquet(cfg["dir"] / "records.parquet")
    txt = {(r.gene, r.id): " ".join([r.title, r.abstract] + list(r.fulltext_sentences)) for r in rec.itertuples()}
    p["relevant"] = [bool(IBD.search(txt[(g, i)])) for g, i in zip(p.gene, p.id)]
    p.to_parquet(PROCESSED / f"{cfg['out']}_papers.parquet", index=False)
    genes = _genes(which)
    rows = []
    for g in genes:
        allx = p[(p.gene == g) & (p.kind != "not_about_gene")]
        x = allx[allx.relevant]
        func = x[x.kind == "functional"]
        ctx = allx[~allx.relevant & (allx.kind == "functional") & allx.direction.isin(["pathogenic", "protective"])]
        context = "; ".join(f"outside IBD: {r.direction} in {r.system.replace('_', ' ')} (PMID {r.pmid})" for r in ctx.itertuples())
        n_path, n_prot = int((func.direction == "pathogenic").sum()), int((func.direction == "protective").sum())
        if x.empty:
            light, verdict = "grey", "no paper about the gene in UC or IBD"
        elif x.kind.isin(["list_mention", "review"]).all():
            light, verdict = "grey", f"only named in lists or reviews ({len(x)} paper{'s' if len(x) > 1 else ''}): nothing specific"
        elif n_path and not n_prot:
            light, verdict = "green", f"{n_path} functional stud{'y' if n_path == 1 else 'ies'} with a pathogenic role"
        elif n_prot and not n_path:
            light, verdict = "red", f"{n_prot} functional stud{'y' if n_prot == 1 else 'ies'} with a protective role: inhibition may harm"
        elif n_prot and n_path:
            light, verdict = "amber", f"functional studies in both directions ({n_path} pathogenic, {n_prot} protective)"
        else:
            kinds = x.kind.value_counts()
            verdict = "only " + ", ".join(f"{v} {k.replace('_', ' ')}" for k, v in kinds.items())
            light = "amber"
        if context:
            verdict += f" ({context})"
        mouse = func[(func.system == "mouse") & func.direction.isin(["pathogenic", "protective"])]
        key = pd.concat([func, x[x.kind.isin(["biomarker", "genetic"])], x[x.kind == "expression"], ctx]).drop_duplicates("id")
        spec = pd.concat([x[~x.kind.isin(["list_mention", "review"])], ctx]).drop_duplicates("id")
        findings = [{"pmid": r.pmid, "pmcid": r.pmcid, "year": None if pd.isna(r.year) else int(r.year), "kind": r.kind,
                     "direction": r.direction, "system": r.system, "summary": r.summary, "ibd_setting": bool(r.relevant)}
                    for r in spec.sort_values(["relevant", "year"], ascending=False).itertuples()]
        rows.append({"gene": g, "light": light, "verdict": verdict, "papers_any_context": len(allx),
                     "list_or_review": int(x.kind.isin(["list_mention", "review"]).sum()),
                     "findings": json.dumps(findings, ensure_ascii=False),
                     "mouse_colitis": "; ".join(f"{r.direction} (PMID {r.pmid})" for r in mouse.itertuples()),
                     "records_retrieved": int((rec.gene == g).sum()),
                     "records_readable": int(((rec.gene == g) & (rec.gene_in_abstract | (rec.fulltext_sentences.map(len) > 0))).sum()),
                     "about_gene": len(x), "functional": len(func), "pathogenic": n_path, "protective": n_prot,
                     "by_kind": json.dumps(x.kind.value_counts().to_dict()),
                     "key_papers": json.dumps([{"pmid": r.pmid, "pmcid": r.pmcid, "year": None if pd.isna(r.year) else int(r.year),
                                                "kind": r.kind, "direction": r.direction, "system": r.system,
                                                "quote": r.quote, "summary": r.summary} for r in key.head(4).itertuples()],
                                              ensure_ascii=False)})
    d = pd.DataFrame(rows)
    d.to_parquet(PROCESSED / f"{cfg['out']}.parquet", index=False)
    return {"genes": len(d), "lights": d.light.value_counts().to_dict(),
            "green": d.gene[d.light == "green"].tolist(), "red": d.gene[d.light == "red"].tolist()}


def load(store=None, which: str = "shortlist") -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    cfg = SETS[which]
    d = pd.read_parquet(PROCESSED / f"{cfg['out']}.parquet")
    n = 0
    with store.batch():
        store.add_node(cfg["dataset"][0], "Dataset", cfg["dataset"][1])
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            ev = store.add_evidence(evidence_id(cfg["out"], r.gene), cfg["out"],
                                    cfg["dataset"][0], date.today().isoformat(), gene=r.gene, light=r.light,
                                    verdict=r.verdict, papers_about_gene=int(r.about_gene), functional=int(r.functional),
                                    pathogenic=int(r.pathogenic), protective=int(r.protective),
                                    key_papers=json.loads(r.key_papers))
            store.add_claim(node["id"], "has_literature_review_in", cfg["dataset"][0], ev)
            n += 1
    return n


if __name__ == "__main__":
    import sys
    step = sys.argv[1] if len(sys.argv) > 1 else "retrieve"
    which = sys.argv[2] if len(sys.argv) > 2 else "shortlist"
    fn = {"retrieve": retrieve, "validate": validate, "summarize": summarize}[step]
    print(json.dumps(fn(which=which), indent=1))
