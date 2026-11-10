"""M3: literature ingestion from PubTator3 (normalized annotations) + extraction batches + validation.

Relation extraction itself is done by the `extractor-relaciones` Claude Code subagent on the batches
written here; `validate_extractions()` keeps only relations whose literal sentence appears in the abstract.
"""
from __future__ import annotations

import json
import re
import time

import pandas as pd
import requests

from .config import load
from .paths import LITERATURE, PROCESSED, RAW

API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
PUBTATOR_DIR = LITERATURE / "pubtator"
BATCH_DIR = LITERATURE / "extraction_batches"
EXTRACTED_DIR = LITERATURE / "extracted"
RELATIONS = ["upregulated_in", "downregulated_in", "biomarker_of", "target_of", "drives", "associated_with", "inhibits"]
for _d in (PUBTATOR_DIR, BATCH_DIR, EXTRACTED_DIR):
    _d.mkdir(parents=True, exist_ok=True)


def _get(url: str, params: dict | None = None) -> requests.Response:
    for attempt in range(5):
        r = requests.get(url, params=params, timeout=60)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
        time.sleep(0.35)  # PubTator3 asks for <= 3 requests/second
        return r
    r.raise_for_status()
    return r


def search(text: str, n: int) -> list[int]:
    cache = PUBTATOR_DIR / f"search_{re.sub(r'[^A-Za-z0-9]+', '_', text)[:80]}_{n}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    pmids: list[int] = []
    page = 1
    while len(pmids) < n:
        res = _get(f"{API}/search/", {"text": text, "page": page}).json()
        hits = [int(h["pmid"]) for h in res.get("results", [])]
        if not hits:
            break
        pmids.extend(p for p in hits if p not in pmids)
        if page >= res.get("total_pages", page):
            break
        page += 1
    pmids = pmids[:n]
    cache.write_text(json.dumps(pmids))
    return pmids


def fetch_biocjson(pmids: list[int]) -> None:
    missing = [p for p in pmids if not (PUBTATOR_DIR / f"{p}.json").exists()]
    for i in range(0, len(missing), 100):
        chunk = missing[i:i + 100]
        data = _get(f"{API}/publications/export/biocjson", {"pmids": ",".join(map(str, chunk))}).json()
        for doc in data.get("PubTator3", []):
            (PUBTATOR_DIR / f"{doc['id']}.json").write_text(json.dumps(doc))


def _parse(doc: dict) -> tuple[dict, list[dict]]:
    pub = {"pmid": int(doc["id"]), "title": "", "abstract": "", "year": None, "journal": ""}
    anns = []
    for p in doc.get("passages", []):
        ptype = p.get("infons", {}).get("type")
        if ptype == "title":
            pub["title"] = p.get("text", "")
            pub["journal"] = p["infons"].get("journal", "")
            y = p["infons"].get("year")
            pub["year"] = int(y) if y and str(y).isdigit() else None
        elif ptype == "abstract":
            pub["abstract"] = (pub["abstract"] + " " + p.get("text", "")).strip()
        for a in p.get("annotations", []):
            inf = a.get("infons", {})
            anns.append({"pmid": pub["pmid"], "type": inf.get("type"), "identifier": inf.get("identifier"),
                         "name": inf.get("name"), "text": a.get("text"), "passage": ptype})
    return pub, anns


SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9(])")


def _relations(doc: dict) -> list[dict]:
    """PubTator3's own machine-extracted relations (entity pairs with a type and a score), kept apart from the
    verbatim-checked relations of the local extractor."""
    out = []
    for r in doc.get("relations") or []:
        inf = r.get("infons", {})
        a, b = inf.get("role1") or {}, inf.get("role2") or {}
        out.append({"pmid": int(doc["id"]), "relation_type": inf.get("type"), "score": float(inf.get("score") or 0),
                    "type1": a.get("type"), "id1": str(a.get("identifier")), "name1": a.get("name"),
                    "type2": b.get("type"), "id2": str(b.get("identifier")), "name2": b.get("name")})
    return out


def build_tables() -> dict:
    pubs, anns, chunks, rels = [], [], [], []
    for f in sorted(PUBTATOR_DIR.glob("[0-9]*.json")):
        doc = json.loads(f.read_text())
        pub, a = _parse(doc)
        rels.extend(_relations(doc))
        pubs.append(pub)
        anns.extend(a)
        text = (pub["title"] + " " + pub["abstract"]).strip()
        for i, s in enumerate(SENT_SPLIT.split(text)):
            if len(s) > 20:
                chunks.append({"chunk_id": f"PMID{pub['pmid']}_S{i}", "pmid": pub["pmid"], "text": s})
    pd.DataFrame(pubs).to_parquet(LITERATURE / "publications.parquet", index=False)
    pd.DataFrame(anns).to_parquet(LITERATURE / "annotations.parquet", index=False)
    pd.DataFrame(chunks).to_parquet(LITERATURE / "chunks.parquet", index=False)
    pd.DataFrame(rels).to_parquet(LITERATURE / "pubtator_relations.parquet", index=False)
    return {"publications": len(pubs), "annotations": len(anns), "chunks": len(chunks), "pubtator_relations": len(rels)}


def top_expression_genes(n: int) -> list[str]:
    """The n strongest meta-analysis changes (half up, half down) that are significant in the meta-analysis and in at
    least two cohorts in the same direction: the genes whose literature is searched individually."""
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    up = m[m["meta_log2fc"] > 0].nlargest(n - n // 2, "meta_log2fc")["gene"].tolist()
    down = m[m["meta_log2fc"] < 0].nsmallest(n // 2, "meta_log2fc")["gene"].tolist()
    return up + down


def ingest(refresh: bool = False) -> dict:
    """Fetch PubTator3 annotations. The PMID lists are frozen in query_log.json (versioned) so that a fresh
    clone reproduces the same corpus; pass refresh=True to search PubTator3 again."""
    log_path = LITERATURE / "query_log.json"
    if log_path.exists() and not refresh:
        log = json.loads(log_path.read_text())
        all_pmids = sorted(set(log["disease_pmids"]) | {p for v in log["gene_pmids"].values() for p in v}
                           | {p for v in log.get("extra_gene_pmids", {}).values() for p in v})
        fetch_biocjson(all_pmids)
        stats = build_tables()
        stats["top_genes"] = len(log["top_genes"])
        stats["source"] = "frozen query_log.json"
        return stats
    cfg = load()["literature"]
    genes = top_expression_genes(cfg["top_genes"])
    disease_pmids = search(cfg["disease_query"], cfg["disease_pmids"])
    # literature may use an older symbol (e.g. IL8 for CXCL8): use the platform-era symbol when it differs
    symmap = pd.read_parquet(PROCESSED / "symbol_map.parquet").dropna(subset=["symbol"])
    symmap = symmap[symmap["dataset_symbol"].str.upper() != symmap["symbol"].str.upper()]
    lit_symbol = symmap.drop_duplicates("symbol").set_index("symbol")["dataset_symbol"].to_dict()
    gene_pmids: dict[str, list[int]] = {}
    for g in genes:
        gene_pmids[g] = search(f"{lit_symbol.get(g, g)} {cfg['disease_query']}", cfg["pmids_per_gene"])
    all_pmids = sorted(set(disease_pmids) | {p for v in gene_pmids.values() for p in v})
    fetch_biocjson(all_pmids)
    (LITERATURE / "query_log.json").write_text(json.dumps(
        {"top_genes": genes, "disease_pmids": disease_pmids, "gene_pmids": gene_pmids}, indent=1))
    stats = build_tables()
    stats["top_genes"] = len(genes)
    return stats


def extend_corpus(genes: list[str], n: int = 10) -> dict:
    """Add up to n PubTator3 hits for '<gene> ulcerative colitis' for genes outside the original top-gene list (report
    candidates, convergence top 40, first-in-class genes; added 2026-11-10, third round). The PMID lists are frozen in
    query_log.json under extra_gene_pmids, so a rebuild fetches the same corpus. No new LLM extraction: the extra
    abstracts contribute PubTator3 annotations and relations and BM25 sentences."""
    log_path = LITERATURE / "query_log.json"
    log = json.loads(log_path.read_text())
    extra = log.setdefault("extra_gene_pmids", {})
    cfg = load()["literature"]
    for g in genes:
        if g in log["gene_pmids"] or g in extra:
            continue
        extra[g] = search(f"{g} {cfg['disease_query']}", n)
    log_path.write_text(json.dumps(log, indent=1))
    fetch_biocjson(sorted({p for v in extra.values() for p in v}))
    return {"extra_genes": len(extra), "extra_pmids": len({p for v in extra.values() for p in v}), **build_tables()}


UC_MESH = {"MESH:D003093", "MESH:D015212", "MESH:D003092"}   # ulcerative colitis, IBD, colitis


def load_pubtator_relations(store=None) -> int:
    """One pubtator_relations record per gene: PubTator3 relations in the local corpus that involve the gene, with the
    relation types to UC/IBD/colitis and the most frequent partners (machine-extracted; not verbatim-checked)."""
    from collections import Counter
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    rel = pd.read_parquet(LITERATURE / "pubtator_relations.parquet")
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t",
                     usecols=["symbol", "entrez_id"], dtype=str).dropna()
    ent = dict(zip(hg.entrez_id, hg.symbol))
    rows = []
    for side, other in (("1", "2"), ("2", "1")):
        d = rel[rel[f"type{side}"] == "Gene"].copy()
        d["gene"] = d[f"id{side}"].map(ent)
        d["partner"] = d[f"name{other}"]
        d["partner_id"] = d[f"id{other}"]
        d["partner_type"] = d[f"type{other}"]
        rows.append(d.dropna(subset=["gene"])[["gene", "pmid", "relation_type", "score", "partner", "partner_id", "partner_type"]])
    g = pd.concat(rows, ignore_index=True)
    n = 0
    with store.batch():
        store.add_node("DATASET:PUBTATOR3_RELATIONS", "Dataset", "PubTator3 machine-extracted relations in the local corpus")
        for gene, d in g.groupby("gene"):
            node = store.resolve(gene)
            if not node or node["type"] != "Gene":
                continue
            uc = d[d.partner_id.isin(UC_MESH)]
            rec = {"gene": gene, "relations": int(len(d)), "publications": int(d.pmid.nunique()),
                   "with_uc_or_ibd": {k: int(v) for k, v in Counter(uc.relation_type).most_common()},
                   "uc_or_ibd_pmids": sorted(set(int(x) for x in uc.pmid))[:15],
                   "top_partners": [{"partner": p, "type": t, "n": int(c)} for (p, t), c in
                                    Counter(zip(d.partner, d.partner_type)).most_common(8)],
                   "method": "PubTator3 relation extraction (machine), not verbatim-checked"}
            ev = store.add_evidence(evidence_id("pubtator_relations", gene), "pubtator_relations",
                                    "DATASET:PUBTATOR3_RELATIONS", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "mentioned_with", "DATASET:PUBTATOR3_RELATIONS", ev)
            n += 1
    return n


def write_extraction_batches() -> list[str]:
    """Pick <= max_extraction_abstracts abstracts: round-robin over gene-specific hits, then disease hits."""
    cfg = load()["literature"]
    log = json.loads((LITERATURE / "query_log.json").read_text())
    pubs = pd.read_parquet(LITERATURE / "publications.parquet").set_index("pmid")
    has_abs = set(pubs.index[pubs["abstract"].str.len() > 100])
    chosen: list[int] = []
    queues = [list(v) for v in log["gene_pmids"].values()] + [list(log["disease_pmids"])]
    while len(chosen) < cfg["max_extraction_abstracts"] and any(queues):
        for q in queues:
            while q:
                p = q.pop(0)
                if p in has_abs and p not in chosen:
                    chosen.append(p)
                    break
            if len(chosen) >= cfg["max_extraction_abstracts"]:
                break
    for old in BATCH_DIR.glob("batch_*.jsonl"):
        old.unlink()
    paths = []
    bs = cfg["extraction_batch_size"]
    for i in range(0, len(chosen), bs):
        path = BATCH_DIR / f"batch_{i // bs + 1:02d}.jsonl"
        with open(path, "w", encoding="utf-8") as f:
            for p in chosen[i:i + bs]:
                f.write(json.dumps({"pmid": int(p), "title": pubs.loc[p, "title"],
                                    "abstract": pubs.loc[p, "abstract"]}, ensure_ascii=False) + "\n")
        paths.append(str(path))
    return paths


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip().lower()


def validate_extractions() -> dict:
    """Keep relations whose literal sentence is verbatim (whitespace/case-insensitive) in the abstract."""
    pubs = pd.read_parquet(LITERATURE / "publications.parquet").set_index("pmid")
    batch_pmids = set()
    for b in BATCH_DIR.glob("batch_*.jsonl"):
        batch_pmids |= {json.loads(line)["pmid"] for line in b.read_text(encoding="utf-8").splitlines() if line}
    kept, rejected = [], []
    for f in sorted(EXTRACTED_DIR.glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                rejected.append({"reason": "invalid_json", "raw": line[:200]})
                continue
            reason = None
            pmid = r.get("pmid")
            if not isinstance(pmid, int) or pmid not in batch_pmids or pmid not in pubs.index:
                reason = "pmid_not_in_batch"
            elif r.get("relation") not in RELATIONS:
                reason = "relation_not_allowed"
            elif not r.get("subject") or not r.get("object"):
                reason = "missing_subject_or_object"
            elif not isinstance(r.get("confidence"), (int, float)) or not 0 <= r["confidence"] <= 1:
                reason = "bad_confidence"
            else:
                text = _norm(pubs.loc[pmid, "title"] + " " + pubs.loc[pmid, "abstract"])
                if not r.get("literal_sentence") or _norm(r["literal_sentence"]) not in text:
                    reason = "sentence_not_in_abstract"
            if reason:
                rejected.append({"reason": reason, **r})
            else:
                kept.append(r)
    df = pd.DataFrame(kept)
    if not df.empty:
        df = df.drop_duplicates(["pmid", "subject", "relation", "object"])
        df.insert(0, "relation_id", [f"REL{i:05d}" for i in range(len(df))])
    df.to_parquet(LITERATURE / "relations.parquet", index=False)
    (LITERATURE / "relations_rejected.json").write_text(json.dumps(rejected, indent=1, ensure_ascii=False))
    return {"kept": int(len(df)), "rejected": len(rejected),
            "reasons": pd.Series([r["reason"] for r in rejected]).value_counts().to_dict() if rejected else {}}


if __name__ == "__main__":
    print(json.dumps(ingest(), indent=1))
    print(write_extraction_batches())
