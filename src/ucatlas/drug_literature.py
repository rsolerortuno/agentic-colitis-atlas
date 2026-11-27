"""What a drug binds and how it did in its trials, from peer-reviewed papers (added 2026-11-26).

Open Targets lists a drug under every gene of the complex or family it acts on, without saying which subunit it binds
(anti-p19 antibodies appear under IL12B, upadacitinib under JAK1, JAK2, JAK3 and TYK2) and without the trial results.
An expert-style review of the top candidates found that the report therefore called anti-p19 antibodies IL12B drugs,
said nothing of JAK1 selectivity and nothing of the neutrophil drop CXCR2 antagonists cause. This layer gives, for a
small curated set of drugs, a verbatim sentence from the PubMed abstract of a peer-reviewed trial paper (phase 3 or
2b/3 in UC, or the trial that measured the effect stated), so statements can cite it like any other record.

Inputs: data/drug_literature.json (kept in the repository). check() re-fetches the abstracts from PubMed (NCBI
E-utilities) and confirms every quote is verbatim. load() writes one drug_literature record per paper, linked to each
gene listed for it (gene -> drug_trial_of -> PUBMED:<pmid>).
"""
from __future__ import annotations

import json
import re

from .paths import ROOT

DATA = ROOT / "data" / "drug_literature.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def records() -> list[dict]:
    return json.loads(DATA.read_text(encoding="utf-8"))["records"]


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def check() -> dict:
    """Every quote must be verbatim in the abstract PubMed returns for its PMID."""
    import requests
    ids = ",".join(str(r["pmid"]) for r in records())
    xml = requests.get(EUTILS, params={"db": "pubmed", "id": ids, "rettype": "abstract", "retmode": "xml"},
                       timeout=60).text
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml)
    text = {}
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//PMID")
        parts = ["".join(a.itertext()) for a in art.iter("AbstractText")]
        text[int(pmid)] = _norm(" ".join(parts))
    bad = [(r["key"], q[:60]) for r in records() for q in r["quotes"] if _norm(q) not in text.get(r["pmid"], "")]
    # synonyms: PubMed must map every code to the same chemical (its query translation names the concept)
    syn_bad = []
    for s in synonyms():
        concepts = []
        for name in s["names"][1:]:
            tr = requests.get(EUTILS.replace("efetch", "esearch"), params={"db": "pubmed", "term": name, "retmode": "json"},
                              timeout=60).json()["esearchresult"].get("querytranslation", "")
            concepts.append(re.findall(r'"([^"]+)"\[Supplementary Concept\]', tr))
        if not concepts or not set.intersection(*map(set, concepts)):
            syn_bad.append(s["key"])
    return {"papers": len(records()), "quotes": sum(len(r["quotes"]) for r in records()), "not_verbatim": bad,
            "synonym_sets": len(synonyms()), "synonyms_not_confirmed": syn_bad}


def synonyms() -> list[dict]:
    return json.loads(DATA.read_text(encoding="utf-8")).get("synonyms", [])


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    n = 0
    with store.batch():
        for r in records():
            pub = f"PUBMED:{r['pmid']}"
            store.add_node(pub, "Publication", r["title"], pmid=r["pmid"], doi=r["doi"])
            for g in r["genes"]:
                node = store.resolve(g)
                if not node or node["type"] != "Gene":
                    continue
                ev = store.add_evidence(evidence_id("drug_literature", r["key"], g), "drug_literature",
                                        f"PMID:{r['pmid']}", date.today().isoformat(), gene=g, pmid=r["pmid"],
                                        doi=r["doi"], journal=r["journal"], year=r["year"], title=r["title"],
                                        trial=r["trial"], drug=r["drug"], binds=r["binds"],
                                        sentence=" ".join(r["quotes"]))
                store.add_claim(node["id"], "drug_trial_of", pub, ev)
                n += 1
        for s in synonyms():
            src = f"MESH:{s['mesh_ui']}"
            store.add_node(src, "Dataset", f"NLM MeSH supplementary concept {s['mesh_ui']} ({s['drug']})")
            for g in s["genes"]:
                node = store.resolve(g)
                if not node or node["type"] != "Gene":
                    continue
                ev = store.add_evidence(evidence_id("drug_literature", s["key"], g), "drug_literature", src,
                                        date.today().isoformat(), gene=g, drug=s["drug"], names=s["names"],
                                        mesh_source=s["source"], sentence=s["sentence"])
                store.add_claim(node["id"], "drug_names", src, ev)
                n += 1
    return n


if __name__ == "__main__":
    print(check())
    print("records loaded:", load())
