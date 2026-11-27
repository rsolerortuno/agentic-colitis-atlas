"""Peer-reviewed drug facts (drug_literature.py): the curated file is complete and its tags rest on its records."""
import json

from ucatlas.drug_literature import DATA, records, synonyms


def test_every_record_has_a_paper_and_verbatim_quotes():
    for r in records():
        assert r["pmid"] and r["doi"] and r["journal"] and r["quotes"] and r["genes"]
        assert all(len(q) > 40 for q in r["quotes"])     # whole sentences, not fragments


def test_tags_name_only_drugs_with_a_record():
    d = json.loads(DATA.read_text(encoding="utf-8"))
    drugs = {r["drug"].lower() for r in records()} | {n.lower() for s in synonyms() for n in s["names"]}
    for gene, tag in d["tags"].items():
        named = [w.strip(" ,;()").lower() for w in tag.replace(":", " ").split()]
        for w in named:
            if w.endswith(("mab", "nib", "ixin")) or w.startswith(("sb-", "sch-", "mk-")):
                assert w in drugs, (gene, w)
