"""Front-page selection rule, PubMed quote checks and the expert review's bookkeeping (2026-11-27). No network."""
import re

from ucatlas import expert_review as X
from ucatlas import frontpage, pubmed_quotes as P


def test_frontpage_rule_flags_unjustified_genes_and_missing_controls(monkeypatch):
    crit = {"IL12B": ["positive control: drug approved in UC on this gene or its class, genetic association 0.89"],
            "JAK2": ["positive control: drug approved in UC on this gene or its class, genetic association 0.85"],
            "CXCR2": ["convergence rank 1 (first 30)"]}
    monkeypatch.setattr(frontpage, "criteria", lambda f: crit)
    f = {"targets": [{"gene": "CXCR2"}, {"gene": "IL12B"}, {"gene": "FOO"}]}
    c = frontpage.check(f)
    assert c["unjustified"] == ["FOO"]
    assert c["missing_controls"] == ["JAK2"]
    assert c["reasons"]["CXCR2"] == ["convergence rank 1 (first 30)"]


def _rec(**kw):
    return {"pmid": 1, "title": "t", "journal": "j", "year": "2020", "doi": None, "retracted": False,
            "publication_types": ["Journal Article"], "abstract": "Drug X binds p40. It worked in UC.", **kw}


def test_pubmed_quotes_must_be_verbatim_and_peer_reviewed():
    assert P.problems(_rec(), ["Drug X binds  p40."]) == []          # whitespace is normalised
    assert P.problems(_rec(), ["Drug X binds p19."])                  # not in the abstract
    assert P.problems(_rec(publication_types=["Preprint"]), [])        # not peer-reviewed
    assert P.problems(_rec(publication_types=["Journal Article", "Letter"]), [])
    assert P.problems(_rec(retracted=True), [])
    assert P.problems(_rec(abstract=""), [])


def test_objections_without_citation_are_discarded_and_framing_kept():
    objs = [{"id": "C01", "kind": "framing", "claim": "a", "objection": "b"},
            {"id": "C02", "kind": "factual", "claim": "a", "objection": "b"},
            {"id": "C03", "kind": "opinion", "claim": "a"}]
    kept, dropped = X.check_objections(objs)
    assert [o["id"] for o in kept] == ["C01"]
    assert {o["id"] for o in dropped} == {"C02", "C03"}
    assert "without a PMID" in next(o for o in dropped if o["id"] == "C02")["reason"]


def test_a_seeded_error_is_caught_only_when_the_objection_points_at_it():
    text = "Positive control. IL12B is the p19 subunit of interleukin-12. It is up."
    s = {"text": text, "inserted": "IL12B is the p19 subunit"}
    assert X._hits({"claim": "IL12B is the p19 subunit of interleukin-12"}, s)
    assert X._hits({"claim": "p19 subunit"}, s)
    assert not X._hits({"claim": "It is up."}, s)
    assert not X._hits({"claim": ""}, s)


def test_every_seed_changes_the_words_it_targets():
    for name, kind, pat, rep in X.SEEDS:
        assert kind in X.KINDS
        m = re.compile(pat)
        assert m.pattern and rep, name
