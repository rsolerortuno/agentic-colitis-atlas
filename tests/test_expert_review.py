"""Decision rule, PubMed quote checks and the expert review's bookkeeping (2026-11-27, 2026-11-28). No network."""
import re

from ucatlas import expert_review as X
from ucatlas import decision, pubmed_quotes as P


def _pool():
    return {"first_in_class": [{"gene": g} for g in ("A1", "A2", "A3", "A4", "A5")],
            "known": [{"gene": g} for g in ("K1", "K2", "K3", "K4")],
            "repurposing": [{"gene": g} for g in ("R1", "R2", "R3", "R4")],
            "genetic": [{"gene": g} for g in ("G1", "G2", "G3", "G4")]}


def test_decision_rule_needs_four_per_list_from_the_pool_and_no_clinic_gene(monkeypatch):
    monkeypatch.setattr(decision, "pool", _pool)
    monkeypatch.setattr(decision, "_clinical_genes", lambda: {"K4"})
    ok = [{"gene": g, "list": k, "evidence_ids": ["e"]} for k, rows in _pool().items() for g in [r["gene"] for r in rows][:4]]
    c = decision.check({"targets": ok})
    assert c["problems"] == ["K4: has a drug tried or approved in UC or IBD"]
    bad = ok[:3] + [{"gene": "X9", "list": "first_in_class", "evidence_ids": ["e"]},
                    {"gene": "A1", "list": "known", "evidence_ids": []}] + ok[4:]
    p = decision.check({"targets": bad})["problems"]
    assert "X9: not in the pool of first_in_class" in p
    assert "A1: not in the pool of known" in p and "A1: chosen twice" in p and "A1: cites no evidence" in p
    assert any(x.startswith("known: 5 genes chosen") for x in p)


def test_verdict_reasons_need_a_citation_unless_they_are_judgement(monkeypatch):
    from ucatlas import expert_review as X
    import ucatlas.pubmed_quotes as PQ
    monkeypatch.setattr(PQ, "fetch_many", lambda pmids: {1: _rec()})
    key = {"C01": "targets.CFB", "C02": "narrative.summary[0]"}
    vs = [{"id": "C01", "verdict": "concerns", "summary": "s", "reasons": [
              {"kind": "judgement", "text": "a red safety light"},
              {"kind": "factual", "text": "binds p40", "pmid": 1, "quote": "Drug X binds p40."},
              {"kind": "factual", "text": "uncited"},
              {"kind": "factual", "text": "wrong quote", "pmid": 1, "quote": "Drug X binds p19."}]},
          {"id": "C02", "verdict": "support", "reasons": []},
          {"id": "C01", "verdict": "maybe", "reasons": []}]
    out, dropped = X.check_verdicts(vs, key)
    assert list(out) == ["CFB"] and out["CFB"]["verdict"] == "concerns"
    assert [r["text"] for r in out["CFB"]["reasons"]] == ["a red safety light", "binds p40"]
    assert len(dropped) == 4


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


def test_clinical_list_marks_approved_and_keeps_its_criteria():
    import pytest
    from ucatlas.clinical_list import MIN_COHORTS, MIN_PAPERS, OUT, view
    if not OUT.exists():
        pytest.skip("clinical list not built")
    v = view()
    rows = v["rows"]
    assert v["approved"] == sum(r["approved"] for r in rows) > 0
    for r in rows:
        if not r["approved"]:
            assert r["cohorts_same_direction"] >= MIN_COHORTS and r["uc_papers"] > MIN_PAPERS
            assert r["meta_padj"] is not None and r["meta_padj"] < 0.05
    scores = [r["score"] for r in rows if r["score"] is not None]
    assert scores == sorted(scores, reverse=True)
    assert {"TNF", "IL12B", "JAK2"} <= set(v["approved_genes"])
