"""M0: the disease-agnostic engine works on a toy graph (no downloaded data needed)."""
import json
import re
from pathlib import Path

import pytest

from ucatlas import coverage, provenance, verify
from ucatlas.graph import GraphStore, evidence_id

ENGINE = Path(__file__).resolve().parents[1] / "src" / "ucatlas"


@pytest.fixture()
def toy(tmp_path, monkeypatch):
    # no HGNC file in a unit test: genes resolve to themselves
    monkeypatch.setattr(verify, "resolve_gene", lambda g: g)
    st = GraphStore(tmp_path / "kg.sqlite")
    with st.batch():
        st.add_node("gene:LCN2", "Gene", "LCN2")
        st.add_node("gene:AQP8", "Gene", "AQP8")
        st.add_node("MONDO:0005101", "Disease", "ulcerative colitis")
        up = st.add_evidence(evidence_id("omics_de", "LCN2", "rna"), "omics_de", "DATASET:TOY", gene="LCN2",
                             layer="rna", log2fc=3.1, padj=1e-9, significant=True)
        dn = st.add_evidence(evidence_id("omics_de", "AQP8", "rna"), "omics_de", "DATASET:TOY", gene="AQP8",
                             layer="rna", log2fc=-4.2, padj=1e-12, significant=True)
        st.add_claim("gene:LCN2", "upregulated_in", "MONDO:0005101", up)
        st.add_claim("gene:AQP8", "downregulated_in", "MONDO:0005101", dn)
    return st, up, dn


def test_evidence_ids_are_deterministic():
    assert evidence_id("omics_de", "LCN2", "rna") == evidence_id("omics_de", "LCN2", "rna")
    assert evidence_id("omics_de", "LCN2", "rna").startswith("EV:omics_de:")


def test_verifier_keeps_linked_and_drops_unlinked_or_missing(toy):
    st, up, dn = toy
    v = verify.Verifier(store=st, semantic_check=False)
    ok, _ = v.check_statement({"gene": "LCN2", "statement": "LCN2 rises.", "evidence_ids": [up]})
    assert ok and ok["evidence_ids"] == [up]
    bad, log = v.check_statement({"gene": "LCN2", "statement": "LCN2 rises.", "evidence_ids": [dn, "EV:x:none"]})
    assert bad is None and log["nonexistent"] == ["EV:x:none"] and log["not_linked_to_gene"] == [dn]


@pytest.mark.parametrize("text,kept", [
    ("LCN2 is higher in inflamed mucosa.", True),
    ("LCN2 is lower in healthy mucosa than in colitis.", True),
    ("LCN2 is lower in inflamed mucosa.", False),
    ("LCN2 is decreased in UC.", False),
])
def test_direction_check_uses_disease_and_control_words(toy, text, kept):
    st, up, _ = toy
    v = verify.Verifier(store=st)
    out, log = v.check_statement({"gene": "LCN2", "statement": text, "evidence_ids": [up]})
    assert (out is not None) == kept, log


def test_non_inflamed_is_the_control_side(toy):
    st, _, dn = toy
    v = verify.Verifier(store=st)
    out, log = v.check_statement({"gene": "AQP8", "statement": "AQP8 is higher in non-inflamed mucosa.",
                                  "evidence_ids": [dn]})
    assert out is not None, log


def test_basis_flags_changed_evidence(toy, tmp_path):
    st, up, _ = toy
    basis = provenance.start_basis(tmp_path, st, [up])
    with st.batch():
        st.add_evidence(up, "omics_de", "DATASET:TOY", gene="LCN2", layer="rna", log2fc=9.9, padj=1e-9,
                        significant=True)
    v = verify.Verifier(store=st, semantic_check=False, basis=basis)
    out, log = v.check_statement({"gene": "LCN2", "statement": "LCN2 rises.", "evidence_ids": [up]})
    assert out is None and log["stale_evidence"] == [up]


def test_coverage_lists_uncited_key_findings(tmp_path):
    (tmp_path / "agent_expression.verified.json").write_text(json.dumps({"findings": [
        {"gene": "LCN2", "key": True, "statement": "x", "evidence_ids": ["EV:omics_de:a"]},
        {"gene": "AQP8", "key": True, "statement": "y", "evidence_ids": ["EV:omics_de:b"]}]}), encoding="utf-8")
    (tmp_path / "narrative.json").write_text(json.dumps({"summary": [{"text": "LCN2", "evidence_ids": ["EV:omics_de:a"]}],
                                                         "sections": {}, "targets": [], "limitations": []}), encoding="utf-8")
    c = coverage.coverage(tmp_path)
    assert c["key_findings"] == 2 and [m["gene"] for m in c["key_findings_not_cited"]] == ["AQP8"]


def test_engine_has_no_oncology_specific_terms():
    bad = re.compile(r"ccrcc|renal|kidney|cptac|tcga|depmap|\bvhl\b", re.I)
    hits = [f"{p.name}:{i}" for p in ENGINE.rglob("*.py") for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1)
            if bad.search(line)]
    assert not hits, hits
