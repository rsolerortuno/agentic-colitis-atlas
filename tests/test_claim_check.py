"""Reader check of claims against their cited evidence (claim_check.py): validation and scoring, no model call."""
import json

import pytest

from ucatlas import claim_check as cc


def _setup(tmp_path, monkeypatch, reader_lines):
    monkeypatch.setattr(cc, "BASE", tmp_path)
    run = tmp_path / "r"
    (run / "batches").mkdir(parents=True, exist_ok=True)
    (run / "read").mkdir(exist_ok=True)
    ev = [{"evidence_id": "EV:1", "record": '{"gene": "CXCL1", "auc": 0.71}',
           "abstract": "CXCL1 (AUC = 0.710) and CYP2R1 (AUC = 0.863) were markers."}]
    rows = [{"id": "S001", "statement": "CXCL1 has a lower AUC than CYP2R1.", "gene": "CXCL1", "evidence": ev},
            {"id": "S002", "statement": "CXCL1 has a higher AUC than CYP2R1.", "gene": "CXCL1", "evidence": ev}]
    (run / "batches" / "batch_01.jsonl").write_text("\n".join(json.dumps(r) for r in rows), encoding="utf-8")
    (run / "key.json").write_text(json.dumps({"S001": "real|x", "S002": "flip|x"}), encoding="utf-8")
    (run / "read" / "batch_01.jsonl").write_text("\n".join(json.dumps(r) for r in reader_lines), encoding="utf-8")
    return "r"


def test_quotes_must_be_verbatim(tmp_path, monkeypatch):
    run = _setup(tmp_path, monkeypatch, [
        {"id": "S001", "claims": [{"claim": "lower AUC than CYP2R1", "verdict": "supported",
                                   "quote": "CXCL1 (AUC  = 0.710)", "evidence_id": "EV:1"}]},
        {"id": "S002", "claims": [{"claim": "higher AUC than CYP2R1", "verdict": "contradicted",
                                   "quote": "CXCL1 had the best AUC", "evidence_id": "EV:1"}]}])
    rows = {r["id"]: r for r in cc.validate(run)}
    assert rows["real|x"]["claims"][0]["verdict"] == "supported"        # whitespace ignored
    assert rows["flip|x"]["claims"][0]["verdict"] == "not_in_evidence"  # invented quote
    assert "not verbatim" in rows["flip|x"]["claims"][0]["reason"]


def test_every_statement_needs_a_reader_line(tmp_path, monkeypatch):
    run = _setup(tmp_path, monkeypatch, [{"id": "S001", "claims": []}])
    with pytest.raises(ValueError):
        cc.validate(run)


def test_flip_is_caught_only_when_the_flagged_claim_holds_the_flipped_word(tmp_path, monkeypatch):
    ok = {"id": "S001", "claims": [{"claim": "lower AUC than CYP2R1", "verdict": "supported",
                                    "quote": "CYP2R1 (AUC = 0.863)", "evidence_id": "EV:1"}]}
    run = _setup(tmp_path, monkeypatch, [ok, {"id": "S002", "claims": [
        {"claim": "higher AUC than CYP2R1", "verdict": "contradicted", "quote": "CYP2R1 (AUC = 0.863)",
         "evidence_id": "EV:1"}]}])
    assert cc.summarize(run, write=False)["sabotage"]["caught_by_reader"] == 1
    run = _setup(tmp_path, monkeypatch, [ok, {"id": "S002", "claims": [
        {"claim": "CXCL1 is a marker in one study", "verdict": "not_in_evidence", "quote": ""}]}])
    assert cc.summarize(run, write=False)["sabotage"]["caught_by_reader"] == 0
