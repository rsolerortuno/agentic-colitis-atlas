"""AlphaGenome layer: validation holds, thresholds come from it, risk-allele signs and graph records are consistent."""
import json

import pandas as pd
import pytest

from ucatlas.paths import PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "alphagenome_uc.parquet").exists(), reason="run ucatlas.alphagenome_uc")


def summary():
    return json.loads((PROCESSED / "alphagenome_summary.json").read_text(encoding="utf-8"))


def test_validation_supports_use():
    s = summary()
    assert s["validation_atlas"]["pairs"] >= 200
    assert s["validation_atlas"]["agreement_top_25pct"] >= 0.85
    m, mm = s["validation_ag2"]["sc_expression_matched"], s["validation_ag2"]["sc_expression_mismatched"]
    assert m["agreement_top_quarter"] >= 0.85 and m["spearman"] > mm["spearman"]


def test_confident_means_above_validation_threshold():
    s = summary()
    d = pd.read_parquet(PROCESSED / "alphagenome_uc.parquet")
    thr = {"Atlas": s["validation_atlas"]["threshold_top_quarter"],
           "AlphaGenome 2": s["validation_ag2"]["sc_expression_matched"]["threshold_top_quarter"]}
    for src, x in d.groupby("source"):
        assert ((x.effect.abs() >= thr[src]) == x.confident).all()


def test_loci_and_credible_sets():
    g = pd.read_parquet(PROCESSED / "alphagenome_uc_loci.parquet")
    assert g.locus.nunique() == summary()["loci"]
    cs = g[g.in_cs].groupby("locus").pip.sum()
    assert (cs >= 0.95 - 1e-6).all()
    assert set(g.risk_sign.unique()) <= {-1, 1}


def test_il10_direction_matches_protective_mr():
    d = pd.read_parquet(PROCESSED / "alphagenome_uc.parquet")
    x = d[(d.gene == "IL10") & d.confident]
    if x.empty:
        pytest.skip("IL10 not confident in this build")
    assert (x.effect < 0).all()   # risk allele lowers IL10; MR: more IL10 protects
