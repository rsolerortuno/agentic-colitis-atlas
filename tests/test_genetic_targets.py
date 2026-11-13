"""Genetics-first list: signals are non-redundant, tiers follow the validation, directions and points are consistent."""
import json

import numpy as np
import pandas as pd
import pytest

from ucatlas.paths import PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "genetic_targets.parquet").exists(), reason="run ucatlas.genetic_targets")


def test_swap_and_support_are_consistent():
    from ucatlas.genetic_targets import _support, _swap
    assert _swap("1_100_A_G") == "1_100_G_A"
    e = pd.DataFrame({"signal": "S1", "gene": "G", "context": "colon tissue", "source": "Atlas", "lineage": "tissue",
                      "tier": "confident", "variant": ["a", "b", "c"], "pip": [0.5, 0.3, 0.1], "effect": [0.2, 0.1, -0.4]})
    s = _support(e, "confident")
    assert len(s) == 1 and s.pip_support.iloc[0] == pytest.approx(0.8)      # majority direction only
    assert s.effect.iloc[0] == pytest.approx((0.5 * 0.2 + 0.3 * 0.1) / 0.8) and s.top_variant.iloc[0] == "a"


def test_signals_do_not_share_variants_and_skip_the_mhc():
    s = pd.read_parquet(PROCESSED / "genetic_signals.parquet")
    assert not s.duplicated(["variant"]).any()
    assert not ((s.chrom == "6") & s.pos.between(25_000_000, 34_000_000)).any()
    assert set(s.risk_sign.unique()) <= {-1, 1}
    d = s.drop_duplicates("signal")
    assert (d.preprint_only == d.source.str.contains("preprint")).all()
    assert s[s.ibd_signal].groupby("signal").uc_nominal.any().all()         # IBD signals need a UC association


def test_gene_assignment_rule_holds_on_eqtl_and_in_the_effects():
    from ucatlas.genetic_targets import RATIO, _top_gene, thresholds
    th = thresholds()
    assert th["check"]["atlas"]["right_gene"] >= 0.8 and th["check"]["ag2"]["right_gene"] >= 0.8
    assert th["check"]["atlas"]["right_gene"] > th["check"]["top_gene_alone"]["atlas"] + 0.3    # the rule adds a lot
    e = pd.read_parquet(PROCESSED / "genetic_effects.parquet")
    assert (e.ratio >= RATIO).all() and not e.duplicated(["variant", "context", "signal"]).any()
    assert (e[e.source == "Atlas"].score.abs() >= th["atlas_confident"]).all()
    assert (e[e.source == "AlphaGenome 2"].score.abs() >= th["ag2_confident"]).all()
    assert np.allclose(e.effect, e.score * e.risk_sign)
    d = pd.DataFrame({"variant": "v", "gene": ["TNF", "IL10", "IL6"], "score": [0.5, -0.2, 0.01]})
    top = _top_gene(d, ["variant"])
    assert top.gene.tolist() == ["TNF"] and top.ratio.iloc[0] == pytest.approx(2.5)


def test_list_rows_are_consistent():
    d = pd.read_parquet(PROCESSED / "genetic_targets.parquet")
    assert d.gene.is_unique
    assert ((d.effect > 0) == (d.hypothesis == "inhibit")).all()
    assert (d[d.tier == "confident"].pip_support >= 0.5).all()
    assert d[d.tier == "likely"].pip_support.between(0.1, 0.5, inclusive="left").all()
    assert (d[d.group == "A"].tier == "confident").all() and (d[d.group == "A"].points >= 5).all()
    assert (d[d.group == "C"].tier == "likely").all() and set(d.tier) <= {"confident", "likely"}
    s = json.loads((PROCESSED / "genetic_targets_summary.json").read_text(encoding="utf-8"))
    assert s["genes"] == len(d) and s["signals_with_gene"] <= s["signals"]
