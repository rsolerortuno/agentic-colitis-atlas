"""UC GWAS against the three QTL layers (qtl_overlap.py)."""
import json

import pandas as pd
import pytest

from ucatlas import qtl_overlap as qo
from ucatlas.paths import PROCESSED


def test_venn_parts_add_up_to_the_union():
    sets = {"plasma_protein": {"A", "B", "X"}, "colon_tissue": {"B", "C", "X"}, "gut_cell": {"C", "D", "X"}}
    v = qo._venn(sets)
    parts = [v[k] for k in ("plasma_only", "colon_only", "cell_only", "plasma_and_colon_only", "plasma_and_cell_only",
                            "colon_and_cell_only", "all_three")]
    assert sum(parts) == v["any"] == 5 and v["all_three_genes"] == ["X"]


def test_shared_variant_and_colocalisation_are_counted_apart(tmp_path, monkeypatch):
    """A QTL that shares a credible-set variant with a UC signal is not colocalised unless Open Targets says so, and a
    colocalising credible set that is not an instrument of a layer is ignored."""
    monkeypatch.setattr(qo, "PROCESSED", tmp_path)
    pd.DataFrame({"signal": ["S001", "S001", "S002"], "cs": ["uc1", "uc1", "uc2"],
                  "variant": ["1_100_A_G", "1_200_C_T", "2_500_G_A"], "preprint_only": [False, False, False]}
                 ).to_parquet(tmp_path / "genetic_signals.parquet")
    inst = {"plasma_protein": pd.DataFrame({"studyLocusId": ["p1"], "gene": ["GENEP"], "studyId": ["ukb"],
                                            "variantId": ["1_100_A_G"], "locus": [json.dumps([["1_100_A_G", 0.9, 0.1]])]}),
            "colon_tissue": pd.DataFrame({"studyLocusId": ["c1"], "gene": ["GENEC"], "studyId": ["gtex"],
                                          "variantId": ["2_500_G_A"], "locus": [json.dumps([["2_500_G_A", 0.8, 0.2]])]}),
            "gut_cell": pd.DataFrame(columns=["studyLocusId", "gene", "studyId", "variantId", "locus"])}
    answers = {"uc1": [{"other": "p1", "h4": 0.3, "clpp": 0.001, "method": "x"},
                       {"other": "not_an_instrument", "h4": 0.99, "clpp": 0.5, "method": "x"}],
               "uc2": [{"other": "c1", "h4": 0.95, "clpp": None, "method": "x"}]}
    monkeypatch.setattr(qo, "_query", lambda cs: answers[cs])
    t, n = qo.signal_table(inst)
    assert n == 2 and len(t) == 2
    p = t[t.layer == "plasma_protein"].iloc[0]
    c = t[t.layer == "colon_tissue"].iloc[0]
    assert p.shares_variant and not p.colocalises
    assert c.shares_variant and c.colocalises and c.gene == "GENEC"


def test_summary_matches_the_mendelian_randomisation_tables():
    p = PROCESSED / "qtl_overlap.json"
    if not p.exists() or not (PROCESSED / "mr_colon.parquet").exists():
        pytest.skip("QTL overlap not built")
    s = json.loads(p.read_text(encoding="utf-8"))
    mrp = pd.read_parquet(PROCESSED / "mr_pqtl.parquet")
    mrc = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    expect = {"plasma_protein": mrp, "colon_tissue": mrc[mrc.layer == "colon_tissue"], "gut_cell": mrc[mrc.layer == "gut_cell"]}
    for layer, r in expect.items():
        v = s["layers"][layer]
        assert v["mr_genes_tested"] == v["matched_genes"] == len(r)
        assert v["mr_significant"] == int((r.q < 0.05).sum())
        assert v["matched_instruments"] <= v["qtl_signals"] and v["matched_genes"] <= v["qtl_genes"]
        assert v["uc_signals_colocalised"] <= s["uc_signals"]
    t = pd.read_parquet(PROCESSED / "qtl_overlap_signals.parquet")
    assert s["uc_signals_colocalised_any_layer"] == t[t.colocalises].signal.nunique()
    assert s["overlap_uc_signals_colocalised"]["any"] == s["uc_signals_colocalised_any_layer"]
