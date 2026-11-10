"""Checks for the layers added on 2026-11-10 to close four documented limits: pQTL Mendelian randomisation, wider
spatial coverage, a wider temporal set and the single-cohort rule of the meta-analysis."""
import json

import numpy as np
import pandas as pd
import pytest

from ucatlas.paths import PROCESSED


# ------------------------------------------------------------------ Mendelian randomisation (no data needed)
def test_wald_ratio_and_ivw():
    from ucatlas.mr import estimate
    h = pd.DataFrame({"gene": ["A", "B", "B"], "studyId": ["s1", "s2", "s2"], "lead_variant": ["1_1_A_G", "1_2_C_T", "1_3_C_T"],
                      "beta_exposure": [0.5, 0.4, 0.2], "p_exposure": [1e-20, 1e-30, 1e-10],
                      "beta_outcome": [0.1, 0.08, 0.04], "se_outcome": [0.02, 0.02, 0.02],
                      "substitute_pp": [None, None, None]})
    r = estimate(h).set_index("gene")
    assert r.at["A", "beta"] == pytest.approx(0.2) and r.at["A", "se"] == pytest.approx(0.04)
    assert r.at["A", "method"] == "Wald ratio" and r.at["B", "method"] == "IVW" and r.at["B", "n_instruments"] == 2
    assert r.at["B", "beta"] == pytest.approx(0.2)          # both instruments give a ratio of 0.2


def test_allele_alignment_flips_the_outcome_sign():
    """Open Targets betas are for the ALT allele; a GWAS row that reports REF as the effect allele must be flipped."""
    from ucatlas.mr import _lookup
    gw = {("1", "100"): [{"effect": "A", "other": "G", "beta": 0.3, "se": 0.1, "p": 0.01, "rsid": "rs1"}]}
    assert _lookup(gw, "1_100_G_A")[0] == pytest.approx(0.3)   # alt A is the GWAS effect allele
    assert _lookup(gw, "1_100_A_G")[0] == pytest.approx(-0.3)  # alt G is the GWAS other allele
    assert _lookup(gw, "1_100_C_T") is None


def test_palindromic_variants_are_detected():
    from ucatlas.mr import _palindromic
    assert _palindromic("A", "T") and _palindromic("G", "C") and not _palindromic("A", "G")


@pytest.mark.skipif(not (PROCESSED / "mr_pqtl.parquet").exists(), reason="MR not run")
def test_mr_results_are_sane():
    r = pd.read_parquet(PROCESSED / "mr_pqtl.parquet")
    assert len(r) > 1000 and r["gene"].is_unique
    assert ((r["q"] >= r["p"] - 1e-12)).all()
    assert (r.loc[r["q"] >= 0.05, "colocalises"].isna()).all()   # colocalisation is only fetched for significant proteins


# ------------------------------------------------------------------ meta-analysis: a pooled call needs two cohorts
@pytest.mark.skipif(not (PROCESSED / "meta.parquet").exists(), reason="expression not run")
def test_single_cohort_genes_have_no_pooled_significance():
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    assert m.loc[m["k"] < 2, "meta_padj"].isna().all()
    assert m["k"].max() == 8
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    assert {"GSE109142", "GSE117993", "GSE165512"} <= set(de["cohort"])


# ------------------------------------------------------------------ spatial
@pytest.mark.skipif(not (PROCESSED / "spatial_visium.parquet").exists(), reason="Visium not run")
def test_visium_compartments_match_known_biology():
    v = pd.read_parquet(PROCESSED / "spatial_visium.parquet").set_index("gene")
    assert v.at["DUOX2", "compartment"] == "epithelial"
    assert v.at["LCN2", "compartment"] == "epithelial"
    assert v.at["S100A8", "compartment"] == "neutrophil"
    assert v.at["MMP3", "compartment"] == "fibroblast"
    for g in v.index:
        assert set(json.loads(v.at[g, "lineage_rho"])) >= {v.at[g, "compartment"], v.at[g, "second_compartment"]}


@pytest.mark.skipif(not (PROCESSED / "spatial_sections.parquet").exists(), reason="spatial not run")
def test_enough_dense_xenium_sections():
    s = pd.read_parquet(PROCESSED / "spatial_sections.parquet")
    good = s[s["assigned_fraction"] >= 0.5]
    assert (good["state"] == "inflamed").sum() >= 2 and len(good) >= 4


# ------------------------------------------------------------------ temporal set
@pytest.mark.skipif(not (PROCESSED.parents[1] / "evals" / "results.json").exists(), reason="evaluation not run")
def test_wider_temporal_set_is_large_enough():
    ev = json.loads((PROCESSED.parents[1] / "evals" / "results.json").read_text(encoding="utf-8"))
    assert len(ev["gold_sets"]["temporal_wide_targets"]) >= 15
