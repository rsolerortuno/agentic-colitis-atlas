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


# ------------------------------------------------------------------ 2026-11-10, second round
@pytest.mark.skipif(not (PROCESSED / "mr_colon.parquet").exists(), reason="layer not built")
def test_gut_mr_layers_and_positive_signals():
    r = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    assert set(r.layer) == {"colon_tissue", "gut_cell"}
    col = r[r.layer == "colon_tissue"].set_index("gene")
    # IRF5 is an established IBD risk gene: genetically higher colon expression raises risk and colocalises
    assert col.at["IRF5", "q"] < 0.05 and col.at["IRF5", "beta"] > 0 and col.at["IRF5", "colocalises"] is True


@pytest.mark.skipif(not (PROCESSED / "spatial_visium_summary.json").exists(), reason="layer not built")
def test_visium_deconvolution_attributes_known_cells():
    import json
    s = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())
    t = {x["gene"]: x for x in s["targets"]}
    assert t["S100A8"]["top_in_situ_expresser"] == "Neutrophil"
    assert t["MMP3"]["top_in_situ_expresser"] == "Inflammatory_fibroblast"
    assert t["LCN2"]["top_in_situ_expresser"] == "Epithelial"
    # candidates are never deconvolution markers
    markers = {g for ms in s["marker_genes"].values() for g in ms}
    assert not markers & {"S100A8", "MMP3", "LCN2", "DUOX2", "CXCR2", "IL11"}
    assert s["mean_share_uc"]["Neutrophil"] > s["mean_share_healthy"]["Neutrophil"]


@pytest.mark.skipif(not (PROCESSED / "structure_summary.json").exists(), reason="layer not built")
def test_structure_controls_have_pockets():
    import json
    st = json.loads((PROCESSED / "structure_summary.json").read_text())
    ctrl = [r for r in st["results"] if r.get("control") and r["modeled"]]
    assert len(ctrl) >= 8 and all(r["pockets"] for r in ctrl)


def test_fingerprint_ignores_float_noise_not_changes():
    from ucatlas.provenance import fingerprint
    r = {"evidence_id": "x", "gene": "G", "retained": 0.7976121951555426}
    assert fingerprint(r) == fingerprint({**r, "retained": 0.7976121951555432})
    assert fingerprint(r) != fingerprint({**r, "retained": 0.7977})


# ------------------------------------------------------------------ 2026-11-10, third round
@pytest.mark.skipif(not (PROCESSED / "proteomics_summary.json").exists(), reason="layer not built")
def test_protein_level_replication():
    s = json.loads((PROCESSED / "proteomics_summary.json").read_text())
    assert s["replication_rate"] > 0.6 and s["spearman_log2fc"] > 0.4
    d = pd.read_parquet(PROCESSED / "proteomics.parquet").set_index("gene")
    assert d.at["S100A8", "log2_ratio"] > 1 and d.at["AQP8", "log2_ratio"] < 0


@pytest.mark.skipif(not (PROCESSED / "mr_pqtl.parquet").exists(), reason="MR not run")
def test_mr_controls_replicate_in_finngen():
    r = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").set_index("gene")
    if "replicates_finngen" not in r:
        pytest.skip("FinnGen not attached")
    for g in ("IL12B", "IL10", "CD274"):
        assert r.at[g, "replicates_finngen"] is True


@pytest.mark.skipif(not (PROCESSED / "mouse_colitis_summary.json").exists(), reason="layer not built")
def test_mouse_models_enrich_clinical_targets():
    s = json.loads((PROCESSED / "mouse_colitis_summary.json").read_text())
    assert s["clinical_targets_rate"] > 5 * s["background_rate"]
    d = pd.read_parquet(PROCESSED / "mouse_colitis.parquet").set_index("gene")
    assert d.at["IL10", "direction"] == "loss gives or worsens colitis"


@pytest.mark.skipif(not (PROCESSED / "cosmx_detection.parquet").exists(), reason="layer not built")
def test_cosmx_detects_epithelial_lcn2():
    d = pd.read_parquet(PROCESSED / "cosmx_detection.parquet")
    r = d[(d.gene == "LCN2") & (d.panel == "6k")].iloc[0]
    assert r.best_lineage == "epithelial" and bool(r.above_background)


@pytest.mark.skipif(not (PROCESSED / "samples.parquet").exists(), reason="samples not built")
def test_new_response_cohorts_have_labels_and_no_controls():
    s = pd.read_parquet(PROCESSED / "samples.parquet")
    for c in ("GSE72819", "GSE23597"):
        d = s[(s.cohort == c) & (s.timepoint == "baseline")]
        assert len(d) > 20 and (d.group == "UC").all() and d.response.notna().sum() > 20


def test_split_site_moves_marked_blocks(tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location("split_site", PROCESSED.parents[1] / "scripts" / "split_site.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    html = ('<html><script data-asset="lib.js">var a=1;</script>'
            '<script id="d" type="application/gzip-base64" data-asset="d.b64">QUJD</script></html>')
    (tmp_path / "r.html").write_text(html, encoding="utf-8")
    mod.split(tmp_path / "r.html", tmp_path / "site")
    out = (tmp_path / "site" / "index.html").read_text(encoding="utf-8")
    assert 'src="assets/lib.js"' in out and 'data-src="assets/d.b64"' in out and "QUJD" not in out
    assert (tmp_path / "site" / "assets" / "d.b64").read_text() == "QUJD"


# ------------------------------------------------------------------ 2026-11-10, fourth round
@pytest.mark.skipif(not (PROCESSED / "de_adjusted_summary.json").exists(), reason="layer not built")
def test_changes_hold_after_covariate_adjustment():
    s = json.loads((PROCESSED / "de_adjusted_summary.json").read_text())
    assert s["spearman_adjusted_vs_unadjusted"] > 0.9 and s["share_holding"] > 0.8


@pytest.mark.skipif(not (PROCESSED / "mr_sensitivity.parquet").exists(), reason="layer not built")
def test_mr_sensitivity_on_il12b():
    d = pd.read_parquet(PROCESSED / "mr_sensitivity.parquet")
    r = d[(d.gene == "IL12B") & (d.layer == "plasma_protein")].iloc[0]
    assert r.steiger_correct_direction and r.wmedian > 0 and r.egger_intercept_p > 0.05


@pytest.mark.skipif(not (PROCESSED / "safety.parquet").exists(), reason="layer not built")
def test_safety_phewas_finds_psoriasis_for_il12b():
    d = pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene")
    assert any("soriasis" in t for t in d.at["IL12B", "phewas_traits"])
    assert d.at["JAK1", "constrained"]


@pytest.mark.skipif(not (PROCESSED / "report_layers.json").exists(), reason="layer not built")
def test_ranking_not_explained_by_expression_level():
    nc = json.loads((PROCESSED / "report_layers.json").read_text()).get("negative_control")
    if not nc:
        pytest.skip("negative control not run")
    for k in ("convergence", "target profile"):
        assert nc[k]["auroc_expression_matched"] > nc[k]["null_expression_matched_p95"] + 0.1


@pytest.mark.skipif(not (PROCESSED / "antitnf_summary.json").exists(), reason="layer not built")
def test_antitnf_atlas_counts():
    s = json.loads((PROCESSED / "antitnf_summary.json").read_text())
    assert s["uc_remission"] == 6 and s["uc_non_remission"] == 14


@pytest.mark.skipif(not (PROCESSED / "validation_shortlist.parquet").exists(), reason="layer not built")
def test_validation_shortlist_covers_first_in_class():
    sl = pd.read_parquet(PROCESSED / "validation_shortlist.parquet")
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    assert set(sl.gene) == set(d.gene[d.status == "first-in-class opportunity"])
    assert list(sl["rank"]) == list(range(1, len(sl) + 1))
    for c in sl.criteria:
        lights = json.loads(c)
        assert set(lights) == {"reproducible", "protein", "tissue", "cell", "genetics", "druggable", "safety", "literature", "mouse"}
        assert all(v[0] in ("green", "amber", "red", "grey") and v[1] for v in lights.values())


@pytest.mark.skipif(not (PROCESSED / "validation_shortlist_calibration.parquet").exists(), reason="layer not built")
def test_validation_shortlist_rules_on_markers_and_controls():
    from ucatlas import validation_shortlist as vs
    sl = pd.read_parquet(PROCESSED / "validation_shortlist.parquet").set_index("gene")
    ccno = json.loads(sl.at["CCNO", "criteria"])
    assert sl.at["CCNO", "regeneration_flag"] and "Ciliary dyskinesia" in ccno["safety"][1]
    assert "no extracellular part" in json.loads(sl.at["OSBPL3", "criteria"])["druggable"][1]
    cal = pd.read_parquet(PROCESSED / "validation_shortlist_calibration.parquet")
    assert all(json.loads(c)["druggable"][0] == "green" for c in cal.criteria)
    sc = vs._SingleCell()
    assert sc.per_cell("DUOX2")[0] > 3 and sc.per_cell("AQP8")[0] < -1 and abs(sc.per_cell("TNF")[0]) < 0.3
    assert sc.crypt["OLFM4"] > vs.CRYPT_SHARE > sc.crypt["EPCAM"]


@pytest.mark.skipif(not (PROCESSED / "shortlist_literature.parquet").exists(), reason="layer not built")
def test_shortlist_literature_review():
    d = pd.read_parquet(PROCESSED / "shortlist_literature.parquet").set_index("gene")
    assert d.at["IKBIP", "light"] == "red" and "31826938" in d.at["IKBIP", "mouse_colitis"]
    assert d.at["RFTN1", "light"] == "green"
    rej = json.loads((PROCESSED / "shortlist_literature_rejected.json").read_text(encoding="utf-8"))
    papers = pd.read_parquet(PROCESSED / "shortlist_literature_papers.parquet")
    assert len(rej) == 0 and len(papers) > 150


@pytest.mark.skipif(not (PROCESSED / "validation_shortlist.parquet").exists(), reason="layer not built")
def test_validation_shortlist_groups_and_experiments():
    from ucatlas import validation_shortlist as vs
    sl = pd.read_parquet(PROCESSED / "validation_shortlist.parquet")
    n = len(vs.SCHEMES)
    assert sl.stability.between(0, n).all()
    assert (sl.group[sl.stability == n] == "A").all() and (sl.group[sl.stability == 0] == "C").all()
    assert (sl.group[(sl.stability > 0) & (sl.stability < n)] == "B").all()
    assert sl.stability.is_monotonic_decreasing                      # table order: groups first
    assert all(vs._points(json.loads(c)) == s for c, s in zip(sl.criteria, sl.score))   # base scheme = stored score
    ex = sl.set_index("gene").experiment
    assert "endothelial" in ex["PCDH17"]                              # a vessel gene is tested in endothelium, not fibroblasts
    assert "regulatory T cells" in ex["LAIR2"] and "organoid" in ex["CCNO"] and "Published lead" in ex["RFTN1"]
    bias = json.loads((PROCESSED / "validation_shortlist_summary.json").read_text())["protein_bias"]
    assert bias["median_rna_percentile_quantified"] > bias["median_rna_percentile_not_quantified"] + 0.2


# ------------------------------------------------------------------ crypt niche and the fourth atlas (2026-11-11)
def test_gse282122_member_names_map_to_sample_titles():
    from ucatlas.sc_gse282122 import _title_of
    titles = {"CID003352-2-reup", "CID003376-1"}
    assert _title_of("filtered_processed_data/CID003352-2/filtered_feature_bc_matrix.h5", titles) == "CID003352-2-reup"
    assert _title_of("filtered_processed_data/CID003376-1/filtered_feature_bc_matrix.h5", titles) == "CID003376-1"
    assert _title_of("filtered_processed_data/CID009999-1/filtered_feature_bc_matrix.h5", titles) is None


@pytest.mark.skipif(not (PROCESSED / "sc_celltype_share.parquet").exists(), reason="single cells not processed")
def test_niche_threshold_separates_markers():
    """The crypt-niche threshold sits halfway between a generic fibroblast gene and the niche markers."""
    from ucatlas.validation_shortlist import NICHE_SHARE, _SingleCell, niche_markers
    m = niche_markers(_SingleCell())
    assert all(v >= NICHE_SHARE for v in m["niche"].values())
    assert all(v < NICHE_SHARE for v in m["generic"].values())
    mid = (m["generic"]["COL1A2"] + np.mean(list(m["niche"].values()))) / 2
    assert abs(mid - NICHE_SHARE) < 0.01


@pytest.mark.skipif(not (PROCESSED / "validation_shortlist.parquet").exists(), reason="shortlist not built")
def test_niche_flag_caps_cell_light():
    d = pd.read_parquet(PROCESSED / "validation_shortlist.parquet")
    for r in d[d.niche_flag].itertuples():
        assert json.loads(r.criteria)["cell"][0] != "green" and "crypt niche" in json.loads(r.criteria)["cell"][1]
    assert not d[d.lineage != "stromal"].niche_flag.any()


@pytest.mark.skipif(not (PROCESSED / "sc4_summary.json").exists(), reason="GSE282122 not processed")
def test_gse282122_typing_on_canonical_markers():
    s = json.loads((PROCESSED / "sc4_summary.json").read_text())
    for lin in ("epithelial", "stromal", "endothelial", "myeloid", "T/NK", "plasma", "mast"):
        chk = s["typing_check"][lin]
        assert chk["present"] and chk["own_top"] >= chk["markers"] - 1, (lin, chk)
    sh = pd.read_parquet(PROCESSED / "sc4_celltype_share.parquet", filters=[("health", "==", "inflamed")])
    top = sh.sort_values("share").drop_duplicates("gene", keep="last").set_index("gene").cluster
    assert top["EPCAM"] == "epithelial" and top["COL1A2"] == "stromal" and top["CD3E"] == "T/NK" and top["MZB1"] == "plasma"


@pytest.mark.skipif(not (PROCESSED / "sc_replication4.parquet").exists(), reason="four-atlas replication not run")
def test_four_atlas_replication_and_cell_light():
    r = pd.read_parquet(PROCESSED / "sc_replication4.parquet")
    agree = sum((r[f"lineage_{k}"] == r.lineage_smillie) for k in ("garrido", "gubatan", "gse282122"))
    assert (agree == r.n_agree).all() and (r.replicated == (r.n_agree >= 2)).all()
    assert r.n_agree.mean() > 1.5   # the attribution mostly replicates
    sl = pd.read_parquet(PROCESSED / "validation_shortlist.parquet")
    for x in sl.itertuples():
        if json.loads(x.criteria)["cell"][0] == "green":
            assert x.atlases_agree >= 2
