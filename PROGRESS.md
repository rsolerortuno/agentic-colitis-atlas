# PROGRESS

Spec: [OBJETIVO.md](OBJETIVO.md). A milestone is closed only when its check passes; evidence is recorded below.

**Current milestone:** M3 (literature)

| Milestone | Status | Evidence |
|---|---|---|
| M0 Environment and generic engine | ✅ passed | uv + Python 3.12; engine copied from the kidney atlas (GraphStore, evidence ids, verifier with numeric/direction/staleness checks, coverage, provenance, reproducibility, BM25 search) and made disease-agnostic: the direction check now reads disease words (inflamed, lesional, colitis, UC, tumour) against control words (healthy, normal, control, non-inflamed). `tests/test_engine.py` 10 passed on a toy graph, including a check that no oncology-specific term remains in the engine |
| M1 Bulk data check and download | ✅ passed | Five cohorts checked against GEO before download (data/README.md) and downloaded with sha256 manifest. Baseline colon samples: MSCCR GSE193677 UC inflamed 293 / non-inflamed 391 / control 340 (135 UC patients with both); UNIFI GSE206285 UC 550 / healthy 18; PURSUIT GSE92415 87 / 21; GSE16879 24 / 12 (8 infliximab responders); GSE73661 67 / 12. Response = earliest post-treatment label. Caveat recorded: infliximab patients of GSE16879 and GSE73661 probably overlap. SCP259 single-cell files checked (dimensions match, 365,492 cells) with manifest. `tests/test_cohorts.py` green |
| M2 Differential expression and meta-analysis | ✅ passed | Per-cohort DE (patients averaged across biopsies; Welch or paired t-test, BH) and DerSimonian-Laird meta-analysis. MSCCR inflamed vs control: 1,086 genes with padj<0.05 and abs(log2FC)>=1; replicated in the same direction in UNIFI 83 %, PURSUIT 78 %, GSE73661 81 %, GSE16879 52 % (Spearman 0.64-0.77). Meta-analysis: 1,256 genes significant with abs(log2FC)>=1 of 24,768. All 12 positive controls significant in the expected direction in >=4 of 5 cohorts (S100A8, LCN2, REG4, DUOX2, CXCL8, MMP3, IL1B, OSM up; AQP8, SLC26A2, HMGCS2, ABCG2 down). Paired inflamed vs non-inflamed: 13,171 genes padj<0.05. `tests/test_expression.py` green |
| M3 Literature | pending | |
| M4 Knowledge graph | pending | |
| M5 MCP servers | pending | |
| M6 Agents and orchestrator | pending | |
| M7 Web report | pending | |
| M8 Evaluation | pending | |
| M9 Reproducibility | pending | |
| M10 Publication | pending | |
| M11 Single cell and spatial | pending | |
| M12 Cell composition | pending | |
| M13 Discovery of new candidates | pending | |
| M14 Genetics and independent layers | pending | |
| M15 Pathway layer | pending | |
| M16 Treatment response | pending | |
| M17 Repeats and agent benchmark | pending | |
| M18 Full report and clean-up | pending | |

## Blockers / notes
- None yet.
