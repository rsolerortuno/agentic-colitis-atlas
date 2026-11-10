# Ulcerative colitis agentic atlas

An agentic, evidence-verified multi-omics GraphRAG for **ulcerative colitis (UC)**: the same system as the
[kidney cancer atlas](https://github.com/rsolerortuno/agentic-kidney-atlas), applied to an immune-mediated disease.

**Question.** Which therapeutic targets, immune pathways and biomarkers of ulcerative colitis have the strongest
convergent evidence?

Report: `reports/uc_report.html` (self-contained; open it in a browser). Plan: [OBJETIVO.md](OBJETIVO.md).
Progress and the check behind every milestone: [PROGRESS.md](PROGRESS.md). Evaluation: [evals/results.md](evals/results.md).

## What it found

- **The strongest convergent signal is the neutrophil chemokine axis (CXCR2, CXCR1, CXCL8)**, but single-cell data
  with neutrophils and a composition adjustment show that most of its rise in biopsies is neutrophils arriving
  (CXCR2 keeps 17 % of its change after adjustment). It is a target on recruited cells, not a change of the tissue.
- **Tissue changes that survive the adjustment**: CHI3L1 and MMP1 in inflammatory fibroblasts, DUOX2 in enterocytes,
  IDO1 in myeloid cells; the loss of epithelial metabolism genes (HMGCS2, SLC26A2, AQP8) is the cleanest biomarker.
- **OSMR** combines genetic support, a clinical antibody and a neutrophil OSM → OSMR circuit to inflammatory
  fibroblasts in the single-cell data.
- **35 first-in-class candidates** pass composition, novelty, tractability and vital-organ filters (for example
  XKR9 and SDR16C5 in epithelium, CPXM1 and CHRDL2 in stroma, LAIR2 and LAX1 in T cells): hypotheses to test.
- **Treatment**: TNF and NF-κB pathway activity is higher before treatment in non-responders to infliximab,
  ustekinumab and vedolizumab, and also in placebo non-responders, so it marks severity more than a drug mechanism.

## How it works

1. **Measure (deterministic Python).** Five public colon-biopsy cohorts (GEO): inflamed UC vs control, non-inflamed
   UC vs control and inflamed vs non-inflamed mucosa of the same patient (Mount Sinai MSCCR, RNA-seq), and UC vs
   healthy in UNIFI, PURSUIT and two Leuven cohorts; random-effects meta-analysis. Three single-cell datasets (one
   with neutrophils), a cell-composition adjustment of the biopsies, Xenium spatial sections, pathway and
   transcription-factor activity, ligand-receptor communication, correlation with the Mayo endoscopic score,
   response to four biologics with placebo arms, and a discovery screen.
2. **Store with receipts.** Every result, literature relation (PubTator3, verbatim-validated) and external fact (Open
   Targets, STRING, UniProt, Reactome, HPA) enters a knowledge graph through an evidence record.
3. **Ask the agents.** Six specialists (Claude Sonnet 5.5: expression, literature, external, cells, pathways,
   treatment and discovery) query the graph through five MCP servers; an orchestrator (Claude Opus 5.5) integrates
   them into a narrative.
4. **Check every sentence.** A verifier drops statements whose evidence does not exist, is not about the gene, does
   not contain the quoted number, contradicts the stated direction or changed since the agents read it; a coverage
   step checks that every key finding reaches the report.

## How it was checked

| Check | Result |
|---|---|
| Verification of the published run | 120 statements, 0 unsupported, 0 stale; 27 of 27 key findings cited |
| Independent repeats | 116 and 119 statements, 0 unsupported; target-list Jaccard 0.65–0.87 |
| Recovery of 59 UC clinical targets | Target profile AUROC 0.84 vs predicted tractability alone 0.83, but 14 % vs 2 % of them in its top 100; convergence 0.86 (partly circular); data only 0.74 |
| Known-answer questions (`claude -p` + MCP) | 10/10 |
| Closed-book LLM without tools | Lists approved UC targets (60–70 % clinical) with nothing to check; the pipeline lists data-supported candidates (100 % in the data-only top 500) with verified evidence |
| Verifier sabotage test | Wrong numbers 78/78, wrong gene 83/83, wrong evidence 83/83 caught; flipped directions 16/25 where bulk evidence is cited |
| One agent vs the team | Same quality (10/12 checklist facts, 0 unsupported); the team is about 27 % cheaper and 2.4x faster |
| Deterministic rebuild | `uv run ucatlas run-all --run 2026-11-09_full` |

## Reproduce

```
uv sync
uv run ucatlas run-all --run 2026-11-09_full   # downloads (sha256), analyses, graph, verification, report
powershell -File scripts/run_informe.ps1      # a new agent run (Claude Code, subscription login)
```

The single-cell files of Smillie et al. 2019 (SCP259) and Gubatan et al. 2026 (SCP3755) need a Broad Single Cell
Portal login and are placed by hand in `data/raw/scp259/` and `data/raw/scp3755/` (sha256 in their manifests).

## Limits

Public data only. No tissue proteomics at the scale of the kidney atlas, and Mendelian randomisation with plasma pQTL
was not done. The composition adjustment can over-correct (lineage scores also rise when resident cells activate).
Two of four Xenium 5K sections had too few typed cells and the panel lacks several candidates. The temporal
evaluation has only 4 targets. Microarray fold changes are compressed. Nothing here validates a target
experimentally.
