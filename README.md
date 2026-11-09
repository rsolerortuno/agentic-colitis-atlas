# Ulcerative colitis agentic atlas

An agentic, evidence-verified multi-omics GraphRAG for **ulcerative colitis (UC)**: the same system as the
[kidney cancer atlas](https://github.com/rsolerortuno/agentic-kidney-atlas), applied to an immune-mediated disease.

**Question.** Which therapeutic targets, immune pathways and biomarkers of ulcerative colitis have the strongest
convergent evidence?

## How it works

1. **Measure (deterministic Python).** Five public colon-biopsy cohorts (GEO) are analysed one by one: inflamed UC vs
   control, non-inflamed UC vs control, and inflamed vs non-inflamed mucosa of the same patient (Mount Sinai MSCCR,
   RNA-seq), and UC vs healthy in four trial and cohort datasets (UNIFI, PURSUIT, two Leuven cohorts). Changes are
   combined by random-effects meta-analysis.
2. **Store with receipts.** Every result, literature relation (PubTator3, verbatim-validated) and external fact (Open
   Targets, STRING, UniProt, Reactome) enters a knowledge graph through an evidence record.
3. **Ask the agents.** Specialists (Claude Sonnet 5.5) query the graph through MCP tools and write findings with
   evidence ids; an orchestrator (Claude Opus 5.5) integrates them into a narrative.
4. **Check every sentence.** A verifier drops statements whose evidence does not exist, is not about the gene, does not
   contain the quoted number, contradicts the stated direction or changed since the agents read it. A coverage step
   checks that the specialists' key findings reach the report.

## Status

See [PROGRESS.md](PROGRESS.md) for the milestones and their checks, [OBJETIVO.md](OBJETIVO.md) for the plan,
[data/README.md](data/README.md) for the cohorts and [evals/results.md](evals/results.md) for the evaluation.

## Reproduce

```
uv sync
uv run ucatlas run-all          # downloads (sha256), analyses, graph, convergence, verification, report
powershell -File scripts/run_informe.ps1   # a new agent run (Claude Code, subscription login)
```

The single-cell files of Smillie et al. 2019 (SCP259) need a Broad Single Cell Portal login and are placed by hand in
`data/raw/scp259/` (sha256 in its manifest).

## Limits

Public data only; bulk biopsy changes mix tissue changes with infiltrating cells until the composition layer;
microarray fold changes are compressed; no tissue proteomics at the scale of the kidney atlas. Nothing here validates a
target experimentally.
