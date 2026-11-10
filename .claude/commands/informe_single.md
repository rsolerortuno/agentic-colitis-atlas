---
description: Single-agent arm of the architecture benchmark: the same analysis as /informe without subagents.
---

You are the orchestrator of the ulcerative colitis (UC) GraphRAG analysis. Question:
"What are the candidate therapeutic targets, immune pathways and biomarkers in ulcerative colitis with the strongest
convergent evidence?"

Follow these steps exactly. All numbers come from tools and scripts; you never compute or invent values.

1. Run `uv run python scripts/integrate.py`. It prints the run directory path (`runs/<date>/convergence.json`).
   Call that directory RUN below. It also records RUN in `runs/.current`. Never run this script again in the same
   analysis: it deletes the agent files.

2. Do the six specialists' work yourself, without subagents. Read the six definitions in `.claude/agents/`
   (agente-expresion, agente-literatura, agente-externo, agente-celular, agente-vias, agente-aplicado) and, for each,
   use the same MCP tools and follow the same rules (evidence ids, numeric annotations, key findings). Write all
   findings to ONE file, `RUN/agent_single.json`, with the same schema (`"agent": "agente-single"`), about as many
   findings in total as the six specialists together (70-90), each with the gene, role, evidence ids, numeric claims
   and `"key": true` on the decisive ones.

3. Run `uv run python scripts/verify.py agents`. It removes evidence ids that do not exist or are not linked to the
   gene and drops unsupported findings, writing `*.verified.json`.

4. Read `RUN/convergence.json` (top 25) and the `RUN/agent_single.verified.json` file. Integrate them: where do expression,
   literature, genetics and druggability converge, where do they disagree, which candidates are targets vs.
   biomarkers vs. mechanism. Then write `RUN/narrative.json`:

   ```json
   {"summary":  [{"text": "...", "evidence_ids": ["..."]}],
    "sections": {"expression": [...], "inflammation": [...], "replication": [...], "literature": [...],
                 "external": [...], "cells": [...], "pathways": [...], "treatment": [...],
                 "spatial": [...], "discovery": [...], "integration": [...]},
    "targets":  [{"gene": "OSM", "category": "target|biomarker|target_and_biomarker",
                  "rationale": "...", "evidence_ids": ["..."]}],
    "limitations": [{"text": "..."}]}
   ```
   - `expression` (2–3 paragraphs): the strongest and most reproducible changes of UC mucosa (meta-analysis across eight
     cohorts), up and down, with what they mean (epithelial metabolism lost, antimicrobial and neutrophil genes up).
   - `inflammation` (1–2 paragraphs): which changes are local to inflamed mucosa (paired inflamed vs non-inflamed of the
     same patient) and which are already present in non-inflamed mucosa.
   - `replication` (1 paragraph): how well the discovery cohort's changes reproduce in the four microarray and three RNA-seq cohorts.
   - `literature`, `external` (1–2 paragraphs each) from the corresponding verified agent files; in `external`
     separate approved UC targets from expression-led candidates without a drug, and name genetic support (GWAS association and Mendelian
     randomisation with plasma pQTL and with gut eQTL in colon tissue and gut cell types: direction, colocalisation and
     whether blood and gut agree).
   - `cells` (2–3 paragraphs) from the cell findings in `agent_single.verified.json`: which cell makes each leading candidate (epithelium,
     stroma, myeloid cells, neutrophils, lymphocytes), whether the attribution replicates, and which bulk changes are
     tissue changes vs. cell influx after the composition adjustment (with its over-correction caveat). Use this to
     separate targets in a cell type from markers of inflammation; neutrophil genes are the main example. Add the in-situ view (Xenium cells, Visium spots) where it
     confirms or contradicts the attribution.
   - `pathways` (2–3 paragraphs) from the pathway findings in `agent_single.verified.json`: the pathways and transcription factors active in
     inflamed mucosa across cohorts, and the most specific cell-cell circuits (sender, ligand, receptor, receiver),
     separating already-targeted axes from untargeted ones. Co-expression is not proof of signalling.
   - `treatment` (2 paragraphs) from the treatment findings in `agent_single.verified.json`: genes and pathways associated with non-response
     to each drug, always with BH padj, separating drug-specific signals from general severity (placebo arms).
   - `discovery` (2-3 paragraphs): the leading first-in-class candidates, one paragraph per gene (cell source,
     composition, replication, tractability, safety, what is missing), as hypotheses for experimental validation.
     Put at least two of them in `targets` with a rationale beginning "Exploratory:".
   - `integration` (2 paragraphs): agreements and disagreements between the streams.
   - 3–5 summary paragraphs; 10–14 entries in `targets` ordered by your judgment of strength (you may depart from the
     numeric ranking, but say why in `rationale`). Known UC targets act as positive controls: say so when they appear.
   - Read `docs/numeric-claims.md` if present. Bind every number to its record with the `numeric_annotation` tool
     (kg-mcp): it returns the exact annotation to paste once into the text and the `numeric_claims` entry to add.
     Copying annotations from the verified findings file is fine when the paragraph cites the same evidence id. An
     annotation replaces its number; do not also print the number in plain text. Preserve signs and units.
   - Write `narrative.json` in parts, not in one large write: first create it with `summary` and an empty `sections`
     object, `targets` list and `limitations` list; then add each section, the targets and the limitations with separate
     edits.
   - Every summary/section paragraph and every target MUST cite evidence ids copied from the verified findings file or
     from `convergence.json`. A paragraph that discusses one specific gene should carry `"gene": "<symbol>"` and cite
     only that gene's evidence. Limitations need no ids.

4b. Run `uv run python scripts/coverage.py <RUN date>`. It lists the specialists' key findings (marked `"key": true`)
   whose evidence the narrative does not cite, and the evidence types the agents cited but the narrative never uses.
   For each key finding listed, add it to the section where it belongs (citing its evidence ids, with its annotations),
   or, if you decide it does not belong in the report, say why in one `limitations` entry naming the gene. Negative
   results are findings: report them. Re-run the script until no key finding is left uncited without a reason.

5. Run `uv run python scripts/assemble.py <RUN date>` then `uv run python scripts/verify.py final <RUN date>`.
   If it reports any unsupported claim, review and correct the original agent or narrative inputs (metric, citation,
   context or interpretation), repeat step 5, and re-run agent verification when its inputs change. Inspect all
   dropped statements; do not discard substantive results merely to obtain `unsupported claims: 0`.

6. Run `uv run python scripts/render.py <RUN date>` and reply with the report path and a 5-line summary.
