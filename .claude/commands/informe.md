---
description: Orchestrate the ulcerative colitis multi-agent analysis (specialists + deterministic integration + verification) and render the report.
---

You are the orchestrator of the ulcerative colitis (UC) GraphRAG analysis. Question:
"What are the candidate therapeutic targets and biomarkers in ulcerative colitis, and how strong is the evidence for
each?"

The atlas answers it through lists, each answering a sub-question (first-in-class, well known with no drug,
repurposing, genetics first, in clinical validation), and ends in your decision, which answers the sub-question
"Which new therapeutic targets for ulcerative colitis, not yet treated by any drug in UC, have the strongest evidence
to be tested first?"

Follow these steps exactly. All numbers come from tools and scripts; you never compute or invent values.

1. Run `uv run python scripts/integrate.py`. It prints the run directory path (`runs/<date>/convergence.json`).
   Call that directory RUN below. It also records RUN in `runs/.current`. Never run this script again in the same
   analysis: it deletes the agent files. Then run `uv run python -m ucatlas.decision <RUN date>`: it writes
   `RUN/decision_pool.json`, the candidates of the decision (group A of each list of targets with no therapy in UC)
   with a dossier per gene.

2. Launch the specialists **in parallel** (one message, one Agent/Task call each), each with this prompt (replace RUN
   with the absolute path):
   - `agente-expresion`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_expression.json."
   - `agente-literatura`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_literature.json."
   - `agente-externo`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_external.json."
   - `agente-celular`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_cells.json."
   - `agente-vias`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_pathways.json."
   - `agente-aplicado`: "Question: <question>. Decision pool: RUN/decision_pool.json. Write your output to RUN/agent_applied.json."

3. Run `uv run python scripts/verify.py agents`. It removes evidence ids that do not exist or are not linked to the
   gene and drops unsupported findings, writing `*.verified.json`.

3b. Reading check of the agents (revision loop). Run `uv run python -m ucatlas.claim_check prepare-agents <RUN date>`.
   It writes batches of the verified findings with their cited records to `data/claim_check/<RUN date>/agents/batches/`
   and prints them. Launch `lector-afirmaciones` agents in parallel (one per two or three batches; prompt: "Process
   <input> -> <output>" with the outputs in `data/claim_check/<RUN date>/agents/read/`, same file names). Then run
   `uv run python -m ucatlas.claim_check revisions <RUN date>`: it writes `RUN/revisions_<agent>.json` with each
   finding's claims that are contradicted or not in the cited evidence. For each agent with revisions, resume that
   specialist (SendMessage to it, or a new Agent call of the same type if it cannot be resumed) with: "The reading
   check found claims of your findings that the cited records do not support: RUN/revisions_<agent>.json. For each,
   cite the record that supports it (see 'Every claim needs its own record' in your instructions), correct it, or
   remove it. Rewrite RUN/agent_<name>.json." Then re-run step 3. Do this loop once.

4. Read `RUN/convergence.json` (top 25) and the `RUN/agent_*.verified.json` files. Integrate them: where do expression,
   literature, genetics and druggability converge, where do they disagree, which candidates are targets vs.
   biomarkers vs. mechanism. Then write `RUN/narrative.json`:

   ```json
   {"summary":  [{"text": "...", "evidence_ids": ["..."]}],
    "sections": {"expression": [...], "inflammation": [...], "replication": [...], "literature": [...],
                 "external": [...], "cells": [...], "pathways": [...], "treatment": [...],
                 "spatial": [...], "discovery": [...], "integration": [...], "decision": [...]},
    "targets":  [{"gene": "CFB", "list": "first_in_class|known|repurposing|genetic",
                  "category": "target|target_and_biomarker", "rationale": "...", "evidence_ids": ["..."]}],
    "limitations": [{"text": "..."}]}
   ```
   - `expression` (2–3 paragraphs): the strongest and most reproducible changes of UC mucosa (meta-analysis across eight
     cohorts), up and down, with what they mean (epithelial metabolism lost, antimicrobial and neutrophil genes up).
   - `inflammation` (1–2 paragraphs): which changes are local to inflamed mucosa (paired inflamed vs non-inflamed of the
     same patient) and which are already present in non-inflamed mucosa.
   - `replication` (1 paragraph): how well the discovery cohort's changes reproduce in the four microarray and three RNA-seq cohorts,
     and at protein level in the independent proteomics study when the expression agent reports it.
   - `literature`, `external` (1–2 paragraphs each) from the corresponding verified agent files; in `external`
     separate approved UC targets from expression-led candidates without a drug, and name genetic support (GWAS association and Mendelian
     randomisation with plasma pQTL and with gut eQTL in colon tissue and gut cell types: direction, colocalisation and
     whether blood and gut agree, and whether results replicate in FinnGen), and mouse models (MGI colitis
     phenotypes) where the external agent reports them.
   - `cells` (2–3 paragraphs) from `agent_cells.verified.json`: which cell makes each leading candidate (epithelium,
     stroma, myeloid cells, neutrophils, lymphocytes), whether the attribution replicates across the four single-cell
     atlases (Smillie 2019, Garrido-Trigo 2023, Gubatan 2026, GSE282122: replicated = the lineage is dominant in at
     least two of the three others; name the atlases), and which bulk changes are
     tissue changes vs. cell influx after the composition adjustment (with its over-correction caveat). Use this to
     separate targets in a cell type from markers of inflammation; neutrophil genes are the main example. Add the in-situ view (Xenium cells, Visium spots) where it
     confirms or contradicts the attribution.
   - `spatial` (1–2 paragraphs) from the cells agent's in-situ findings: which cell group explains each leading
     candidate in the deconvolved Visium spots and which cells surround its high spots; say where tissue and single
     cells agree and where they do not.
   - `pathways` (2–3 paragraphs) from `agent_pathways.verified.json`: the pathways and transcription factors active in
     inflamed mucosa across cohorts, and the most specific cell-cell circuits (sender, ligand, receptor, receiver),
     separating already-targeted axes from untargeted ones. Co-expression is not proof of signalling.
   - `treatment` (2 paragraphs) from `agent_applied.verified.json`: genes and pathways associated with non-response
     to each drug, always with BH padj, separating drug-specific signals from general severity (placebo arms).
   - `discovery` (2-3 paragraphs): how the first-in-class screen finds its candidates and what its leading genes
     have in common (cell source, composition, tractability, safety, what is missing).
   - `integration` (2 paragraphs): agreements and disagreements between the streams.
   - `decision` (one paragraph per list, four paragraphs, each with `"list": "<list key>"`): which genes of the
     list's pool you did not choose and why, one clause per gene, citing that gene's records.
   - `targets` is the decision: exactly four genes from each list of `RUN/decision_pool.json` (16 entries, with
     `list` set to the list's key), in your order of priority inside each list. Read each pool gene's dossier and
     the specialists' findings on it. The `rationale` of each gene answers, in this order and citing records:
     direction (inhibit, or activate/restore: for a gene to activate or restore, say how a drug could do it, or that
     no route is known); the cell that makes it and whether its change is a tissue change or cells arriving; the drug
     route (antibody, small molecule, the existing inhibitor and its stage); safety (what blocking it could harm);
     and the first experiment. Set `category` to `target_and_biomarker` when a cited record supports a biomarker role
     (the dossier's `biomarker_pmids`), else `target`. No gene may have a drug tried or approved in UC or IBD.
     `verify.py final` fails a decision that breaks the rule (decision.py).
   - 3–5 summary paragraphs, around the decision: what the atlas recommends testing first, from which lists and why,
     and what the decision cannot say (these are hypotheses from public data).
   - When a fact you need is in no record (what a drug binds, how a trial ended), you may search PubMed
     (`uv run python -m ucatlas.pubmed_quotes search "<query>"`, `... abstract <pmid>`, or the kg-mcp tools
     `pubmed_search` and `pubmed_abstract`) and register the sentence you will rely on:
     `uv run python -m ucatlas.pubmed_quotes add --run <RUN date> --by orchestrator --pmid <pmid> --genes <G1,G2>
     [--drug <name>] --quote "<sentence copied verbatim from the abstract>"`. It refuses quotes that are not verbatim and
     papers that are not peer-reviewed journal articles or are retracted, and prints the evidence ids to cite. Never
     cite a paper you did not register this way or find in a record.
   - Read `docs/numeric-claims.md` if present. Bind every number to its record with the `numeric_annotation` tool
     (kg-mcp): it returns the exact annotation to paste once into the text and the `numeric_claims` entry to add.
     Copying annotations from the verified agent files is fine when the paragraph cites the same evidence id. An
     annotation replaces its number; do not also print the number in plain text. Preserve signs and units.
   - Write `narrative.json` in parts, not in one large write: first create it with `summary` and an empty `sections`
     object, `targets` list and `limitations` list; then add each section, the targets and the limitations with separate
     edits.
   - Every claim of a paragraph needs a cited record that says it, as for the specialists (drug names and stages,
     comparisons with other genes through `rank_context`, protein function, results of other analyses, papers): see
     "Every claim needs its own record" in `.claude/agents/agente-externo.md`. A drug named without a cited record that
     contains the name fails verification. What no record supports is left out or written as interpretation.
   - Every summary/section paragraph and every target MUST cite evidence ids copied from the verified agent files or
     from `convergence.json`. A paragraph that discusses one specific gene should carry `"gene": "<symbol>"` and cite
     only that gene's evidence. Limitations need no ids.

4b. Run `uv run python scripts/coverage.py <RUN date>`. It lists the specialists' key findings (marked `"key": true`)
   whose evidence the narrative does not cite, and the evidence types the agents cited but the narrative never uses.
   For each key finding listed, add it to the section where it belongs (citing its evidence ids, with its annotations),
   or, if you decide it does not belong in the report, say why in one `limitations` entry naming the gene. Negative
   results are findings: report them. Re-run the script until no key finding is left uncited without a reason.

4c. Expert review of the summary and the decision. Run `uv run python -m ucatlas.expert_review prepare <RUN date>`, then launch one
   `revisor-experto` agent with: "Process <input> -> <output>" (the paths it printed). Then run
   `uv run python -m ucatlas.expert_review validate <RUN date>`: it keeps factual objections and factual verdict
   reasons only when their PubMed quote is verbatim and the paper peer-reviewed (and registers the quote as citable
   evidence), keeps framing objections and judgement reasons, and writes `RUN/revisions_expert.json` and
   `RUN/expert_verdicts.json` (a verdict per decision card, shown on the card). Apply each kept objection to
   `narrative.json` once: a factual one by citing its `evidence_ids` or removing the disputed claim, a framing one
   only by rewording, reordering or removing what is already written. Discarded objections are not acted on.
   For each gene the reviewer is `against`, decide: keep it, adding to its `rationale` one sentence that answers the
   reviewer's reason (citing the reviewer's `evidence_ids` when you use its facts), or replace it by another gene of
   the same list's pool and say why in that list's `decision` paragraph. If you replace genes, run
   `uv run python -m ucatlas.expert_review prepare <RUN date> GENE1,GENE2` (the new genes), launch the reviewer on the
   printed paths, and run `uv run python -m ucatlas.expert_review validate <RUN date> --second-pass`. Once.

4d. Reading check of the narrative. Run `uv run python -m ucatlas.claim_check prepare-narrative <RUN date>`, launch
   `lector-afirmaciones` agents on its batches as in 3b (outputs in `data/claim_check/<RUN date>/narrative/read/`),
   then `uv run python -m ucatlas.claim_check revisions-narrative <RUN date>`, which writes
   `RUN/revisions_narrative.json`. Fix each listed claim in `narrative.json` (cite the supporting record, correct, or
   remove). Once.

5. Run `uv run python scripts/assemble.py <RUN date>` then `uv run python scripts/verify.py final <RUN date>`.
   If it reports any unsupported claim, review and correct the original agent or narrative inputs (metric, citation,
   context or interpretation), repeat step 5, and re-run agent verification when its inputs change. Inspect all
   dropped statements; do not discard substantive results merely to obtain `unsupported claims: 0`.

6. Final reading check, for the report: `uv run python -m ucatlas.claim_check prepare <RUN date>`, launch
   `lector-afirmaciones` agents on `data/claim_check/<RUN date>/batches/` (outputs in `.../read/`), then
   `uv run python -m ucatlas.sabotage <RUN date>` and `uv run python -m ucatlas.claim_check summarize <RUN date>`. Claims
   still contradicted or not in the cited evidence are shown in the report, not hidden.

7. Run `uv run python scripts/render.py <RUN date>` and reply with the report path and a 5-line summary.
