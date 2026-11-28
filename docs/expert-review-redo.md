# Redoing the expert review of a run

The expert review (step 4c of `/informe`) can be redone on its own, for example with another model, without
re-running the specialists, the decision pool or the orchestrator's choice. The run of 2026-11-28 was reviewed by
Claude Sonnet because the Opus 5.5 reviewer was stopped twice by the model's safeguards; these are the steps to redo
it.

## What the reviewer reads

`data/expert_review/<run>/input.json`: the question, the decision rule, the four lists, and 20 items (the summary
paragraphs and the 16 decision cards, each with its list, the gene's dossier and the orchestrator's rationale). It is
written by `uv run python -m ucatlas.expert_review prepare <run>` from `runs/<run>/narrative.json`. Its opaque item
ids (C01, C02, ...) map back to the narrative through `key.json`.

The instructions are `.claude/agents/revisor-experto.md` (everything after the front matter). The reviewer may only
search PubMed: with the kg-mcp tools `pubmed_search` and `pubmed_abstract`, or, when those are not loaded, with

```
uv run python -m ucatlas.pubmed_quotes search "<PubMed query>" -n 8
uv run python -m ucatlas.pubmed_quotes abstract <pmid>
```

It writes `data/expert_review/<run>/review.json`.

## Prompt for the reviewer

> Read the file `.claude/agents/revisor-experto.md` and follow its instructions exactly (everything after the front
> matter is your role and task). Search the literature only with `uv run python -m ucatlas.pubmed_quotes search
> "<query>" -n 8` and `uv run python -m ucatlas.pubmed_quotes abstract <pmid>` (set PYTHONIOENCODING=utf-8). Do not
> read any other file of the repository and do not run any other command of the project.
> Input: `data/expert_review/<run>/input.json`. Output: write `data/expert_review/<run>/review.json`.
> When done, reply with one line: the verdicts by kind and the number of objections by kind.

## After the review

1. Keep the previous review for comparison (for 2026-11-28 it is already saved as `review_sonnet.json`,
   `runs/2026-11-28/expert_verdicts_sonnet.json` and `runs/2026-11-28/revisions_expert_sonnet.json`).
2. `uv run python -m ucatlas.expert_review validate <run>`: keeps factual objections and factual verdict reasons only
   with a verbatim abstract sentence of a peer-reviewed, non-retracted paper, registers those quotes as citable
   evidence, and writes `runs/<run>/revisions_expert.json` and `runs/<run>/expert_verdicts.json` (the verdicts shown on
   the decision cards).
3. Apply each kept objection to `runs/<run>/narrative.json` once (a factual one by citing its `evidence_ids` or
   removing the disputed claim; a framing one only by rewording, reordering or removing). For each gene the reviewer is
   `against`, keep it with one sentence answering the reviewer, or replace it by another gene of the same list's pool
   and say why in that list's `decision` paragraph; replaced genes get a second review
   (`expert_review prepare <run> GENE1,GENE2`, the reviewer, `expert_review validate <run> --second-pass`).
4. Re-run the steps after the review: `claim_check prepare-narrative`, the readers and `claim_check revisions-narrative`
   (fix what they list), `scripts/assemble.py`, `scripts/verify.py final`, `claim_check prepare`, the readers,
   `python -m ucatlas.sabotage`, `claim_check summarize`, then `scripts/render.py` and `scripts/validate_report.py`.
5. To measure the new reviewer: `expert_review seed <run>`, the reviewer on `seeded/A/input.json` and
   `seeded/B/input.json` (outputs `review.json` in each), then `expert_review score <run>` (writes
   `evals/expert_review.json`).
6. Update the limitation of the narrative that names the model of the review, and `evals/results.md`.
