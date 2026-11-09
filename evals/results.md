# Evaluation (M8)

Rebuild: `uv run python -m ucatlas.evaluate` (results.json) and `uv run python -m ucatlas.qa` (qa_answers.json).
Decisions were fixed before measuring ([history.md](history.md)).

## (a) Recovery of ulcerative colitis clinical targets

Gold set: 59 genes targeted by drugs tested in UC in phase 2 or later (≥ 2 trials or an approval; Open Targets 26.09
clinical precedence). Universe: 24,788 genes. AUROC with a bootstrap 95 % interval over the gold genes.

| Score | AUROC | 95 % CI | Recall@100 | Median rank |
|---|---|---|---|---|
| Convergence v1 (includes literature and drug history) | 0.86 | 0.82–0.91 | 17 % | 1,052 |
| Convergence without literature and tractability (data only) | 0.74 | 0.67–0.82 | 10 % | 2,348 |
| Expression alone | 0.72 | 0.66–0.79 | 7 % | 4,003 |
| Genetics alone | 0.63 | 0.55–0.71 | 7 % | 7,967 |
| Tractability alone | 0.81 | 0.73–0.87 | 61 % | 78 |

**Reading.**
- The data-only score recovers clinical targets well above chance (AUROC 0.74), but the convergence score is
  partly circular: its tractability stream scores 1 for any gene with a clinical drug, and UC clinical targets by
  definition have one. That is why tractability alone reaches 61 % recall@100.
- UC clinical targets are mostly cytokines, receptors and integrins (TNF, IL-23, α4β7, JAK, S1P) whose bulk mRNA
  changes little in biopsies. Expression-led rankings therefore put them low, as the kidney atlas found for
  mutational drivers. The target and biomarker profiles of M14 (without literature and drug history, and with
  genetics weighted) are where this should improve.

## (b) Temporal subset

Only 4 UC clinical targets had their first UC trial on or after 2018-01-01, so the intervals are very wide
(data-only score AUROC 0.72, 0.53–0.88; recall@100 0). This subset is too small to support a conclusion; it is
reported for completeness.

## (c) Known-answer questions (`claude -p`, Sonnet 5.5, project MCP servers)

**10/10 correct.** One grading regex was widened after the run: q07 expected "13" and the answer said "13th"
(correct, but the pattern `\b13\b` did not match). Logged in [history.md](history.md).

## (d) History of decisions

See [history.md](history.md).

## (e) M14 profiles (decided before measuring; same gold set)

| Score | AUROC | 95 % CI | Recall@100 |
|---|---|---|---|
| Target profile (expression up, tissue change, genetics, predicted tractability) | 0.80 | 0.75–0.86 | 12 % |
| Target profile without tractability | 0.71 | 0.65–0.79 | 10 % |
| Predicted tractability alone | 0.81 | 0.74–0.88 | 14 % |
| Biomarker profile | 0.76 | 0.69–0.82 | 7 % |

**Reading.** Unlike the kidney atlas, the target profile does **not** beat predicted tractability alone on UC
clinical targets: approved UC targets are druggable cytokines, receptors and integrins whose mRNA changes little in
biopsies. The data-only signal (0.71 without tractability) is real but modest. The temporal subset (4 genes) stays
uninformative. Disease activity: many genes track the Mayo endoscopic score in MSCCR (one value per patient, n = 375;
for example S100A8 rho 0.49, LCN2 0.33) but none reach significance in GSE73661, where baseline scores are mostly 2–3
and n = 67. Not done: Mendelian randomisation with plasma pQTL (summary statistics need a separate download and a
locus-matching step; left as a limit).
