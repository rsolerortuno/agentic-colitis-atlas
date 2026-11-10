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

## (e) M14 profiles (decided before measuring; same gold set; corrected in the review of 2026-11-10)

| Score | AUROC | 95 % CI | Recall@100 |
|---|---|---|---|
| Target profile (expression up, tissue change, genetics, predicted tractability) | 0.84 | 0.80–0.89 | 14 % |
| Target profile without tractability | 0.71 | 0.65–0.79 | 10 % |
| Predicted tractability alone | 0.83 | 0.78–0.88 | 2 % |
| Target profile, tractability without ligand history | 0.81 | 0.77–0.87 | 12 % |
| Predicted tractability without ligand history, alone | 0.75 | 0.69–0.82 | 3 % |
| Biomarker profile | 0.76 | 0.69–0.82 | 7 % |

**Correction.** The first version of this table took predicted tractability only for the ~1,100 genes enriched
through the Open Targets API, which are mostly genes already associated with UC; every other gene counted as not
tractable. That favoured clinical targets and made "tractability alone" look as good as the profile. With the Open
Targets bulk table for every gene (same buckets as the kidney atlas), the target profile matches tractability on
AUROC (0.84 vs 0.83) but recovers far more clinical targets in its top 100 (14 % vs 2 %), and without the
medicinal-chemistry ligand bucket it beats tractability clearly (0.81 vs 0.75). As in the kidney atlas, the gain is
in early enrichment. Convergence v1 keeps its original tractability stream (enriched genes only, a known bias that
inflates its 0.86); it is not changed because the published narratives were written from it.

Disease activity: many genes track the Mayo endoscopic score in MSCCR (one value per patient, n = 375; for example
S100A8 rho 0.49, LCN2 0.33) but none reach significance in GSE73661, where baseline scores are mostly 2–3 and n = 67.
Not done: Mendelian randomisation with plasma pQTL (summary statistics need a separate download and a locus-matching
step; left as a limit). The temporal subset (4 genes) stays uninformative.

## (f) M17: repeats and one agent vs the team

Three independent team runs (Opus orchestrator + six Sonnet specialists) and one single-agent run (one Opus agent
doing the six jobs), same question, data, tools and verifier. Checklist of 12 facts fixed before reading the
single-agent report (`src/ucatlas/agent_benchmark.py`).

| Set-up | Run | Statements | Unsupported | Checklist | Specialist evidence cited | Cost | Minutes | Turns |
|---|---|---|---|---|---|---|---|---|
| Team | full | 120 | 0 | 10/12 | 83 % | $4.73 | 7.8 | 92 |
| Team | repeat 1 | 116 | 0 | 10/12 | 80 % | $4.10 | 6.3 | 55 |
| Team | repeat 2 | 119 | 0 | 10/12 | 77 % | $4.32 | 6.5 | 66 |
| One agent | single | 157 | 0 | 10/12 | 90 % | $6.07 | 16.8 | 167 |

- **Quality ties.** Every run passes verification with no unsupported statement and states the same 10 of 12
  checklist facts. Both set-ups miss the neutrophil OSM → OSMR circuit (they discuss OSMR, not its neutrophil
  ligand) and the published infliximab non-response genes.
- **The team is about 27 % cheaper and 2.4 times faster** ($4.38 mean vs $6.07; 6.9 vs 16.8 minutes). The single agent
  writes more statements and cites more of the evidence it collected.
- **Stability.** Target-list Jaccard between team runs 0.65–0.87; CXCR2 first in all three; the single agent's list
  overlaps the team's by 0.59–0.80.

## (g) Verifier sabotage test (review of 2026-11-10)

Every statement of the published run that passes verification was copied and corrupted one way at a time
(`uv run python -m ucatlas.sabotage`, evals/verifier_sabotage.json): a wrong digit in a number 78/78 caught, the
statement given to another gene 83/83, another gene's evidence 83/83, the direction word flipped 16/55 overall and
16/25 when the statement cites significant bulk expression evidence for its gene. The direction check is narrow by
design (it reads only significant per-cohort and meta-analysis changes); a flipped direction in a sentence about a
paper or a cell type is not detected. Two fixes came from this test: the direction check now also reads the
meta-analysis records (before, 1 of 55 flips was even eligible), and ordering words ("moved down") and "up to" are no
longer read as expression directions.
