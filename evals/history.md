# Evaluation history (decisions are written here BEFORE results are measured)

## 2026-11-09, M8 v1 (decided before the first evaluation run)
- Gold set: Open Targets 26.09 clinical precedence, disease text matching "ulcerative colitis", phase 2 or later,
  >= 2 trials or an approval, drugs with <= 5 targets. Temporal subset: first UC trial on or after 2018-01-01.
- Scores: convergence v1 (config.yaml), the same without literature and tractability, and single streams.
- No weight of convergence v1 is changed after seeing these results; any later change gets a new version and its
  own entry here.

## 2026-11-09, M8 v1 results and one grading fix
- Results in results.md / results.json. No weight was changed after seeing them.
- QA q07: the answer "13th" was correct but the regex `\b13\b` did not match; the pattern was widened to
  `\b13(th)?\b` and the saved answers were regraded (10/10). Logged here because it was changed after the run.
