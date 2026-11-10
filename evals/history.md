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

## 2026-11-09, M12 composition verdict (changed after the first run)
- First rule: 'intrinsic' if the composition-adjusted change stays significant in the same direction. The first run
  showed CXCR2 keeping significance with only 17 % of its change (2.60 -> 0.44 log2FC), which that rule called
  intrinsic. Added: 'intrinsic' also needs >= 50 % of the change retained; otherwise 'partly_composition'.
  Decided after seeing results, so logged here.

## 2026-11-09, M14 convergence v2 and profiles (decided before measuring)
- Convergence v1 is kept unchanged (evaluated in M8). v2 adds three streams with weight 1 each: tissue change
  (composition-retained share of the change, M12), disease activity (max |Spearman rho| with the Mayo endoscopic
  score in the two cohorts that have it, significant only) and cell attribution replication is not scored.
- Target profile (no literature, no drug history): expression up (only increases count), tissue change, genetics,
  predicted tractability (SM or AB tractable, regardless of clinical drugs). Weights 2 / 1 / 2 / 0.5.
- Biomarker profile: |expression| either direction, replication, disease activity, inflamed vs non-inflamed.
  Weights 2 / 1 / 1 / 1.
- Weight sensitivity: 1,000 random reweightings (each weight multiplied by a log-uniform factor in [1/2, 2]),
  5th-95th percentile rank and probability of staying in the top 10 / 25.
- The profiles are evaluated against the same M8 gold set; no weight is changed after seeing the result.

## 2026-11-10, review: tractability for every gene
- Bug: predicted tractability (target profile and the 'tractability alone' baseline) came only from the ~1,100
  genes enriched through the Open Targets API; all other genes scored 0. Fixed with the Open Targets 26.09 bulk
  tractability table (same strong buckets as the kidney atlas). Target profile AUROC 0.80 -> 0.84, predicted
  tractability alone 0.81 -> 0.83, recall@100 14 % vs 2 %. Weights unchanged. Convergence v1 keeps its biased
  tractability stream (documented) because the published narratives were written from it.
