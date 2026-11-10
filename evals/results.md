# Evaluation (M8)

Rebuild: `uv run python -m ucatlas.evaluate` (results.json) and `uv run python -m ucatlas.qa` (qa_answers.json).
Decisions were fixed before measuring ([history.md](history.md)). **Numbers updated on 2026-11-10** after three RNA-seq
cohorts were added to the meta-analysis (eight cohorts) and single-cohort genes lost their pooled significance;
the previous values are in the git history.

## (a) Recovery of ulcerative colitis clinical targets

Gold set: 59 genes targeted by drugs tested in UC in phase 2 or later (≥ 2 trials or an approval; Open Targets 26.09
clinical precedence). Universe: 24,788 genes. AUROC with a bootstrap 95 % interval over the gold genes.

| Score | AUROC | 95 % CI | Recall@100 | Median rank |
|---|---|---|---|---|
| Convergence v1 (includes literature and drug history) | 0.88 | 0.84–0.92 | 17 % | 1,119 |
| Convergence without literature and tractability (data only) | 0.78 | 0.72–0.85 | 8 % | 2,742 |
| Expression alone | 0.76 | 0.70–0.83 | 5 % | 3,318 |
| Genetics alone | 0.55 | 0.47–0.64 | 7 % | 12,514 |
| Tractability alone | 0.82 | 0.75–0.89 | 59 % | 75 |

**Reading.**
- The data-only score recovers clinical targets well above chance (AUROC 0.78), but the convergence score is
  partly circular: its tractability stream scores 1 for any gene with a clinical drug, and UC clinical targets by
  definition have one. That is why tractability alone reaches 61 % recall@100.
- UC clinical targets are mostly cytokines, receptors and integrins (TNF, IL-23, α4β7, JAK, S1P) whose bulk mRNA
  changes little in biopsies. Expression-led rankings therefore put them low, as the kidney atlas found for
  mutational drivers. The target and biomarker profiles of M14 (without literature and drug history, and with
  genetics weighted) are where this should improve.

## (b) Temporal subsets

**Strict** (as decided first): only 4 UC clinical targets had their first UC trial on or after 2018-01-01, so the
intervals are very wide (target profile AUROC 0.75; recall@100 0). Too small to conclude; kept for continuity.

**Wider** (added 2026-11-10 after the strict subset proved too small; history.md): targets of a drug with ≤ 5 targets
that reached phase 2+ in UC (one trial is enough), first UC trial of any phase on or after 2016-01-01, complete dates;
earlier UC targets removed from the universe. 20 genes (CD40LG, EPCAM, GIPR, GLP1R, HMGCR, HRH1, IL1A, IL1B, IL1RAP,
IL36R, IL4R, MAP3K8, NR1H4, OSMR, PPARA, PTAFR, SELPLG, SLC5A2, SLC6A2, XDH).

| Score | AUROC | 95 % CI | Recall@100 |
|---|---|---|---|
| Convergence v1 | 0.81 | 0.72–0.90 | 15 % |
| Convergence, data only | 0.75 | 0.63–0.86 | 15 % |
| Target profile | 0.81 | 0.73–0.88 | 15 % |
| Target profile, tractability without ligand history | 0.82 | 0.74–0.89 | 15 % |
| Predicted tractability alone | 0.79 | 0.67–0.88 | 0 % |
| Tractability alone (clinical drug history) | 0.61 | 0.48–0.73 | 20 % |
| Expression alone | 0.71 | 0.57–0.82 | 5 % |

Reading: on targets that entered UC trials recently, the data-driven profiles keep their AUROC and put three of the 20
in their top 100 (MAP3K8, OSMR, IL1B), where predicted tractability alone puts none. Retrospective: the data
are not frozen at 2016, and several of these genes (GLP1R, SLC5A2, HMGCR) come from repurposing trials of metabolic
drugs that no expression signal would predict.

## (c) Known-answer questions (`claude -p`, Sonnet 5.5, project MCP servers)

**10/10 correct.** One grading regex was widened after the run: q07 expected "13" and the answer said "13th"
(correct, but the pattern `\b13\b` did not match). Logged in [history.md](history.md).

## (d) History of decisions

See [history.md](history.md).

## (e) M14 profiles (decided before measuring; same gold set; corrected in the review of 2026-11-10)

| Score | AUROC | 95 % CI | Recall@100 |
|---|---|---|---|
| Target profile (expression up, tissue change, genetics, predicted tractability) | 0.85 | 0.80–0.89 | 12 % |
| Target profile without tractability | 0.75 | 0.68–0.81 | 14 % |
| Predicted tractability alone | 0.83 | 0.77–0.87 | 2 % |
| Target profile, tractability without ligand history | 0.83 | 0.78–0.87 | 12 % |
| Predicted tractability without ligand history, alone | 0.74 | 0.67–0.81 | 2 % |
| Biomarker profile | 0.79 | 0.73–0.85 | 7 % |

**Correction.** The first version of this table took predicted tractability only for the ~1,100 genes enriched
through the Open Targets API, which are mostly genes already associated with UC; every other gene counted as not
tractable. That favoured clinical targets and made "tractability alone" look as good as the profile. With the Open
Targets bulk table for every gene (same buckets as the kidney atlas), the target profile matches tractability on
AUROC (0.84 vs 0.83) but recovers far more clinical targets in its top 100 (14 % vs 2 %), and without the
medicinal-chemistry ligand bucket it beats tractability clearly (0.81 vs 0.75). As in the kidney atlas, the gain is
in early enrichment. Convergence v1 keeps its original tractability stream (enriched genes only, a known bias that
inflates its 0.88); it is not changed because the published narratives were written from it.

Disease activity: many genes track the Mayo endoscopic score in MSCCR (one value per patient, n = 375; for example
S100A8 rho 0.49, LCN2 0.33) but none reach significance in GSE73661, where baseline scores are mostly 2–3 and n = 67.
Mendelian randomisation with plasma pQTL: see (i). The temporal question is answered by the wider set in (b).

## (f) M17: repeats and one agent vs the team

Run on the five-cohort data of 2026-11-09 (before the layers of 2026-11-10); kept as the architecture comparison.

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

Every statement of the published run (runs/2026-11-10, 117 statements) that passes verification was copied and
corrupted one way at a time (`uv run python -m ucatlas.sabotage 2026-11-10`, evals/verifier_sabotage.json): a wrong
digit in a number 76/76 caught, the statement given to another gene 81/81, another gene's evidence 81/81, the
direction word flipped 23/72 overall and 23/32 (72 %) when the statement cites significant bulk expression evidence
for its gene. (On the 2026-11-09 run: 78/78, 83/83, 83/83, 16/25.) The direction check is narrow by
design (it reads only significant per-cohort and meta-analysis changes); a flipped direction in a sentence about a
paper or a cell type is not detected. Two fixes came from this test: the direction check now also reads the
meta-analysis records (before, 1 of 55 flips was even eligible), and ordering words ("moved down") and "up to" are no
longer read as expression directions.

## (h) What the agentic layer adds (closed-book baseline, review of 2026-11-10)

Claude Opus 5.5 without tools or data, asked the same question three times (`uv run python -m ucatlas.baselines`,
evals/baselines.json), against the three team runs:

| | Closed-book LLM | GraphRAG pipeline |
|---|---|---|
| Mean pairwise Jaccard of the lists | 0.72 | 0.76 |
| Candidates that are UC clinical targets | 70 %, 60 %, 70 % | 21 %, 21 %, 21 % |
| Candidates in the data-only top 500 (target or biomarker profile) | 50 %, 60 %, 50 % | 100 %, 100 %, 100 % |
| Statements citing checkable evidence | 0 | 120, 116, 119 checked; 0 unsupported |

The closed-book model lists the drugs already used in UC (TNF, IL-23, α4β7, JAK, S1P, TL1A); the pipeline lists
candidates the data support, most of them not yet drug targets in UC, each with checkable evidence. The same pattern
as the kidney atlas.

## (i) Mendelian randomisation with plasma pQTL (added 2026-11-10)

`uv run python -m ucatlas.mr` (data/processed/mr_pqtl.parquet). UKB-PPP cis-pQTL instruments (Open Targets 26.09
fine-mapped credible sets) on the de Lange 2017 UC GWAS; design fixed before measuring (history.md).

- 1,693 proteins have at least one usable cis instrument (8,013 instruments); **64 are significant (BH q < 0.05)** and
  17 of those colocalise with a UC, IBD or Crohn's GWAS signal in Open Targets (coloc H4 ≥ 0.8 or CLPP ≥ 0.01).
- **Positive controls.** IL12B, the target of ustekinumab, is the strongest result: OR 1.54 per SD of plasma protein
  (q 1e-21, 7 instruments, colocalises with UC). IL10 is strongly protective (OR 0.07), as expected from monogenic
  IL-10 deficiency colitis. CD274 (PD-L1) is protective (OR 0.60, colocalises with UC), consistent with colitis under
  checkpoint inhibitors. CD6 (target of itolizumab) and IL1RL1 (IL-33 receptor) are risk-raising and colocalise.
- **Candidates of this atlas.** CD274 and TNFRSF6B (OR 0.58, colocalises) are in the target profile's top 10; both
  are protective, so for them activation, not inhibition, is what genetics supports. OSMR has a borderline protective
  estimate without colocalisation (OR 0.88, q 0.04). CXCL8, IDO1, MMP1, CHI3L1 and ITGB7 are tested and null; CXCR2,
  CXCR1, OSM, DUOX2, CCL2, JAK2, TNF and IL23R have no cis-pQTL in UKB-PPP (not tested, not evidence against them).
- **Against UC clinical targets.** Only 15 of the 59 are tested; IL12B is the only one significant. AUROC of
  -log10 p among tested proteins 0.60. MR is a strong check for the few targets it can test, not a ranking.
- Limits: plasma protein is not mucosal protein; one assay platform (Olink); European ancestry only; a significant
  estimate without colocalisation can come from a nearby causal variant.

## (j) Spatial coverage (added 2026-11-10)

- Xenium 5K: two denser sections from GSE339323 (UC resection, most and least inflamed region; 283,086 and 188,644
  cells, median 293 and 259 transcripts per cell, 72 % and 70 % typed) join GSE312420. Four of six sections now pass
  the 50 %-typed rule (healthy, non-inflamed and two inflamed); the two sparse GSE312420 sections (7 and 44 transcripts
  per cell) stay out and are named in each record.
- Visium (GSE189184, whole transcriptome): 7 inflamed UC and 2 healthy colon sections cover the candidates that are
  not on the Xenium panel. The compartment that co-varies most with each gene matches the single-cell attribution for
  the known cases: S100A8, CXCR1, CXCR2, OSM, AQP9 neutrophil; LCN2, DUOX2, CXCL1 epithelial; MMP1, MMP3, IL11, OSMR,
  CHI3L1 fibroblast (tests/test_limits.py). Spots hold several cells, so this is a compartment, not a cell.

## (k) Meta-analysis coverage (added 2026-11-10)

With PROTECT, RISK and GSE165512 added, 13,840 genes are measured in all eight cohorts. Genes measured in one cohort
fell from 4,948 to 3,610 and no longer get a pooled significance (before: 2,522 of them were called significant on
MSCCR alone); the rest are mostly pseudogenes, small RNAs and lowly expressed genes. The discovery changes of MSCCR
replicate in the same direction in 95 % (PROTECT), 94 % (RISK) and 50 % (GSE165512) of tested genes.

## (l) Gut eQTL Mendelian randomisation (added 2026-11-10, second round)

`uv run python -m ucatlas.mr_colon` (data/processed/mr_colon.parquet). Colon tissue (GTEx v10 transverse and sigmoid
colon, CEDAR transverse colon and rectum): 5,852 genes tested, 135 significant (BH q < 0.05), 38 colocalise with an
IBD/UC GWAS signal, among them IRF5, GSDMB, RGS14, NXPE1, IRF8, GNA12 and the RTEL1-TNFRSF6B locus. Gut cell types
(IBDverse): 6,895 genes, 161 significant, 60 colocalise. CHI3L1, IDO1 and MMP1 are null in colon tissue; CXCL8, CXCR2,
OSM, DUOX2 and TNF have no gut instrument. For the 509 genes tested in plasma and colon, the estimates have the same
sign in 61 %.

## (m) Visium deconvolution (added 2026-11-10, second round)

`uv run python -m ucatlas.spatial_deconv`. 16,699 spots (12,989 UC). Held-out candidates (never markers) correlate
with the share of their single-cell group: LCN2 0.61, DUOX2 0.46, MMP3 0.23, CHI3L1 0.24, S100A8 0.18, IL11 0.17,
CXCR1 0.08, CXCR2 0.07. In-situ expressers: S100A8, CXCL8, AQP9, IL1B neutrophil; MMP1, MMP3, MMP10 inflammatory
fibroblast; LCN2, DUOX2, REG1A epithelial (tests/test_limits.py). Neutrophil share 1.7 % of transcripts in UC spots vs
0.9 % in healthy colon.

## (n) Pocket score against UC small-molecule targets (added 2026-11-10, second round)

46 AlphaFold models. The 9 modelled controls rank 4 (JAK1) to 25 (PDE4B) of 45 proteins by best pocket score.
Detected cavities cover a median 53 % of experimentally ligand-bound residues (PDBe-KB, 28 proteins), the top three
pockets 24 %: the detector finds real sites better than the score ranks them first.
