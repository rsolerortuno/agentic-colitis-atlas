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

## (o) Independent checks (added 2026-11-10, third round)

- Protein (Schniers 2019): 373 clear RNA changes testable, 281 replicate as protein (75 %; 77 % for |log2FC| >= 2);
  rho 0.61 over 6,112 genes.
- FinnGen R12: 34/64 significant plasma proteins and 168/296 significant gut estimates replicate (same sign, p < 0.05);
  IL12B, IL10, CD274, TNFRSF6B, IRF5, GSDMB, NXPE1 and RGS14 all replicate.
- Mouse models (MGI): 236 genes with a colitis phenotype (1.2 % of genes with an ortholog); UC clinical targets 24 %,
  convergence top 100 15 %, target-profile top 100 12 %.
- CosMx 6K: 22/108 candidates above background; 1K: 3/50.
- Treatment: ACT1 infliximab 11 genes at q < 0.05 (7 non-responders vs 25 responders), placebo 0; EUCALYPTUS
  etrolizumab 0; no drug x placebo interaction passes correction (ustekinumab, golimumab, infliximab).
- Weight stability: biologics 7/10 top-10 stable (rho 0.97), immune 9/10 (0.99), pockets 8/10 (0.99), novelty 1/10
  (0.91).

## (p) Closing checks (added 2026-11-10, fourth round)

- Covariates: rho 0.989 between adjusted and unadjusted MSCCR changes; 1,086/1,086 clear changes hold.
- MR: Liu 2023 replicates 58/64 plasma and 276/296 gut estimates. Plasma sensitivity: 48 with >= 3 instruments,
  Egger intercept p >= 0.05 in 43, weighted median agrees in 45, heterogeneous in 23 of 64; Steiger correct 64/64,
  135/135 colon, 161/161 gut cell.
- Safety: strongly constrained candidates INHBA and GJC1; clinical targets constrained 5/59 (JAK1, S1PR1, NR3C1,
  PTGS2, NCBP1). PheWAS: IL12B psoriasis and mouth ulcers; ITGA4 white cell counts; SLC26A3 triglycerides.
- Negative control: expression-matched AUROC convergence 0.856, target profile 0.815, biomarker 0.764 (random sets
  95th percentile 0.56-0.57); literature-matched 0.729, 0.702, 0.711.
- PROTECT: corticosteroids 1,263 genes at q < 0.05 (78 NR / 75 R), 5-ASA 0 (23 / 30); calprotectin (145 children):
  AQP8 rho -0.29 and HMGCS2 -0.26 significant, candidates up in UC ~0.2 and not significant.
- Anti-TNF atlas: 6 remission / 14 non-remission UC patients; DUOX2 and LCN2 lower in non-remitter colonocytes before
  adalimumab and higher after treatment.
- Mouse literature: 29 of 69 genes tested in mouse colitis; 0 of the 36 first-in-class genes.

## (q) Final agent benchmark (2026-11-10, closing data)

| Set-up | Run | Statements | Unsupported | Checklist | Evidence cited | Cost | Minutes |
|---|---|---|---|---|---|---|---|
| Team | published | 119 | 0 | 10/12 | 90 % | $4.61 | 7.5 |
| Team | rep1 | 116 | 0 | 11/12 | 96 % | $4.76 | 7.4 |
| Team | rep2 | 124 | 0 | 10/12 | 98 % | $4.99 | 8.0 |
| One agent | single | 141 | 0 | 10/12 | 82 % | $7.52 | 19.3 |

Target-list overlap: team runs 0.87-1.0; single agent vs team 0.56 (it leads with IL12B, TNF and ITGB7).

## (r) Validation shortlist (added 2026-11-11, revised the same day)

The 36 first-in-class genes, each with a light per criterion (rules in `src/ucatlas/validation_shortlist.py`).
Structure now covers every first-in-class gene (72 proteins modelled).

**First version and why it changed.** It put CCNO (cyclin O) first. Read gene by gene, three criteria were too
lenient: (1) CCNO sits in transit-amplifying crypt cells and rises per cell there, like OLFM4 and REG1A: it is part of
the regenerating crypt, which the lineage-level composition adjustment cannot see; (2) its loss causes primary ciliary
dyskinesia (UniProt), which neither LOEUF (recessive disease) nor vital-organ RNA (few ciliated cells) catches;
(3) "surface or secreted" counted any membrane location, so cytoplasmic or ER proteins (OSBPL3, EHD3, RFTN1, SGIP1,
SDR16C5) got an antibody route. Also, a null genetic test scored like a partial success.

**Revisions, each checked on markers or controls before use.**
- Per-cell change inside the gene's own single-cell clusters (Smillie 2019). Markers: DUOX2 +7.1, MMP3 +4.9,
  CHI3L1 +4.1 log2 (state changes), TNF +0.2, IL1B -0.5, OSM -0.6 (influx), AQP8 -2.7. CALHM6, CLEC1A, IKBIP,
  RILPL2, FBXO16 and FILIP1L do not rise per cell (red).
- Crypt-zone regeneration flag: >= 76 % of epithelial RNA in stem/TA/progenitor clusters (EPCAM 58 %, LGR5, OLFM4,
  MKI67, REG1A 92-97 %) and rising there: CCNO and CYP4X1.
- Mendelian disease (UniProt, only "caused by variants", no susceptibility loci): CCNO, TPK1 and CSGALNACT1 among the
  new genes; 23 of 59 clinical targets (mostly immunodeficiencies, the intended effect of immunosuppression). Kept as
  a separate evidence record so the target-safety records cited by the agents stay unchanged (0 stale).
- Antibody route from UniProt topology (extracellular domain, GPI anchor, secretion). Controls: all 16 clinical
  targets in the table keep a green druggable light; 58 % of the new genes do.

**Calibration on the 16 clinical targets** (share of green lights, first-in-class vs clinical): reproducible 94 vs
100 %, protein 14 vs 19 %, tissue 42 vs 12 %, cell 56 vs 69 %, genetics 0 vs 0 %, druggable 58 vs 100 %, safety 89 vs
19 %; median points 8.5 vs 8. Safety is an ascertainment effect (side effects and diseases are known for studied
genes), and tissue change describes the kind of target the screen looks for, not its odds; both are said in the tab.

**Result.** First ten: ERO1A, CHRDL2, CPXM1, PCDH17, LYPD5, LAIR2, OLFML2B, ADGRE3, OSBPL3, FIBIN; CCNO falls to 23rd.
Eight of the ten are secreted or surface proteins of stroma, endothelium, epithelium or immune cells that rise per
cell; ERO1A (ER oxidoreductase with known inhibitors) leads on protein replication and per-cell rise. CCNO's
"high-quality ligand" in Open Targets comes from the CDK2-cyclin O complex, not from CCNO itself. The table is a way
to read the other layers, not a new estimate; genes with equal points are equally supported.

**Literature review of the candidates (added 2026-11-11).** Europe PMC: every record naming the gene (symbol or an
unambiguous alias) with UC anywhere in the text, plus records with the gene in the title or abstract and colitis, IBD
or Crohn's disease: 475 records for the 36 genes. Only 62 name the gene in the title or abstract and 4 together with
UC there; for open-access papers the full-text sentences naming the gene were added, leaving 183 readable records
(292 name the gene only in a full text that is not open). Eight classifier subagents (Sonnet, subscription) labelled
each record (kind, system, direction) with a verbatim sentence: 183 kept, 0 rejected by the check. A paper counts
only in a colitis or IBD setting (deterministic text check): 63 of 183, 60 of them about the gene. Result: 24 genes have nothing specific in
the colitis literature (only gene lists or no paper), 10 a specific expression, biomarker or genetic result, one a
functional study supporting inhibition (RFTN1: Raftlin1 promotes pathogenic Th17 cells in C. rodentium colitis,
PMID 37591835) and one a functional study against it (IKBIP: IKIP-deficient mice develop worse DSS colitis, PMID
31826938; the earlier title search for mouse studies missed it because it used the symbol only). The paper count of
the discovery screen (median 8 per gene) therefore overstates what is known: most hits name the gene in a
supplementary table or are about another disease. New order: ERO1A, CPXM1, ADGRE3, RFTN1, CHRDL2, PCDH17, LYPD5,
LAIR2, OLFML2B, SGIP1; IKBIP falls to 32nd.

**Groups instead of positions (added 2026-11-11).** Review of the order: many genes tie on points, and the order inside
a tie came from tie-breakers (greens, fold change), not from evidence. The table now gives the stability of each gene
across six scorings (base; without the protein light; without the literature light; without both; red lights costing
2; an unmeasured protein counted as amber), as the number of them that place it among the first ten, ties at the border
included. The rule was fixed before it was applied. Group A, all six: ERO1A, CPXM1, ADGRE3, CHRDL2, PCDH17, LYPD5, LAIR2,
OLFML2B, SGIP1. Group B, one to five: 13 genes, among them RFTN1 (5 of 6, but its position ranges from 2nd to 23rd:
it leans on the protein light and on a single functional paper), OSBPL3 (3) and CCNO (1). Group C, none: 14 genes.
**Protein bias:** 7 of the 36 candidates have a quantified protein, and across all genes those quantified sit at the
78th percentile of RNA abundance in inflamed mucosa against the 37th for the rest, so a measured protein earns points an
unmeasured one cannot; the tab says so and the "without protein" schemes show each gene's dependence.
**Experiments by concrete cell system:** the broad lineage sent vessel, macrophage and Treg genes to the wrong culture.
The first experiment now follows the cell groups that hold the gene's RNA in inflamed single cells (a second system
when it holds at least 20 %): PCDH17 endothelium (87 %), OLFML2B macrophages and DC (44 %) plus fibroblasts (29 %),
LAIR2 regulatory T cells (48 %) plus T and NK cells (41 %), RFTN1 B and plasma cells (49 %) with the published Th17
lead, ERO1A organoids (33 %) plus B and plasma cells (23 %), CCNO organoids with a regrowth assay (92 % in crypt
progenitors).

## (s) Crypt-niche flag and a fourth single-cell atlas (added 2026-11-11)

**Crypt-niche flag (stroma).** The fibroblasts under the crypt base (Smillie WNT2B+ Fos-hi, WNT2B+ Fos-lo 1 and 2,
RSPO3+; not the Inflammatory Fibroblasts) hold 3,594 of 7,239 fibroblasts in inflamed mucosa, so a generic fibroblast
gene has about half of its fibroblast RNA there: COL1A2 52 %, COL3A1 50 %, VIM 50 %, DCN 59 %. Niche markers: WNT2B
62 %, GREM1 70 %, RSPO3 83 %. Threshold 62 %, halfway between COL1A2 and the niche-marker mean (72 %), set as the
crypt threshold was (EPCAM vs LGR5/OLFM4). Flag = niche share >= 62 % and a per-cell rise >= 1 log2 inside the niche
clusters. **No shortlist gene is flagged.** Group A stromal genes: CHRDL2 48 %, CPXM1 52 %, SGIP1 54 %, OLFML2B 50 %
(all rising per cell in the niche, +2.0 to +3.4 log2, but with the niche share of a generic fibroblast gene) and
PCDH17 13 % (an endothelial gene). The flag would cap the cell light at amber and add an organoid-fibroblast
co-culture read-out to the experiment.

**GSE282122 as a fourth atlas.** The filtered 10x matrices (2.8 GB tar, sha256 in data/raw/gse282122/manifest.json)
have no cell labels on GEO. 112 colon biopsies (103 UC before and after adalimumab, 9 healthy; ileum excluded),
694,801 cells after QC, typed by k-means on a canonical marker panel (60 clusters) into nine lineages. Check: every
marker of the eight lineages found is highest in its own lineage; no neutrophil cluster (droplet filtering loses
them); 51,036 cells in four unclear clusters left out. The series mixes two Cell Ranger references, matched by Ensembl
id. Tables per gene and lineage in the sc_celltype_share format, also by treatment time and remission
(sc4_celltype_share, sc4_response_share).

**Replication across four atlases** (five shared lineages, inflamed UC, 1,003 genes): the Smillie lineage is dominant
in all three other atlases for 512 genes, in two for 299, in one for 156, in none for 36. Per atlas: Garrido-Trigo 804
of 995, Gubatan 635 of 878, GSE282122 851 of 997. The cell-source light now needs >= 2 of 3 for green. Effect on the
shortlist: one change, OLFML2B (stromal in Smillie and GSE282122, myeloid in Garrido-Trigo and Gubatan) goes from
green to amber and from group A to B. New graph evidence types sc_lineage_gse282122 and sc_replication_4atlas; the
records the agents cite are unchanged (verify final: 0 unsupported, 0 stale).

**Agent run with four atlases (runs/2026-11-11).** The cell tool now returns the lineage in each atlas and the
two-of-three rule, so the six specialists and the orchestrator were re-run on it. Verification: 116 statements,
0 unsupported, 0 stale, 28 of 28 key findings cited (first pass dropped 33 of 73 agent statements for formatting of
numbers; fixed by the agents from the verifier output). Content changes: CXCR2's Smillie lineage (T/NK) is confirmed
by only the anti-TNF atlas, Garrido-Trigo and Gubatan give myeloid; MMP1, DUOX2 and HMGCS2 replicate in all three
other atlases, CHI3L1 and IDO1 in two. Leading targets unchanged (CXCR2 first, IL12B as positive control); CHRDL2 and
CPXM1 are the exploratory first-in-class targets.

## (t) Neutrophil-made genes in the shortlist (added 2026-11-12)

Review of group A found ADGRE3 scored as a macrophage gene: Smillie 2019 has no neutrophils, so its per-cell rise
(+2.4 log2) was measured in macrophages and its lineage "replicated" in three atlases that also lack neutrophils, while
Gubatan 2026 puts 82 % of its RNA in neutrophils. Rule (fixed before applying): a gene with more than half of its
inflamed-mucosa transcripts in neutrophils (Gubatan) is neutrophil-directed; its tissue light uses the per-cell change
inside the Gubatan neutrophil states (inflamed vs uninflamed UC, states with >= 50 cells, same thresholds), its cell
light is amber (one atlas only, cannot replicate) and its first experiment is in blood neutrophils. Checks: S100A8
+1.1, CXCR2 +0.7, FCGR3B +0.7 log2 per neutrophil; ADGRE3 -0.3, MMP9 -0.2, PTGS2 -0.2 (red: rise by arrival).
Shortlist genes affected: ADGRE3 (A to C), CEACAM4 (B to C), CLEC1A, CSGALNACT2, RILPL2 (stay C). With ADGRE3 out of
the top ten, the 9-point ties fit inside it under all six scorings, so group A grows from 8 to 13 genes (FIBIN, XKR9,
FJX1, OLFML2B, CYP4X1, SDR16C5 join without any change of their own points). The cell family column also shows
"endothelial" for vessel genes (PCDH17). These genes remain valid hypotheses of another kind: blocking neutrophil
arrival, like CXCR2.
