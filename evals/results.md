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

Update of 2026-11-23 (run 2026-11-14, 115 statements): numbers 84/84, gene 80/80, evidence 80/80, but direction only
14/53 (26 %), 14/21 with bulk expression evidence. Of the 39 flips missed, about 28 cited evidence with a sign that the
check did not read (pathway and transcription-factor activity, treatment response, Mendelian randomisation) or were
hidden by an ordering word elsewhere in the sentence ("tops the ranking: log2FC is down"); about 11 were literature
sentences. The check now reads all of these (the literature through the direction of the verbatim sentence cited, when
it speaks of the gene's expression and not of an experiment on it): direction 44/53 (83 %), 18/21 with bulk evidence,
and none of the 115 real statements is newly flagged. A negation now counts only before the direction word ("most
reduced ... without an inhibitor" is not negated), and a clause naming several genes is checked when it contradicts
the evidence of every gene it names; a clause naming none is about the gene of the previous clause: direction 46/53
(87 %), 19/21 with bulk evidence, still no real statement flagged. Still missed: flips of ordering words ("moved
down", "lower-ranked"), a sentence made self-contradictory ("increased (lower)"), a pathway the statement cites no
evidence for, a comparison with the rest of the shortlist ("fewest cohorts"), and comparisons that are not
expression directions (AUC, mouse versus human).

## (g2) Reading check (2026-11-23)

The deterministic checks cannot tell whether a claim is contained in what is cited. claim_check.py gives every
statement that passes verification, with the full records it cites and the abstracts of cited papers, to a reader
subagent (.claude/agents/lector-afirmaciones.md; this first time run as general-purpose agents on the same model and
instructions, because a new agent type loads only in a new session). The reader splits each statement into claims and
calls each supported, contradicted or not in the evidence, with a quote that must be verbatim in the evidence (all
were). The 53 direction-flipped copies were mixed in under opaque ids.

- Real statements: 606 claims, 3 contradicted, 130 not in the cited evidence. The 3 were checked by hand and are real:
  targets.IDO1 calls the IDO1 paper "a single mouse study" although it measured IDO1 in UC patients (GEO and western
  blot); narrative.sections.spatial says the neighbourhood of MMP1-high spots "was not reported in this run" although
  the cited Visium record gives it (inflammatory fibroblasts, log2 enrichment 1.13, and neutrophils, 1.06); and
  narrative.sections.pathways names "Neutrophil CXCR4" as receiver of CXCL8 to CXCR2, where the record says
  "Neutrophil". The 130 unsupported claims are about 23 comparisons across genes backed by one gene's record ("the
  largest", "the only"), 33 drug names or approval status, 9 gene functions, 58 results of other analyses of this atlas
  that the sentence does not cite (cell source, CosMx, MR sensitivity, other pathways), and 7 literature claims taken
  from an abstract or paper not cited (CXCL1 "five-gene model" and "several studies", S100A9 "nanopore RNA-seq",
  MMP10 "biomarker panel", IDO1 "single experimental study"). These are citation gaps, not wrong facts, except the
  three above.
- Sabotage copies: 51 of 53 flips caught (a flagged claim containing the flipped word); the two missed are the
  ranking-word flips ("moved up", "lower-ranked") that the instructions exclude. The deterministic checks catch 46.

The reader is a language model; its verdicts are shown in the report as a review aid and do not change the verified
evidence. Correcting the three statements needs a new agent run.

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

## (u) Drugs by indication and repurposing leads (added 2026-11-12)

The Candidates tab listed each gene's drugs without saying what they are for, and the discovery screen sets aside every
gene with a clinical drug. New layer drug_indications.py (Open Targets 26.09 clinical precedence; drugs with more than
five targets set apart): 111 genes have a drug in UC, 40 in Crohn's disease or IBD only, 874 only in other diseases,
490 only broad drugs. In the Top 25: CXCR2, CXCR1 and MMP1 have UC trials (phase 2); CXCL8, CD274 and IDO1 have drugs
only for other diseases. In the discovery screen, 52 genes set aside as "known druggable target" are repurposing
leads; since they rise in UC, an inhibitor is what counts: 35 have one in the clinic (Open Targets directionOnTarget
LoF), 17 have only activators or drugs whose action is not stated (often antibody-drug conjugates, e.g. TACSTD2).
Checks: an approval resting only on an ATC code is reported as not confirmed (epacadostat, IDO1, was never approved);
4,922 of 101,815 approval rows are ATC-only. New evidence type clinical_indications; verify final 0 unsupported, 0 stale.

## (v) Repurposing shortlist (added 2026-11-12)

The 35 discovery genes with an inhibitor in the clinic only for other diseases, scored with the shortlist lights.
Two lights read differently: the inhibitor's stage replaces "druggable" (16 approved), and safety adds a flag for
inhibitors made only for cancer (8 genes). Immune checkpoints were found in the first draft: CTLA4 reached group A on
the strength of ipilimumab, a drug that causes colitis; CTLA4 and CD274 are now red on both lights (they rise in UC as
brakes, and for CD274 genetics says the protein is protective). A cancer flag first caught GUCY1A1 because methylene
blue lists "colorectal neoplasms" (it is a colonoscopy dye); the flag now needs every indication to be a cancer.
Group A: CFB (iptacopan), SELP (crizanlizumab), GUCY1A1, NOS2, PCSK9, INHBA, LOXL2, C3 (pegcetacoplan), KDR, BIRC3,
ICOS, ANGPT2. The literature light is grey for all (not reviewed for this list), so it does not change the order.

**Literature review of the repurposing genes (added 2026-11-12).** The same classifier and checks as the first-in-class
review, run on the 35 repurposing genes: the 25 most cited papers per gene with the gene in the title or abstract and
colitis, IBD or Crohn's disease (828 records, 818 readable, 33 batches). 816 calls kept, 2 rejected (quote not in the
text). Lights: 9 green (functional evidence that blocking helps: CD40, F2R, INHBA, KDR, NAMPT, PCSK9, PLA2G2A, SELE,
TREM1), 6 red (protective: BIRC3, CTLA4, CTSK, IDO1, LYN, MME), 17 amber, 3 grey. Effect: BIRC3 leaves group A (its
loss worsens colitis), C3, LOXL2, ICOS and ANGPT2 move to B, SELE joins A. Group A: CFB, SELP, GUCY1A1, PCSK9, NOS2,
SELE, INHBA, KDR. KDR stays in A on a pathogenic colitis study; its inhibitors are cancer drugs with known GI toxicity
(safety amber), which the list states but does not score further.

## (w) Neutrophil genes without the composition adjustment; tops read for an expert (added 2026-11-12)

**Neutrophil-made genes.** The bulk composition adjustment removes the neutrophil influx that defines these genes, so
using it on top of the per-neutrophil change penalised them twice (and raised their green threshold to 1.5 log2). It
is no longer scored for them. A first version used the per-neutrophil change alone and moved CLEC1A from group C to A,
although 47 % of its RNA is in stroma and endothelium where it does not rise per cell; the change is now weighted by
where the RNA is (neutrophil share x per-neutrophil change + rest x per-cell change in Smillie 2019). Result: TREM1
green (+1.3, repurposing C to B), CEACAM4 green, CLEC1A amber (+0.6), ADGRE3 and RILPL2 red; MMP9 and PTGS2 red among
the clinical targets. First-in-class group A unchanged.

**Discovery spotlight.** It showed the top first-in-class genes by novelty score alone, including genes the
shortlist puts in group C (TPK1, a thiamine enzyme with a Mendelian disease; FBXO16, no papers and 26 % kept; ENKUR, a
ciliary protein). It now shows shortlist groups A then B per cell family: epithelial ERO1A, LYPD5, XKR9, CYP4X1;
stromal CPXM1, CHRDL2, PCDH17, SGIP1; immune LAIR2, RFTN1, OSBPL3, LAX1.

**Candidates Top 25.** A convergence ranking recovers known UC biology (CXCL8-CXCR2, MMPs, calprotectin, the lost
epithelial programme): that is its check. Each row now says from the data what kind of gene it is (made by
neutrophils, lost in UC, rise explained by arriving immune cells, immune brake, or tissue change, with where its
drugs are), so it is not read as 25 targets to block. Of the 25: 9 neutrophil-made, 4 lost in UC, 1 immune brake
(CD274), 1 explained by immune cells (CXCL5), 10 tissue changes (MMP1 and MMP7 already tried in UC; IDO1, CCL2 and
INHBA with drugs for other diseases).

## (x) Three target lists; neutrophil genes through the discovery screen (added 2026-11-12)

**Discovery.** The tissue-change filter removed every neutrophil-made gene before the other filters (CXCL8 retained
14 % of its change). Genes with more than half of their inflamed-mucosa RNA in neutrophils now pass it as
neutrophil-directed: 88 genes released; first-in-class 36 to 46 (C19orf38, CMTM2, CYTH4, GLT1D1, LILRA1, LILRA6, MCTP1,
PLEKHO2, RIPOR2, TMTC1), repurposing 35 to 42 (C5AR1, CCR1, CXCL8, IL1R1, LY96, SELL, TLR8). The neutrophil_directed
column is kept out of the discovery graph records; verify final stays at 0 stale.

**Well known in UC, no drug yet.** Genes of the Candidates top 100 that rise as a tissue change or are
neutrophil-directed, have more than 15 UC papers, no clinical drug and a predicted route (37 genes; vital-organ
expression is a safety light, not a filter, so AQP9 and CHI3L1 are in). Literature review as for the other lists
(the 25 most cited colitis/IBD papers per gene, 821 records, 810 calls kept, 1 rejected): 9 green, 7 red (PI3, OLFM4,
SOCS3, FPR1, FPR2, DMBT1, SAA2), 15 amber, 6 grey. Group A: SERPINB5, DUOXA2, CHI3L1, DUOX2, REG3A, KLK10.

**Literature for the added genes.** 10 first-in-class (173 records, only 54 readable: these genes appear mostly in
full text) and 7 repurposing genes (175 records) reviewed; all quotes but 2 pass the verbatim check. First-in-class
group A gains LILRA1; repurposing group A gains C5AR1 and LY96.

**Europe PMC retrieval.** Pages sized to the cap (a 100-record "core" page for 25 records), one HTTP session, the
HGNC table read once, eight genes at a time, and full-text articles fetched by one shared pool of eight instead of
one after another per gene: the 10 added first-in-class genes took 385 s, against more than 28 min unfinished before.
Ambiguous aliases found by the classifiers (LIND, TC1, TCI, ACT, LI2, MIF, PI5, GP39, FSP) are excluded from future
searches.

**Mucosal defence flag (added 2026-11-12).** Read against UC knowledge, the top of the well-known list held defence
effectors (DUOX2 and DUOXA2, the hydrogen-peroxide system whose loss-of-function variants predispose to very-early-onset
IBD; REG3A, an antimicrobial and repair peptide): the lights rewarded their per-cell rise, which is the protective
response itself. Flag from UniProt GO biological process, narrow terms only (hydrogen peroxide production,
antibacterial response, killing of microbes, defence response to bacterium, wound healing); broader terms were rejected
because they mark IL23A, IL12B and CXCL8. Effect on every list: a safety flag and the tissue light capped at amber.
Calibration: no UC clinical target is flagged. Flagged: 9 well-known genes (DUOX2, DUOXA2 and REG3A leave group A),
NOS2 and CASP1 (repurposing), none first-in-class. SERPINB5 stays first on the well-known list: the review found two
2026 functional studies in DSS colitis in which its knock-down protected (not yet independently replicated).

**Safety coverage (fixed 2026-11-12).** The safety step covered 130 genes (report candidates, convergence top 40,
first-in-class at the time, clinical targets), so 23 of 37 well-known, 10 of 46 first-in-class and 35 of 42
repurposing genes showed "not constrained, no Mendelian disease" without having been checked (ABCA12, whose loss
causes congenital ichthyosis, among them). It now covers 545 genes (every discovery gene the lists can hold and the
convergence top 100); a gene outside it would read "not checked". Effect: ABCA12 leaves well-known group A; CFB
(complement factor B deficiency), KDR and C5AR1 leave repurposing group A, PLA2G2A and IL1R1 join it. The target_safety
records cited by the agents are unchanged (0 stale).

**Two errors in repurposing group A (fixed 2026-11-12).** (1) PCSK9 carried "loss causes familial hypercholesterolaemia
3": that disease is caused by gain-of-function variants (loss lowers LDL, and PCSK9 inhibitors are approved), but
UniProt does not say so in the text the parser reads. A curated list of gain-of-function diseases (PCSK9, LYN, TEK
venous malformations, TLR8) and of a somatic cancer listed as a disease (IGHG1, multiple myeloma) now corrects the
Mendelian flag. (2) IL1R1 was listed as never tried in UC, but anakinra was tested in acute severe colitis (IASO,
protocol PMID 30772849), a trial in neither Open Targets nor ClinicalTrials.gov. Drugs of list genes are now also
searched in ClinicalTrials.gov and in Europe PMC trial papers. The first version found six genes, five of them false
(ipilimumab in immunotherapy-induced colitis, methylene blue as a chromoendoscopy dye, a suPAR biomarker paper, other
IL-18 antibodies matched by a fuzzy registry search); with stricter filters (the drug among the study's
interventions, a UC/Crohn's/IBD condition that is not cancer, not a diagnostic study; trial titles naming UC, Crohn's,
IBD or severe colitis, immunotherapy colitis and biomarker papers excluded) only IL1R1 remains.

## (y) AlphaGenome: predicted effect of UC risk variants (added 2026-11-12)

**Variants.** Liu et al. 2023 UC GWAS (GRCh38), SNVs with p < 5e-8 outside the MHC, 500 kb clumping: 142 loci.
Wakefield approximate Bayes factors (prior SD 0.2, one causal variant per locus): 2,381 variants in the 95 % credible
sets. Per locus and gene, effect = PIP-weighted mean of the predicted effects, signed for the risk allele.

**Validation (fine-mapped single-variant eQTL, PIP >= 0.9, as truth).**

| Check | Pairs | Same sign, all | Same sign, top quarter of abs(score) | Spearman |
|---|---|---|---|---|
| Atlas colon tissue vs GTEx v10 colon | 284 | 64 % | 94 % (threshold 0.033) | 0.41 |
| AlphaGenome 2, matched colon cell type vs IBDverse | 236 | 70 % | 100 % (threshold 0.169) | 0.60 |
| AlphaGenome 2, mismatched cell type | 236 | 56 % | 75 % | 0.26 |

**Result.** 7 loci and 17 genes with a confident effect (Atlas: colon tissue, blood T cells, blood monocytes).
Genes of interest (three lists, Candidates top 100, clinical targets): CCL20 (+0.24 T cells, +0.14 monocytes, top
variant rs1811711) and IL10 (-0.09 T cells; consistent with the protective MR estimate). Others: TNFSF15 (-0.25
monocytes), TNRC18 (-0.50 T cells), IL19, SLC29A4 and non-coding genes. AlphaGenome 2 in 12 colon cell types (137
variants near genes of interest, 1,644 predictions): no gene of interest above threshold (largest 0.069).

**Decision.** Shown in the Genetics tab with its evidence records (alphagenome_variant_effect); it does not change any
score or list light, because only two genes of interest are covered and both already have genetic evidence.

## (z) Genetics-first list (added 2026-11-13; corrects (y))

**Risk signals.** Open Targets 26.09 credible sets (p < 5e-8, SNVs, MHC excluded) of five UC GWAS and, when one of
their variants has p < 1e-4 for UC, three IBD GWAS; cross-study duplicates merged (shared variant, or leads within
250 kb in different studies). Loci of the Lyu 2025 European UC meta-analysis (preprint; 136 loci, 32 outside the 142
of Liu 2023) that no credible set covers are fine-mapped approximately. Total 263 signals, 2,375 variants: Liu 2023
UC 148, Liu 2023 IBD 52, MVP 22, Lyu 2025 only 22, FinnGen UC 9, de Lange IBD 4, de Lange UC 3, FinnGen IBD 3.

**Predictions.** Atlas for every variant with PIP >= 0.01 (up to 50 per signal); AlphaGenome 2 for the three most
likely variants with PIP >= 0.1 (352 variants) in 15 conditions: 5,280 predictions.

**What the size threshold does not test.** Around the eQTL variants of the validation, 2.4 % of the other
protein-coding genes pass the cell-type threshold 0.026 (0.6 % at 0.169): about 0.35 genes per variant and cell type.
A list built on size alone had 316 genes (28 confident); its direction agreed with significant MR in 3 of 7 confident
and 12 of 25 probable genes.

**Gene-assignment rule and its check on fine-mapped eQTL.**

| Rule | Tissue (GTEx colon) | Cell type (IBDverse, matched) |
|---|---|---|
| Top-scoring gene of the variant | measured eQTL gene in 38 % | 47 % |
| Top gene, at least twice the second, above the size threshold | 88 % of 48 (direction right in 98 %) | 93 % of 76 (direction right in 94 %) |

**Result.** 6 of 263 signals get a gene. Group A: TUFM (risk allele raises it, strongest in neutrophils and
monocytes; MR in sigmoid colon agrees, OR 1.19) and RORC (risk allele lowers it in monocytes; also lower in inflamed
biopsies). Group B: IL19 (lowered; same variant as IL10, rs3024493) and MAP3K8 (lowered in T cells). Likely, with an
uncertain causal variant: CSPG5 (IBD signal) and LSP1 (MR disagrees). None is on the other three lists. IL10 and
CCL20, named in (y), are not assigned under this rule.

**Decision.** The list is published as a fourth tab with this result stated first; it changes no score. The per-gene
alphagenome_variant_effect records of (y) are no longer loaded.

## (g3) Every claim cited, and the reading check inside the agent run (2026-11-23)

The reading check of run 2026-11-14 found 130 of 606 claims not in the evidence their statement cites. Three changes
followed: (1) citable records for what was missing: `rank_context` (2,655 records: a gene's rank within named sets,
for comparisons across genes), UniProt function and location and Reactome pathway names inside their records, and
the sentence text inside each abstract-sentence record; (2) a deterministic rule: a drug named in a statement must
appear in a cited record (it flagged 15 statements of run 2026-11-14); (3) the specialists' and the orchestrator's
instructions ask for a record per claim, and /informe runs the reader on the specialists' verified findings and on
the narrative, sending what it flags back once before verification.

Run 2026-11-23 (all six specialists and the narrative redone): the first verification dropped 30 of 80 findings
(unbound numbers, uncited drug names, numbers bound to another gene's record) and the reader flagged 28 of 182 claims
of the kept findings; after one revision 79 of 80 findings pass (a CD274 literature finding whose citations are not
linked to the gene was removed). The reader flagged 10 claims of the narrative, two of them wrong (HMGCS2 called the
largest loss although it ranks third; CD274 called lineage-less although the record names myeloid), corrected before
assembly. Final: 124 statements, no failures detected by the automatic checks; 473 claims read, none contradicted,
27 (6 %) not in the cited evidence (23 in specialist findings, 4 in the narrative), listed in the report. The new
`rank_context` records also showed that two earlier claims were wrong: XKR9 is not the only intrinsic epithelial
first-in-class candidate (four are) nor the one with the fewest cohorts (second fewest). Sabotage on the new run:
numbers 77/77, gene 94/94, evidence 94/94, direction 62/73 deterministic (33/34 with bulk evidence) and 72/73 by the
reader. The candidate SLC26A3 is new in this run, so its AlphaFold model and mouse literature were added (87
proteins; pocket controls rank 7th to 46th of 85; ligand-site recall median 53 % of 36 proteins).

## (g4) The top of the candidate list read as a UC clinician would (2026-11-26)

Asked whether IL12B, CXCR2 and JAK2 at the top of the list would look wrong to UC experts. Each text was verified,
but the presentation would: (1) the drug tag showed the most advanced drug for any disease, so JAK2 carried
fedratinib (a myelofibrosis drug) and IL12B guselkumab (an anti-p19 antibody that Open Targets lists under both
subunits of IL-23); (2) the list opened with a neutrophil receptor and a JAK2 card read as proposals; (3) three
sentences were wrong for a specialist: "IL23A, the target ustekinumab blocks" (ustekinumab binds IL12B), "the UC
drugs listed for IL12B" without saying they act on IL-23 as a complex, and CXCR2 antagonists "none for UC", whereas
the clinical-indications record lists SB-656933 and SCH-527123 (navarixin) in phase 2 in UC.
Changes: the drug tag shows the UC drugs first, else Crohn's disease or IBD, else the drug with its disease, and names
the other genes a UC drug is listed for (IL23A for IL12B; JAK1, JAK3 and TYK2 for JAK2); the candidates are shown in
four blocks set by rules on the data (UC-approved drug; most RNA from neutrophils; a tissue change; biomarkers); the
narrative calls IL12B the positive control of the IL-12/IL-23 axis and JAK2 the positive control of the JAK class,
with the reason (the human genetic support of the family is on JAK2, upadacitinib is listed for all four kinases,
and JAK2 carries erythropoietin and thrombopoietin receptor signals), uses IL12B, the drug's own target, for the
ustekinumab null (padj 0.84), and gives the CXCR2 phase 2 trials in UC (narrative and the external specialist's
finding). Changed statements were verified and read again (reader lines reused for unchanged statements): 125
statements, no failures detected by the automatic checks; 499 claims, none contradicted, 16 not in the cited
evidence. The specialists' own findings about the three genes were also cited or reworded where the reader found a
claim without its record (ranks, drug records; "Smillie 2019 captured no neutrophils" became the record's own
figure, at most two percent of the cells of its leading clusters express CXCR2): no claim about IL12B, JAK2 or
CXCR2 is left contradicted or uncited.

## (g5) Independent expert-style reviews of the top candidates (2026-11-26)

Three reviewer subagents, each prompted as a senior IBD gastroenterologist and each new to the report, read the first
three candidate cards (IL12B, JAK2, CXCR2) and every report sentence about them. The first two found real objections,
all fixed before the next reading:
- first review: anti-p19 antibodies shown as IL12B drugs, upadacitinib's JAK1 selectivity unsaid (and tofacitinib and
  filgotinib missing), the positive-control logic partly circular (the Open Targets overall score includes drug
  evidence: the cards now use its genetic component), CXCR2 trials and the neutrophil drop of the class unsaid, the
  thin CXCR2 literature read as a field gap, and Smillie 2019 not explained (it has no neutrophil cell type);
- second review: risankizumab and guselkumab missing, "reached phase 2" read as progress (now "went no further than
  phase 2"), navarixin under three names, why only two positive controls, the TYK2-selective deucravacitinib failure in
  UC, and a published Mendelian randomisation in which higher CXCR1/CXCR2 expression went with lower UC risk;
- third review: "acceptable to UC experts", with five minor points, applied (IL-23 sufficiency rather than IL-12
  irrelevance, tofacitinib also listed against JAK2, JAK2 agreement only at class level, the CXCR2 genetic estimate
  against antagonism, and "strongest among this atlas's genes");
- a fourth reviewer on the final text found no error of biology or drug mapping and asked for rewordings, applied:
  the JAK argument from the atlas's own genetics (JAK1-selective drugs work although JAK1 has no genetic association
  here), the positive controls shown as where genetics recovers known targets rather than as a test, the antagonists
  also listed against CXCR1, one bioinformatic MR study not read as a verdict, and ustekinumab blocking IL-23 through
  IL12B in the applied specialist's finding. A fifth reviewer on the published text judged it acceptable with no
  change required; its one low point (tofacitinib's selectivity) and an optional one (a TYK2/JAK1 inhibitor,
  brepocitinib, effective in UC phase 2b, so deucravacitinib's failure is not read as TYK2 failing) were added from
  two more peer-reviewed papers. A sixth review asked only that tofacitinib not be described by its JAK2 activity on
  the JAK2 card; it is now described, as in the cited selectivity study, as inhibiting JAK1-dependent pathways most
  potently, like the others. Final reading: 523 claims, none contradicted, 14 not in the cited evidence, none in
  a statement about or naming IL12B, JAK2 or CXCR2.
The drug facts come from a new curated layer (drug_literature.py, data/drug_literature.json): verbatim sentences from
ten peer-reviewed trial papers, every one checked against PubMed, and the NLM MeSH record that ties navarixin to the
codes SCH-527123 and MK-7123. Smillie 2019's cell types are a new atlas record (51 types, no neutrophil type). A false
positive of the direction check found on the way ("the antagonist lowered blood neutrophil counts" read as CXCR2 going
down) was fixed: a clause where a drug acts before the direction word is not read as the gene's expression.

## (g6) Front-page rule, drugs only from papers, PubMed for the orchestrator, a cited expert review (2026-11-27)

Four changes, so that the review done by hand in (g4) and (g5) is part of every run:
- Front-page rule (`frontpage.py`). A gene may be prioritised when it meets a fixed criterion, and the card says which:
  a drug approved in UC on the gene or its class with an Open Targets genetic association of at least 0.8 (a positive
  control, always shown: IL12B and JAK2 in this data, which is why JAK2, rank 36, is on the page), one of the first 30
  of the convergence ranking, the leading first-in-class candidate by pooled increase or by retained share, or a
  significant Mendelian randomisation estimate. `verify.py final` fails a front page with a gene that meets none or
  without a mandatory control. All 14 genes of run 2026-11-23 met a criterion. The control block is now titled "a
  gene or drug class with a drug approved in UC" and the JAK2 card is subtitled "JAK inhibitor class".
- Drugs only from papers or datasets (`semantic.uncited_drugs`). A drug named in a statement must now appear in a
  cited peer-reviewed paper record or a treatment arm, not only in Open Targets or the clinical-indications table,
  which say a drug is listed for a target but not what it binds or how it did. Run 2026-11-23 fails this rule in 12
  statements (reparixin, ladarixin, SX-682, ABX-IL8, adakitug, atezolizumab, canakinumab, nangibotide, marimastat,
  doxycycline, SB-656933). The drug layer gained 11 peer-reviewed papers for these drugs (phase 3 of reparixin and
  atezolizumab, CANTOS for canakinumab, phase 2a and 2b of nangibotide, marimastat in pancreatic cancer, ladarixin in
  type 1 diabetes, ABX-IL8 and HuMax-IL8 trials, SX-682 preclinical, subantimicrobial doxycycline as an MMP inhibitor)
  and the MeSH record that makes adakitug HuMax-IL8 (BMS-986253): 23 papers, 40 quotes, all verbatim in PubMed.
- PubMed during a run (`pubmed_quotes.py`, kg-mcp `pubmed_search` and `pubmed_abstract`). The orchestrator may search
  PubMed for a fact no record holds and cite it only after registering a sentence: the paper must be a peer-reviewed
  journal article (not a preprint, comment, letter, editorial or erratum), not retracted, and the sentence verbatim in
  its abstract. Registered sentences are stored in `data/pubmed_quotes.json` and become pubmed_quote evidence.
- Expert review (`expert_review.py`, `.claude/agents/revisor-experto.md`, step 4c of /informe). A reviewer prompted
  as an IBD gastroenterologist and translational scientist reads the front page (summary and cards, as rendered) and
  writes objections: factual ones must give a PMID and a verbatim abstract sentence of a peer-reviewed paper, checked
  as above, and are discarded otherwise; framing ones may only reword, reorder or remove. The orchestrator applies
  them once and the corrected text goes through the reading check and verification.
Measured with seeded errors on run 2026-11-23: 12 errors of the kinds the manual reviews found (8 factual: the p40
and p19 subunits swapped, anti-p19 antibodies said to bind IL12B, upadacitinib called JAK2-selective, the failed
deucravacitinib trial said to succeed, CXCR2 antagonists said approved in UC, MK-7123 said to raise neutrophil counts,
nangibotide said approved for UC, "blocking interleukin-12 is required"; 4 framing: the JAK2 control presented as a
new JAK2-drug candidate, genetics said to support PD-L1 inhibition, an exploratory candidate called validated, a
neutrophil marker called a T-cell target) were written into copies of the cards, split over two blind batches so that
each reviewer saw a card either altered or unaltered, never both. The reviewers caught 12 of 12, 6 of the 8 factual
errors with a PubMed citation that passed the check and 2 as contradictions of the card's own text, and raised no
objection on the 28 unaltered items. The errors were written by the author of the reviewer's instructions, so this
measures the kinds of error already known, not unknown ones.
Run 2026-11-27, the whole /informe redone with these rules: the six specialists' first drafts lost 24 of 66
findings at verification, mostly to the drug rule (names cited only to Open Targets) and to pathway names put in the
gene field; after one revision round with the reading check's list, and a second for bare digits, all 67 findings
pass. Front page: the same 14 genes, each meeting the rule. The expert reviewer read the 18 front-page items and raised
3 factual objections, all with a verbatim quote from a peer-reviewed paper, all kept and applied: CD274 now says that
enterocolitis is among the leading adverse events of checkpoint inhibitors, so its drugs are a safety warning; MMP1
that the anti-MMP9 antibody andecaliximab did not induce remission or response in UC; S100A8 that with S100A9 it forms
calprotectin, an established stool test. Verification: 115 statements, no failures. Reading check: 509 claims, none
contradicted, 16 not in the cited evidence (none on the front page); on 73 direction-flipped copies the reader caught
70 and the deterministic checks 54 (29 of 31 where bulk expression is cited). Two of the missed flips sit on a
meta-analysis record marked not significant although its padj is tiny (JAK2), and on pathway findings without a gene.

## (g7) No ranking on the front page; a fifth list for targets already in the clinic (2026-11-27)

Which gene comes first depends on what the reader is after, so the front page no longer shows a ranked selection. It
keeps the summary and links to five lists, one per question: first-in-class, well known without a drug, repurposing,
genetics first, and a new one, targets in clinical validation (`clinical_list.py`). The orchestrator's selection, with
its selection rule and blocks, moved to the Ranking tab; the one summary sentence that cited ranking places was
removed (re-verified and re-read). The new list holds genes with a drug tried in UC or IBD that the omics back
(meta-analysis BH < 0.05, same direction in at least 5 of 8 cohorts, more than 100 Europe PMC UC records) and, at the
user's request, every gene with a drug approved in UC or IBD, marked approved, on the same scores: 77 genes, 22
approved. Approved targets fall low: median convergence rank 763 and median data-only rank 2,578 (expression,
replication, inflammation and genetics only) of 26,551 genes; TNF is 222 on the full score and 1,888 on data only,
ITGB7 3,611 and 7,539. Genes in trials lead the list (CXCR2, CXCR1, MMP1, MMP7, IL1B). Rows whose only UC drug is
doxycycline (MMP1, MMP7, MMP13) are flagged as not gene-directed.

## (g8) One tab for new targets: Discovery merged into the First-in-class shortlist (2026-11-28)

The Discovery tab and the shortlist answered the same question in two places, and the screen table (top 40 genes per
cell family by novelty score) read as a ranking of new targets although its first rows were genes the screen sets
aside (SERPINB5, DEFA5, ABCA12: well studied in UC) and it held only 18 of the 49 first-in-class genes. The two tabs
are now one, "First-in-class shortlist": the shortlist first (lights, groups, cards with the proposed experiment, the
papers, the calibration), then "How these 49 genes were found" (the screen's counts, where the set-aside genes are,
the novelty and safety charts, the caveats). The screen table and the Spotlight cards (a repeat of groups A and B) were
removed; the list of genes with drugs only for other diseases moved to the Repurposing tab as its full list before
scoring. No data changed; the report has 796 claims (132 fewer: 120 table rows and 12 cards).

## (g9) The orchestrator's decision: four targets to test first per list (2026-11-28)

The question changed to "Which new therapeutic targets for ulcerative colitis, not yet treated by any drug in UC, have
the strongest evidence to be tested first?". The orchestrator's selection (positive controls, the expected neutrophil
signal, new candidates and biomarkers mixed) is replaced by a decision with a fixed rule (`decision.py`): four genes
from group A of each list of targets with no therapy in UC, none with a drug tried or approved in UC or IBD. The pool
(32 genes: 14 first-in-class, 8 well known with no drug, 5 repurposing, 5 genetics first) and a dossier per gene were
given to the six specialists, who reported on the pool genes their tools cover (86 findings, all verified after one
round of fixes and one reading-check revision). The orchestrator chose SDR16C5, CYP4X1, CPXM1, CHRDL2; SERPINB5, CXCL1,
CHI3L1, REG1B; INHBA, CFB, LY96, SELP; TUFM, CEBPB, FIGNL1, CFAP90, and wrote why each was or was not chosen (XKR9 up in
four cohorts only; LILRA1, CEACAM4 and PCDH17 mostly cell arrival; MMP10 knockout mice worse; GUCY1A1 a Moyamoya risk;
RORC to be restored with no route and its regulon up). Positive controls (IL12B, JAK2) left the decision and stay as
reference in the clinical list. No biomarker is filtered out: genes with a biomarker record are marked "also a
biomarker" (CHI3L1, REG1B, INHBA, CFB, SERPINB5).

The expert reviewer gives a verdict per card. It was run with Claude Sonnet: the Opus 5.5 reviewer was stopped twice by
the model's safeguards (flagged as biology research). Verdicts: SERPINB5 support; 11 concerns; SELP, TUFM, FIGNL1 and
CFAP90 against, each kept with an answer on its card (SELP last of its list because the alternative carries a vascular
risk; TUFM because it holds the only causal genetic support; FIGNL1 and CFAP90 because the rule takes four and RORC has
data against its direction). Its three cited factual objections were applied (maspin overexpressed in IBD-associated
neoplasia; follistatin attenuated experimental colitis; an anti-P-selectin antibody did not improve experimental
colitis). Seeded measurement (12 new errors written for the decision cards): 10 caught, all with a valid objection,
24 objections on 28 unaltered items; missed: a reversed mouse result for SERPINB5 and "a transcription factor is easy
to drug". `docs/expert-review-redo.md` redoes the review on the saved input.

The data sections of the narrative were carried over from run 2026-11-27 and verified and read again; the summary, the
decision paragraphs and the 16 rationales are new. Final: 140 statements, 0 failures, decision rule 0 problems;
reading check 674 claims, 0 contradicted, 7 not in the cited evidence; reader caught 47 of 48 flips, deterministic
checks 30. The report moves the convergence ranking to Trust and methods ("Convergence score") and adds the decision as
the last Findings tab; the Agent benchmark tab now shows the decision rule and the expert review results.
