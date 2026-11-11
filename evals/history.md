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

## 2026-11-10, closing four documented limits (decided before measuring unless marked)
- **Meta-analysis, single-cohort genes.** About 4,900 genes were measured only in MSCCR (k = 1) and 2,522 of them
  were "significant" on that one cohort. Two changes: (1) three RNA-seq cohorts were added as replication cohorts,
  UC vs healthy controls, from NCBI-generated raw counts (one pipeline for every series): PROTECT GSE109142 (206 UC,
  20 controls, paediatric rectum at diagnosis), RISK GSE117993 (43 UC, 55 controls, paediatric rectum; Crohn's left
  out) and GSE165512 (40 UC, 35 healthy adult colon; ileum and Crohn's left out). Every HGNC locus group is mapped, as
  for the MSCCR Ensembl ids. (2) A pooled estimate needs at least two cohorts: genes with k = 1 keep their single
  estimate but get no pooled adjusted p-value (BH over k >= 2 genes only), so they cannot enter the expression stream
  on one cohort. Remaining k = 1 genes are mostly pseudogenes, small RNAs and lowly expressed genes. The replication
  stream is a share of cohorts, so its scale is unchanged. Weights unchanged.
- **Mendelian randomisation with plasma pQTL (new layer).** Exposure: UKB-PPP (European) fine-mapped cis-pQTL credible
  sets from Open Targets 26.09; instrument = lead variant of each cis credible set with p < 5e-8, outside the MHC,
  not palindromic; when the lead is absent from the GWAS, the variant of the same 95 % credible set with the highest
  posterior probability that the GWAS has, with its own pQTL beta. Outcome: de Lange 2017 UC GWAS (GCST004133,
  harmonised GRCh38). Wald ratio or fixed-effect IVW; one Olink assay per protein (most instruments, then strongest);
  BH over proteins; colocalisation (coloc H4 >= 0.8 or CLPP >= 0.01 with any UC/IBD/Crohn's GWAS credible set in
  Open Targets) fetched for q < 0.05 proteins only. MR is reported and given to the agents; it is NOT added to the
  convergence score or the profiles, whose weights were fixed earlier.
- **Spatial.** Two GSE312420 Xenium sections had too few transcripts per cell to type most cells. Added (a) the
  Xenium Prime 5K sections of GSE339323 (one UC resection, most and least inflamed region) and (b) whole-transcriptome
  Visium sections of GSE189184 (7 inflamed UC, 2 healthy colon; checkpoint-inhibitor colitis left out) for genes
  that are not on the 5K panel. Visium compartment = lineage whose marker score has the highest median Spearman rho
  with the gene across spots of the UC sections (markers never include the gene; a section counts only if >= 20
  spots express the gene). The 50 %-typed rule for Xenium sections is unchanged.
- **Temporal subset (changed after seeing the 4-gene result).** The strict subset (first UC trial >= 2018, >= 2
  trials or approval) is kept for continuity. A wider set is added: targets of a drug with <= 5 targets that reached
  phase 2+ in UC (one trial is enough) and whose first UC trial of any phase started on or after 2016-01-01 with
  complete dates; earlier UC trial targets are removed from the universe. Decided after the strict subset proved too
  small, so it is logged here; the cutoff was not tuned on the result.

## 2026-11-10, new published run
- After the rebuild with eight cohorts, MR and the new spatial layers, the saved run of 2026-11-09 cited evidence that
  had changed (the verifier flagged it as stale, as designed). A new team run (runs/2026-11-10) is the published one:
  117 statements, 0 unsupported, 0 stale, 29 of 29 key findings cited, $4.08, 6.6 minutes. A first attempt the same
  day ran without the two MR tools (they were registered after `mcp.run()` in the server file) and was discarded.
- The benchmark runs of 2026-11-09 (three team runs, one single agent) stay as they were, on the earlier data.

## 2026-11-10, second round: gut genetics, spot deconvolution and the kidney-atlas views
Designs and weights below were written into the code before the step was run; none was changed after seeing results.
- **Gut eQTL Mendelian randomisation (`mr_colon.py`).** Same instruments, outcome and estimator as the plasma analysis,
  with Open Targets 26.09 credible sets of GTEx v10 gene-level eQTL (transverse and sigmoid colon), CEDAR eQTL
  (transverse colon, rectum) and IBDverse single-cell eQTL. eQTL credible sets carry no z-score, so |z| is derived from
  the p-value (mantissa and exponent). Colon layer: per gene the study with the most instruments; cell layer: the
  strongest cell type with its p multiplied by the number of cell types tested. BH within each layer.
- **Visium deconvolution (`spatial_deconv.py`).** Kidney-atlas design: SCP3755 profiles folded into 12 groups, up to
  40 markers per group (log2 ratio > 1, >= 20 CPM), candidates and the discovery pool excluded from the markers, NNLS
  per spot, NNLS of gene CPM on shares for the in-situ expresser, ring-1 neighbourhood with 1,000 section-stratified
  permutations (seed 20261110), BH over genes x groups.
- **Protein structure (`structure.py`).** Kidney-atlas pipeline without the molecular-glue components (that question
  is not asked here): pocket score = 0.5 confidence + 0.4 geometry + 0.1 hydrophobicity; controls are UC
  small-molecule targets. A PDBe 504 is retried rather than cached.
- **Report layers (`report_layers.py`).** Subgroups (paediatric vs adult cohorts, rectum vs proximal colon, neutrophil
  activity axis), biologic suitability (surface 0.30, increase 0.25, vital organs 0.25, patients 0.20), immune-target
  score (immune share 0.3, coupling 0.25, increase 0.2, tractability 0.15, response 0.1), a novelty score that only
  orders the discovery table, and rank stability of the convergence score. Known UC drug targets rank low on the
  biologic and immune scores (none of the immune controls in the top 100); this is reported, not tuned away.
- **New published run.** The new evidence and tools changed what the agents can cite, so the team ran again
  (runs/2026-11-10, replacing the morning run): 126 statements, 0 unsupported, 0 stale, 31 of 31 key findings cited,
  $5.57, 6.4 minutes. The orchestrator returned failing findings to the specialists once (numbers printed twice) and
  sent the cells specialist back for the in-situ findings of the new spatial section.
- **Report.** Rebuilt on the kidney-atlas layout: 16 tabs, a "?" help on every chart and counter, tissue viewer on the
  H&E sections, 3D protein viewer, candidate evidence cards, profile tables, rank-stability and outcome charts, agent
  benchmark tab.
- **Fingerprints ignore float noise.** The clean rebuild flagged 4 discovery records cited by the new run as stale.
  They differed only in the last binary digits of the composition-retained fraction (at most 174 units in the last
  place; a least-squares fit repeated in another process). Evidence fingerprints now round floats to 10 significant
  digits; the basis of runs/2026-11-10 was migrated (2,525 records unchanged, 4 float-noise only, 0 left stale; listed
  in its basis.json). A real change of a value is still caught (tests/test_limits.py).

## 2026-11-10, third round: independent checks
Designs fixed before running; none changes the convergence score or the profiles (they are reported as checks).
- **Proteomics** (Schniers 2019, additional file 2): authors' ratios and p as published, BH over proteins in > 70 %
  of samples; replication = same sign and padj < 0.05 for genes with meta padj < 0.05 and |log2FC| >= 1.
- **FinnGen R12 UC (strict)** as a second outcome for every MR estimate (same instruments; FinnGen beta for ALT);
  fixed-effect combination with de Lange 2017; replication = same sign and p < 0.05. Colocalisation API answers are
  now cached and retried (an empty API answer stopped one run).
- **Treatment cohorts**: EUCALYPTUS etrolizumab (GSE72819, remission week 10) and ACT1 infliximab and placebo
  (GSE23597, response week 8, baseline biopsies). GSE12251 was downloaded and dropped: its 22 patients are in GSE23597.
- **CosMx (GSE312415)**: slides hold two biopsies without a cell-to-biopsy key, so only the two slides whose biopsies
  are both inflamed UC are used. The Xenium marker-score typing gave lineages of equal size (background counts); a
  z-score rule typed 34 % with implausible proportions; final typing is by correlation with the single-cell reference
  profiles (r >= 0.1, margin 0.02). Both earlier rules are reported here because they were tried. Detection above
  background (best lineage >= 2x the median other lineage in % expressing cells) is the reported result.
- **Mouse models (MGI)**: single-gene genotypes, first with one marker (205 genes); then counting only markers with a
  human ortholog so that conditional knockouts (Cre transgene as second marker) count (236 genes; STAT3 enters).
- **Literature**: ten PubTator3 hits per report candidate, convergence top-40 and first-in-class gene (frozen in
  query_log.json, extra_gene_pmids); PubTator3's own relations summarised per gene. The 100-abstract cap of the
  checked extraction is unchanged.
- **Weight stability** of the biologic, immune, novelty and pocket scores; **thread pinning** (one BLAS thread) after
  the float-noise difference between two rebuilds; **licence** MIT + CC BY 4.0.
- **Report**: the published site is split (scripts/split_site.py): libraries cached, 3D structures and tissue images
  fetched when their tab opens; the single file stays self-contained for offline use.
- **Published run redone** with the new evidence and tools (protein_level, mouse_phenotype, pubtator_relations, FinnGen
  fields, CosMx): 125 statements, 0 unsupported, 0 stale, 31 of 31 key findings cited, $5.18, 8.0 minutes. The
  orchestrator changed CD274 from target to biomarker (genetically higher PD-L1 protects) and cites the protein-level
  replication of HMGCS2. The clean rebuild before it flagged 30 sentences of the previous run as stale (the treatment
  records gained the ACT1 and EUCALYPTUS arms), as designed.

## 2026-11-10, fourth round: closing checks
Designs fixed before running. Nothing below changes the convergence score or the profiles.
- **Covariates** (expression_adjusted.py): biopsy-level OLS with colon segment, age (per 10 years) and sex, patient-
  clustered SE; a change holds with the same sign, >= half the size and padj < 0.05.
- **MR sensitivity** (mr_sensitivity.py): Cochran Q, MR-Egger and weighted median (>= 3 instruments), Steiger with
  r2 approximated from z and N. A first run reported wrong counts because Cochran's Q overwrote the q-value column in
  the record; fixed before any result was used.
- **Liu 2023** (GCST90446794) as third MR outcome, reported apart because its samples overlap the other two.
- **Safety** (safety.py): gnomAD v4.1 LOEUF < 0.35; Open Targets safety liabilities; PheWAS from the cached
  colocalisation answers (H4 >= 0.8). The first PheWAS listed protein-level GWAS; molecular traits are now excluded.
- **PROTECT** outcomes: initial treatment and week-4 remission as two new arms; faecal calprotectin correlation.
- **Negative control** (report_layers.negative_control): expression- and literature-matched AUROC, 200 random sets.
- **Mouse literature** (mouse_models.literature_check): Europe PMC, title and abstract only (a first version searched
  full text and inflated counts, e.g. 11 papers for CHRDL2; replaced before use).
- **Anti-TNF atlas** (antitnf_atlas.py, GSE282122): the authors' per-cell-state tables, because GEO matrices carry no
  annotation; gene names mangled into dates by Excel in the published table are dropped.

## 2026-11-10, final benchmark (one agent vs the team) on the closing data
- Runs: 2026-11-10 (published), 2026-11-10_rep1, 2026-11-10_rep2 (team) and 2026-11-10_single. Same 12-fact
  checklist as the first benchmark (fixed on 2026-11-09). Two runs first stopped at the usage limit (API error, no
  tokens) and were redone. The single agent wrote a helper script into its run folder (to turn placeholders into
  numeric annotations); it is not an output and was removed.
- Team: 119 / 116 / 124 statements, 0 unsupported, checklist 10 / 11 / 10, cost $4.61 / $4.76 / $4.99, 7.5 / 7.4 /
  8.0 minutes. One agent: 141 statements, 0 unsupported, checklist 10, $7.52, 19.3 minutes, 197 turns. Target lists:
  team 0.87-1.0 overlap, single vs team 0.56. Conclusion unchanged from the 2026-11-09 benchmark (team: same quality,
  cheaper, faster), with the team now citing more of its evidence than the single agent (90-98 % vs 82 %).

- **Published run (final)**: runs/2026-11-10 redone on the closing graph: 119 statements, 0 unsupported, 0 stale,
  24 of 24 key findings cited, $4.61, 7.5 minutes. Targets: CXCR2, IL12B, JAK2, IDO1, CXCR1, MMP1 (targets or both);
  DUOX2, CHI3L1, CXCL8, HMGCS2, CD274, S100A8 (biomarkers); CHRDL2 and XKR9 (exploratory).

- **Validation shortlist tab (2026-11-11)**: new tab under Findings → New targets. First version scored the tissue
  criterion on the retained share with its own thresholds and the cell criterion on a different single-cell atlas
  than the discovery screen; both disagreed with the screen (19 reds on genes the screen had passed; EHD3 labelled
  stromal but "B cell" as source). Fixed before publishing by reusing the screen's composition verdict and the
  Smillie 2019 lineage shares.
- **Shortlist revision (2026-11-11)**: CCNO came first in the first version. Reviewed gene by gene and fixed four
  criteria (per-cell change, crypt regeneration flag, Mendelian disease, antibody route from topology; null genetics
  scores 0), each checked on markers and on 16 clinical targets. Adding the Mendelian field to the target-safety
  records made 3 agent statements stale (CXCR2, CD274 and one more); moved it to its own evidence record instead and
  the run is back to 0 stale.
- **Literature review of the shortlist (2026-11-11)**: new layer shortlist_literature.py and classifier agent. First
  retrieval run was slow (shared request helper retried 404s of closed full texts for a minute each); replaced with a
  short-retry cached GET that caches 404s, and four genes in parallel. Classifiers flagged their own judgement calls
  (tumour outcomes, glioma, aliases of other genes); papers outside a colitis or IBD setting were then excluded by a
  text check rather than by hand, and list-only mentions were made grey instead of amber so they earn no points.
- **Literature table (2026-11-11)**: the per-gene conclusions of the review were only in tooltips and in the ten
  cards; added a table to the shortlist tab with every candidate, its light and a one-line summary per specific paper
  (PubMed links), with what each kind of evidence can and cannot say.

