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
- **Groups and cell-system experiments (2026-11-11)**: a check of the shortlist order under six scorings showed that
  only nine genes are in the top ten under all of them and that RFTN1 swings between 2nd and 23rd; the table now reports
  groups and stability instead of positions. The experiment selection was found half written on disk (a truncated
  string, the module did not import) when the session was forked; it was completed with nine cell systems, a second
  system from 20 % of the RNA and the published functional lead.
- **Crypt-niche flag and fourth atlas (2026-11-11)**: a stromal counterpart of the crypt regeneration flag, calibrated
  on WNT2B, GREM1 and RSPO3 against COL1A2; no shortlist gene is flagged. GSE282122 full matrices downloaded and typed
  into lineages. Two problems on the way: the GEO sample titles carry a "-reup" suffix that the tar folders lack (mapped
  explicitly), and the series mixes two Cell Ranger references (33,538 and 33,808 features), which stopped the first
  build; genes are now matched by Ensembl id. The four-atlas replication is written to a new file and new evidence
  types, so the sc_summary records cited by the agents stay as they were (0 stale). OLFML2B leaves group A.
- **Agent run with four atlases (2026-11-11)**: runs/2026-11-11 replaces runs/2026-11-10 as the published run. The
  sc_summary records now carry the lineage in all four atlases and replicated = two of three, so the cell tool shows the
  agents the new replication; the old run would read as stale against them, which is why the agents were re-run rather
  than patched. First verification dropped 33 of 73 agent statements, all for form (numbers printed next to their
  annotation, unbound digits such as "CXCR1/2", "Phase 2-3", PMIDs; pathway findings carrying a gene field); the agents
  fixed their files from the verifier output and four wording fixes were made by the orchestrator. Final: 116
  statements, 0 unsupported, 0 stale, 28 of 28 key findings cited. Report: the old second-atlas block is merged into a
  four-atlas block with a new figure, the cell count, help texts, footer and data appendix cover all four atlases, and
  the anti-TNF section adds a descriptive remission check on the count matrices.
- **Clean rebuild check (2026-11-12)**: `run-all` from frozen inputs, twice. The first attempt stopped because an
  idle MCP server of another session held the graph file open, after the processed data had been deleted; run-all now
  clears the graph first and empties its tables when the file is locked. The first full rebuild then exposed an order
  bug: the four-atlas replication read the discovery gene list from discovery_scored.parquet, which run-all writes
  later, so it covered 390 genes instead of 1,003 (OLFML2B fell back to the two-atlas rule). Fixed to read
  discovery.parquet; the second rebuild reproduces the published results (1,003 genes, OLFML2B in group B, 0
  unsupported, 0 stale, report valid). The mouse-colitis literature table now follows the targets of runs/2026-11-11.
- **Neutrophil-made genes (2026-11-12)**: a review of group A on request showed ADGRE3 scored on macrophages
  (Smillie has no neutrophils). First fix capped the lights at amber as "not measured"; asked whether that was unfair,
  the per-cell change was then measured inside the Gubatan neutrophil states, which showed no rise for ADGRE3 (red), so
  the amber had favoured it. Group A grew to 13 through ties at the top-ten border.
- **Drugs by indication (2026-11-12)**: on request, a column saying whether a gene's drugs are for UC, and a
  repurposing table for discovery genes drugged only for other diseases. Two problems found while building it: Open
  Targets lists epacadostat as approved from its ATC code alone (now "not confirmed"), and many cancer drugs act by
  activating the target or by carrying a toxin to cells, which does not fit genes that rise in UC; the table now says
  whether an inhibitor exists.
- **Repurposing shortlist (2026-11-12)**: kept apart from first-in-class (different question, and only GUCY1A1 of
  the 35 would pass the first-in-class paper filter). First draft put CTLA4 in group A; immune checkpoints are now red.
  A stray control character (a backspace written by a shell escape in a patch script) broke the salt-name regex; found
  by checking the output, removed, and changed files were scanned for control characters.
- **Repurposing tab (2026-11-12)**: on request the repurposing list moved from the shortlist tab to its own tab
  (Findings, New targets: Discovery, Validation shortlist, Repurposing); the shortlist tab is back to its earlier
  content. The proposed test moved from a table column to cards, which had made the rows very wide; table headers of
  the shortlist-style tables may now wrap, so both tables fit the page (measured in the browser: 975 px in a 975 px
  box, rows 58 px), and the long drug texts of the Discovery leads table share one bounded column.
- **Repurposing literature (2026-11-12)**: shortlist_literature.py generalised to two gene sets with separate
  folders and outputs (first-in-class output unchanged: 24 grey, 10 amber, 1 red, 1 green). Classification by 11
  subagents in parallel took far longer than expected (about 1 h 50 min for three batches each, likely throttling);
  two agents were cut off by a session restart and their three batches were rerun. 816 of 818 calls passed the
  verbatim check.
- **Hidden regex bug (2026-11-12)**: the colitis-setting test of the literature review (shortlist_literature.IBD)
  had backspace characters where word boundaries were meant, so "IBD" and "UC" as abbreviations never matched. Fixed;
  no literature light changes in either list (papers using the abbreviations also name colitis or Crohn's disease).
- **Tops read for an expert (2026-11-12)**: composition adjustment dropped for neutrophil-made genes and replaced by
  a weighting by where the RNA is (the per-neutrophil change alone moved CLEC1A to group A, rejected); Discovery
  spotlight follows the shortlist groups; Top 25 gets a "How to read it" column. The Top 25 table did not fit the page
  (1,777 px in 975); the duplicate drug column was merged, headers wrap and three citations are shown per row.
- **Three lists (2026-11-12)**: "Well known in UC, no drug yet" added before the first-in-class shortlist and the
  repurposing list (tab order as requested); neutrophil-made genes no longer removed by the discovery tissue-change
  filter; literature reviewed for 37 + 17 genes with subagents. Europe PMC retrieval sped up (parallel full texts); the
  first attempt with sequential full texts ran for 28 minutes without finishing.
- **Defence flag (2026-11-12)**: first applied as a safety flag only; DUOXA2 then stayed second while DUOX2 fell,
  so the flag also caps the tissue light, as the crypt regeneration flag caps the cell light. My earlier remark that
  SERPINB5 was an odd first gene was wrong for current evidence (two 2026 colitis knock-down studies).
- **AlphaGenome (2026-11-12)**: new layer alphagenome_uc.py on request. Validated before use against fine-mapped
  eQTL (GTEx colon for the Atlas, IBDverse cell types for AlphaGenome 2); weak predictions are near chance, so only the
  top quarter of |score| counts. Atlas ontology filters reject composite cell-type terms and gene filters need
  versioned IDs, so tracks are grouped by biosample name and genes matched after the query. Result: CCL20 and IL10
  are the only genes of interest with a confident effect; nothing enters the target lists' scores.
- **Genetics-first list (2026-11-13)**: requested as a fourth list from AlphaGenome over all UC risk signals. Two
  errors found and fixed on the way. (1) The Lyu 2025 files list alleles alphabetically, not as reference and
  alternative; the Atlas rejected the variants and half of the Open Targets variants did not match. (2) My first
  rule kept every gene whose predicted change was large enough for the direction of a known eQTL to be right, and I
  had lowered the cell-type threshold on that basis. It gave 316 genes that agreed with Mendelian randomisation at
  chance (3 of 7 confident), because the check measured direction, not which gene: 2.4 % of the other genes near an
  eQTL variant pass. The rule now requires the top gene at its variant, twice the next one; 6 genes remain. The
  CCL20 and IL10 calls published on 2026-11-12 used the first rule and are withdrawn in the Genetics tab.
- **Review of the whole atlas (2026-11-13)**: requested to find what was still wrong. Found and fixed:
  (1) GSE16879 compared UC colon with twelve controls of which six are ileal biopsies; only the six colon controls are
  used now. (2) The infliximab patients and six controls of GSE73661 are the people of GSE16879 measured on another
  array (matching expression profiles pairs 22 of 24 patients and controls CO2-CO5; biopsy numbers below 600), so the
  meta-analysis counted them twice. GSE73661 now enters UC vs control with its 44 vedolizumab and placebo patients and
  its 6 later controls; response and disease activity are unchanged. Effect: leading genes barely move (CXCR2 2.81 to
  2.84, S100A8 4.77 to 4.79), two-fold changes 1,625 to 1,657, screening pool 1,111 to 1,127, first-in-class 49 with
  six genes replaced at the margin. (3) sc_celltype_de.parquet was a kept file built before the symbol mapping of the
  pseudobulk and for an older gene list: 2,695 listed genes, CXCL8 among them, had no inflamed vs healthy test per
  cluster. It is rebuilt on every run (vectorised, 12 s instead of over 20 minutes; same values on 306 comparisons).
  (4) The published report printed "cost $0.00, 0.0 minutes" because the 2026-11-11 run had no session.json.
  (5) Typed numbers that had drifted: the literature-matched AUROC range, the pocket-score control ranks (4th to 25th
  of 45 proteins; now 7th to 45th of 83) and about twenty README figures. Measured and documented, not changed: the
  Visium, Xenium and CosMx neighbourhood tests draw random spots or cells as their null, so their p-values are too
  small (51 % of Visium pairs and 40 % of Xenium pairs pass); the report and README now say to read the enrichment
  size. Because every meta-analysis value moved, the agent run was redone (runs/2026-11-13: 122 statements, 0
  unsupported, 29 of 29 key findings cited) after the literature of the eight genes that entered the lists was
  retrieved and classified. Sabotage on the new run: numbers 77/77, gene 82/82, evidence 82/82, direction 22/28 with
  expression evidence. Checked without finding errors: the sample tables of the other cohorts, the differential
  expression and meta-analysis code, claim directions in the graph, evidence descriptions for all 54 record types,
  composition adjustment, discovery screen, profiles, treatment response, disease activity, the adjusted model,
  proteomics and the Mendelian-randomisation estimators. Not read line by line: validation_shortlist.py, the
  literature modules and the spatial viewers.
