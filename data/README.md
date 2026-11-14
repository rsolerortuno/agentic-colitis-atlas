# Data and their checks (M1)

Downloaded by `uv run python -m ucatlas.download` (manifest: `data/raw/manifest.json`, with URL, date, sha256 and
size). The sample table is built by `uv run python -m ucatlas.cohorts` → `data/processed/samples.parquet`.
Only colon and rectum biopsies of ulcerative colitis (UC) patients and non-IBD controls are used.

## Bulk cohorts, checked against GEO on 2026-11-09

| Cohort | Role | Platform | GEO samples | Used (baseline, colon) | Checked against the expectation |
|---|---|---|---|---|---|
| GSE193677 (Mount Sinai MSCCR, PMID 36109152) | Discovery: inflamed vs non-inflamed vs control | RNA-seq, GPL16791 (raw counts) | 2,490 biopsies, all regions, UC/CD/control | UC inflamed 293, UC non-inflamed 391, control 340; 135 UC patients have both an inflamed and a non-inflamed colon biopsy | As expected. Mayo endoscopic score, Nancy and GHAS histology, calprotectin and CRP are available |
| GSE206285 (UNIFI, ustekinumab; PMIDs 36192482, 36170836) | Replication; response | Affymetrix HT HG-U133+ PM, GPL13158 | 568 | UC 550 (ustekinumab 364, placebo 186), healthy 18 | As expected (550 + 18). Endpoints: mucosal healing and clinical remission at week 8 |
| GSE92415 (PURSUIT-SC, golimumab) | Replication; response | GPL13158 | 183 | UC 87 at week 0 (golimumab 59, placebo 28), healthy 21; 75 at week 6 | As expected (87 + 21). Endpoint: clinical response at week 6; Mayo score |
| GSE16879 (infliximab; PMID 19956723) | Replication; response | Affymetrix HG-U133 Plus 2, GPL570 | 133 (UC, CD, control; before and after) | UC 24 before infliximab (8 responders, 16 non-responders), colon control 6 (CD and the 6 ileal controls excluded) | Matches the paper (24 UC, 8 responders) |
| GSE73661 (infliximab and vedolizumab; PMID 27802155) | Replication; response | Affymetrix HuGene 1.0 ST, GPL6244 | 178 | UC 67 at week 0 (infliximab 23, vedolizumab 41, placebo 3), control 12 | Response is the label of each patient's earliest post-treatment biopsy (W4–6 infliximab, W6 or W12 vedolizumab): infliximab 8 R / 15 NR, vedolizumab 9 R / 32 NR. A third-party summary reported 15 R / 8 NR for infliximab; GEO titles say 8 R / 15 NR, and GEO is used |

| GSE109142 (PROTECT; PMID 30604764) | Replication (added 2026-11-10) | RNA-seq, NCBI-generated raw counts (GRCh38.p13) | 226 rectal biopsies | UC 206 (treatment-naive, paediatric), control 20 | As expected (206 + 20) |
| GSE117993 (RISK) | Replication (added 2026-11-10) | RNA-seq, NCBI-generated raw counts | 190 rectal biopsies | UC 43, control 55 (Crohn's disease 92 excluded) | As expected |
| GSE165512 (IL23-IL17 axis in IBD) | Replication (added 2026-11-10) | RNA-seq, NCBI-generated raw counts | 170 (colon and ileum) | UC colon 40, healthy colon 35 (31 in the count file); ileum and Crohn's excluded | Four control samples of the series are missing from the NCBI count file |

**The infliximab patients are the same people.** GSE16879 and GSE73661 come from the same Leuven group. Matching the
expression profiles of the two series pairs 22 of the 24 GSE16879 patients with an infliximab patient of GSE73661, and
its colon controls CO2-CO5 with four GSE73661 controls: the 23 infliximab patients and six controls of GSE73661
(biopsy numbers below 600) were measured again on another array. Each person is used once: GSE16879 for infliximab
response and for UC vs control (24 UC, 6 colon controls); GSE73661 for vedolizumab response and, in UC vs control, only
its 44 vedolizumab and placebo patients and its 6 later controls. Disease activity uses all 67 GSE73661 patients.

## Single cell (for M11)

SCP259 (Smillie et al. 2019), downloaded by the user into `data/raw/scp259/` (portal login required). Manifest with
sha256: `data/raw/scp259/manifest.json`. Checked: the matrix headers match the barcode and gene files
(epithelial 20,028 genes × 123,006 cells; stromal 19,076 × 31,872; immune 20,529 × 210,614), and `all.meta2.txt`
describes all 365,492 cells (the paper reports 366,650 before its own filtering). The metadata give the subject,
health state (healthy, non-inflamed, inflamed), compartment and cell cluster.

## Single cell with neutrophils (for M11-M12)

SCP3755 (Gubatan et al. 2026, Front Immunol, doi:10.3389/fimmu.2026.1705328), downloaded by the user into
`data/raw/scp3755/`: a Seurat object (.RDS, 1.9 GB) and its cell metadata (.csv). 408,930 cells from 20 UC patients
(inflamed and uninflamed biopsies, BD Rhapsody; no healthy controls), 47,516 neutrophils in five states. No R is
installed: the RDS is parsed in Python (`rdata`); the RNA count matrix is the only dgCMatrix with one column per cell
(32,365 genes, integer counts). Second dataset for lineage replication: Garrido-Trigo et al. 2023 (CZ CELLxGENE,
46,700 cells, 6 UC / 6 CD / 6 healthy donors, five coarse lineages), downloaded with sha256 in the manifest.

## Disease identifiers

MONDO:0005101 (cross-referenced to EFO_0000729 and DOID:8577). Open Targets (2026) uses MONDO_0005101; a query
with EFO_0000729 returns no disease, checked 2026-11-09.

## Mendelian randomisation (added 2026-11-10)

- `data/raw/gwas/GCST004133_uc_delange2017.h.tsv.gz`: UC GWAS of de Lange et al. 2017 (Nat Genet), GWAS Catalog
  harmonised file (GRCh38), 300 MB.
- `data/raw/extra/ot_study.parquet`: Open Targets 26.09 study index (to map UKB-PPP pQTL studies to genes).
- `data/raw/extra/ot_pqtl_credible_sets.parquet`: extract of the Open Targets 26.09 credible sets (200 parts, about
  1.8 GB streamed and filtered by `uv run python -m ucatlas.mr download`) keeping the UKB-PPP European pQTL credible
  sets, with each 95 % credible set's variants and betas. The manifest records the source and the filter.

## Spatial (Xenium and Visium)

- `data/raw/spatial/GSE312420/`: Xenium 5K and 337-gene panel expression files (M11).
- `data/raw/spatial/GSE339323/` (added 2026-11-10): Xenium Prime 5K, one UC resection, most and least inflamed region;
  cell x gene h5 and cell table only (the archive is 25 GB, mostly images and transcripts).
- `data/raw/spatial/GSE189184/` (added 2026-11-10): 10x Visium, 7 inflamed UC and 2 healthy colon sections; only the
  filtered spot x gene matrix, the spot positions and the scale factors are read out of each sample zip with HTTP range
  requests (`uv run python -c "from ucatlas.spatial_download import run_extra; run_extra()"`).

## Added on 2026-11-10 (third and fourth rounds)

| Data | Source | What is used | Note |
|---|---|---|---|
| ACT1 infliximab and placebo (GSE23597, Toedter 2011) | GEO series matrix, GPL570 | baseline biopsies, clinical response at week 8 (infliximab 25 R / 7 NR, placebo 5 R / 8 NR) | GSE12251 (Arijs 2009) was checked and dropped: its 22 patients are in GSE23597 |
| EUCALYPTUS etrolizumab (GSE72819) | NCBI-generated raw counts | 70 baseline biopsies with remission at week 10 (12 R / 58 NR) | treatment-response cohort only (no controls) |
| PROTECT outcomes (GSE109142) | series matrix | initial treatment (corticosteroids or 5-ASA), remission at week 4, baseline faecal calprotectin (145 children) | |
| Mucosal proteomics (Schniers et al. 2019) | Clin Proteomics 16:4, additional file 2 | UC/healthy ratios and p of 8,562 protein groups | patients independent of every RNA cohort |
| FinnGen R12 UC (strict) | FinnGen public summary statistics | second, independent MR outcome | 7,220 cases |
| Liu et al. 2023 UC (GCST90446794) | GWAS Catalog harmonised | third MR outcome, largest (23,252 cases, European + East Asian) | includes de Lange 2017 and FinnGen samples: reported, not combined |
| CosMx (GSE312415) | GEO tar archive (6.5 GB, streamed) | expression matrices and cell metadata of 5 slides; 2 used (both biopsies inflamed UC) | the other slides mix biopsy states without a cell-to-biopsy key |
| Anti-TNF single-cell atlas (GSE282122, Thomas 2024) | authors' supplementary tables (Europe PMC) | per-cell-state remission and longitudinal models for 20 UC patients | GEO matrices carry no cell annotation |
| AlphaGenome Atlas and AlphaGenome 2 (Google DeepMind) | API (free key in ALPHAGENOME_API_KEY, not stored); every answer cached in raw/alphagenome with a sha256 manifest | predicted RNA effect of UC GWAS credible-set variants in colon tissue, blood immune cells and 12 colon cell types | outputs for non-commercial research only; not under CC BY |
| UC and IBD GWAS credible sets | Open Targets 26.09 credible_set dataset, filtered to eight studies (raw/extra/ot_uc_ibd_gwas_credible_sets.parquet) | fine-mapped risk signals for the genetics-first list | SuSiE-inf, SuSiE and PICS |
| Lyu et al. 2025 UC meta-analysis (medRxiv 10.1101/2025.09.23.25336483) | authors' GitHub repository, pinned to a commit (MIT licence), 22 files in raw/gwas/lyu2025_uc | European meta-analysis of MVP, FinnGen R12, Pan-UK Biobank and IIBDGC, GRCh38 | preprint, not peer reviewed; alleles in MarkerName are alphabetical, REF and ALT columns give the reference |
| Mouse models | MGI phenotype reports; Europe PMC search | colitis phenotypes of single-gene models; mouse colitis papers per candidate | |
| Constraint and safety | gnomAD v4.1 constraint; Open Targets safety liabilities (API, cached) | LOEUF, curated safety events | |
