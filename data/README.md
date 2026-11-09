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
| GSE16879 (infliximab; PMID 19956723) | Replication; response | Affymetrix HG-U133 Plus 2, GPL570 | 133 (UC, CD, control; before and after) | UC 24 before infliximab (8 responders, 16 non-responders), control 12 (CD excluded) | Matches the paper (24 UC, 8 responders) |
| GSE73661 (infliximab and vedolizumab; PMID 27802155) | Replication; response | Affymetrix HuGene 1.0 ST, GPL6244 | 178 | UC 67 at week 0 (infliximab 23, vedolizumab 41, placebo 3), control 12 | Response is the label of each patient's earliest post-treatment biopsy (W4–6 infliximab, W6 or W12 vedolizumab): infliximab 8 R / 15 NR, vedolizumab 9 R / 32 NR. A third-party summary reported 15 R / 8 NR for infliximab; GEO titles say 8 R / 15 NR, and GEO is used |

**Caveat: infliximab patients probably overlap.** GSE16879 and GSE73661 come from the same Leuven group, with
nearly the same infliximab split (8 R / 16 NR and 8 R / 15 NR). They are not treated as independent replication
cohorts for infliximab response: GSE16879 is used for infliximab, GSE73661 for vedolizumab.

## Single cell (for M11)

SCP259 (Smillie et al. 2019), downloaded by the user into `data/raw/scp259/` (portal login required). Manifest with
sha256: `data/raw/scp259/manifest.json`. Checked: the matrix headers match the barcode and gene files
(epithelial 20,028 genes × 123,006 cells; stromal 19,076 × 31,872; immune 20,529 × 210,614), and `all.meta2.txt`
describes all 365,492 cells (the paper reports 366,650 before its own filtering). The metadata give the subject,
health state (healthy, non-inflamed, inflamed), compartment and cell cluster.

## Disease identifiers

MONDO:0005101 (cross-referenced to EFO_0000729 and DOID:8577). The Open Targets id is confirmed in M4 when the
association query is first run.
