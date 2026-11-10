# Reproducing the ulcerative colitis atlas

## What you need

- [uv](https://docs.astral.sh/uv/) (Python 3.12 and every package are pinned in `uv.lock`).
- About 6 GB of disk for the raw data, 16 GB of RAM.
- For a new agent run only: Claude Code logged in with a subscription (no API key is used anywhere).

## 1. Data

```bash
uv sync --frozen
uv run python -m ucatlas.download          # bulk cohorts, platforms, HGNC, HPA, Open Targets tables, UC GWAS
uv run python -m ucatlas.mr download       # Open Targets credible sets filtered to UKB-PPP pQTL (~1.8 GB streamed)
uv run python -m ucatlas.spatial_download  # Xenium GSE312420 expression files
uv run python -c "from ucatlas.spatial_download import run_extra; run_extra()"   # Xenium GSE339323 + Visium GSE189184 (with H&E)
uv run python -m ucatlas.mr_colon download # Open Targets credible sets filtered to gut eQTL (GTEx colon, CEDAR, IBDverse)
```

Two single-cell datasets need a Broad Single Cell Portal login and are placed by hand (see `data/README.md`):
Smillie et al. 2019 (SCP259) in `data/raw/scp259/` and Gubatan et al. 2026 (SCP3755) in `data/raw/scp3755/`.
Every downloaded file is listed with its URL, date, size and SHA-256 in a manifest next to it.

## 2. Deterministic rebuild

```bash
uv run python -c "from ucatlas.reproducibility import verify; print(verify()['ok'])"   # check the frozen inputs
uv run ucatlas run-all                     # analyses, knowledge graph, verification of the saved agent run, report
uv run pytest -q                           # tests (data-free tests also run in GitHub Actions)
uv run python -m ucatlas.evaluate          # evaluation against UC clinical targets
```

`run-all` refuses to start when a frozen input, `config.yaml` or `uv.lock` differs from
`reproducibility/inputs.json`. It deletes the processed tables and the graph first and keeps only the single-cell
pseudobulk files, which come from long streaming steps. AlphaFold models and PDBe-KB ligand sites are fetched once
into `data/cache/alphafold` (not part of the frozen inputs). It does not call a language model: the agent outputs of the
published run are inputs, re-verified against the rebuilt graph.

A fresh download is not guaranteed to give the same bytes: GEO, NCBI and Open Targets update their files. The
manifests record what was used.

## 3. A new agent run

```bash
./scripts/run_informe.sh        # Windows: powershell -File scripts/run_informe.ps1
```

The orchestrator (Claude Opus) writes `runs/<date>/convergence.json`, launches six Sonnet specialists in parallel,
verifies their findings, writes the narrative, checks that every key finding reaches it, verifies the narrative and
renders `reports/uc_report.html`. Text from a language model is not bit-for-bit repeatable; the saved run is what the
report shows. Then refresh the frozen inventory with
`uv run python -c "from ucatlas.reproducibility import freeze; freeze('<run folder>')"`.

## 4. Publication

`.github/workflows/pages.yml` publishes `reports/uc_report.html` to GitHub Pages whenever it changes on `main`.
