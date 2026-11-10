"""M1: download the bulk cohorts, HGNC and platform annotations, with a manifest (URL, UTC date, sha256, bytes).

Files are fetched once; a file already present is kept and re-hashed. The manifest is the record of what was used.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import requests

from .paths import RAW

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"


def _series(gse: str, *suppl: str) -> list[tuple[str, str]]:
    stub = gse[:-3] + "nnn"
    out = [(f"{GEO}/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz", f"{gse}/{gse}_series_matrix.txt.gz")]
    out += [(f"{GEO}/series/{stub}/{gse}/suppl/{f}", f"{gse}/{f}") for f in suppl]
    return out


def _platform(gpl: str) -> tuple[str, str]:
    stub = gpl[:-3] + "nnn" if len(gpl) > 6 else "GPLnnn"
    return f"{GEO}/platforms/{stub}/{gpl}/annot/{gpl}.annot.gz", f"platforms/{gpl}.annot.gz"


def _ncbi_counts(gse: str) -> tuple[str, str]:
    f = f"{gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    return f"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={gse}&format=file&file={f}", f"{gse}/{f}"


OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"

SOURCES: list[tuple[str, str]] = [
    *_series("GSE193677", "GSE193677_MSCCR_Biopsy_counts.txt.gz"),           # Mount Sinai MSCCR, RNA-seq (discovery)
    *_series("GSE206285", "GSE206285_additional_sample_metadata.tsv.gz"),     # UNIFI, ustekinumab
    *_series("GSE92415"),                                                     # PURSUIT-SC, golimumab
    *_series("GSE16879"),                                                     # infliximab, UC + CD
    *_series("GSE73661"),                                                     # infliximab + vedolizumab, UC
    # RNA-seq replication cohorts (added 2026-11-10 so that genes not on the arrays meet more than one cohort);
    # NCBI-generated raw counts (GRCh38.p13, one pipeline for every series)
    *[x for g in ("GSE109142", "GSE117993", "GSE165512") for x in (_series(g) + [_ncbi_counts(g)])],
    _platform("GPL13158"), _platform("GPL570"), _platform("GPL6244"),
    ("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt", "hgnc_complete_set.txt"),
    # second single-cell dataset for replication of the cell attribution (Garrido-Trigo et al. 2023, CZ CELLxGENE)
    ("https://datasets.cellxgene.cziscience.com/5e75dd8f-5c5f-4056-ade2-117748429022.h5ad", "sc_garrido2023/ibd_hc_10x.h5ad"),
    # Human Protein Atlas RNA tissue consensus (vital-organ expression for the discovery safety filter, M13)
    ("https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip", "extra/hpa_rna_tissue_consensus.tsv.zip"),
    # Open Targets 26.09 tractability for every target (two parts): predicted SM / AB tractability for all genes
    *[("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/target_tractability/"
       f"part-0000{i}-e2201205-916c-4b28-971d-6511ae9ea691-c000.zstd.parquet", f"extra/ot_target_tractability_{i}.parquet")
      for i in (0, 1)],
    # Open Targets 26.09 clinical precedence (all diseases, one part): gold set of the M8 evaluation
    ("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/evidence_clinical_precedence/"
     "part-00000-c18648ed-1b6a-4490-aecf-bfadcf459516-c000.zstd.parquet", "extra/ot_clinical_precedence.parquet"),
    # Mendelian randomisation (added 2026-11-10): UC GWAS summary statistics (de Lange et al. 2017, GWAS Catalog
    # harmonised, GRCh38) and the Open Targets 26.09 study index (to find the UKB-PPP plasma pQTL studies)
    ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/"
     "28067908-GCST004133-EFO_0000729.h.tsv.gz", "gwas/GCST004133_uc_delange2017.h.tsv.gz"),
    (f"{OT}/study/part-00000-90d9a8f0-e00a-492e-b961-40b56a21dbd0-c000.zstd.parquet", "extra/ot_study.parquet"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            shutil.copyfileobj(r.raw, f, length=1 << 20)
    if dest.suffix == ".gz":  # check the archive is complete before accepting it
        with gzip.open(tmp, "rb") as g:
            while g.read(1 << 22):
                pass
    tmp.replace(dest)


def run() -> dict:
    manifest_path = RAW / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    for url, rel in SOURCES:
        dest = RAW / rel
        if not dest.exists():
            print("downloading", rel, flush=True)
            fetch(url, dest)
            manifest[rel] = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        entry = manifest.setdefault(rel, {"url": url, "retrieved_utc": None})
        entry.update(sha256=sha256(dest), bytes=dest.stat().st_size)
    manifest_path.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    return manifest


if __name__ == "__main__":
    for k, v in run().items():
        print(f"{k}: {v['bytes']:,} bytes  sha256 {v['sha256'][:12]}")
