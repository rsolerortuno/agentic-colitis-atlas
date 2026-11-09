"""M11 spatial: download only the expression files of the Xenium IBD series GSE312420 (Garrido-Trigo/consortia 2026).

The 34 GB archive is mostly morphology images; per sample, the gene x cell matrix, barcodes, features and the cell
table (centroids, areas) take about 160 MB in total. Sample metadata (tissue, disease, inflammation, panel) come from
the series matrix. Manifest with sha256 in data/raw/spatial/GSE312420/manifest.json.
"""
from __future__ import annotations

import gzip
import json
import re
from datetime import datetime, timezone

import requests

from .download import fetch, sha256
from .paths import RAW

GSE = "GSE312420"
DEST = RAW / "spatial" / GSE
KEEP = ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz", "cells.parquet.gz")


def file_list() -> list[tuple[str, int]]:
    txt = requests.get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/{GSE}/suppl/filelist.txt", timeout=60).text
    out = []
    for line in txt.splitlines()[2:]:
        parts = line.split("\t")
        if len(parts) >= 4 and parts[1].endswith(KEEP):
            out.append((parts[1], int(parts[3])))
    return out


def run() -> dict:
    DEST.mkdir(parents=True, exist_ok=True)
    stub = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/{GSE}"
    fetch(f"{stub}/matrix/{GSE}_series_matrix.txt.gz", DEST / f"{GSE}_series_matrix.txt.gz")
    mpath = DEST / "manifest.json"
    manifest = json.loads(mpath.read_text()) if mpath.exists() else {}
    for name, size in file_list():
        gsm = name.split("_")[0]
        dest = DEST / name
        if not dest.exists():
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
            fetch(url, dest)
            manifest[name] = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        manifest.setdefault(name, {}).update(sha256=sha256(dest), bytes=dest.stat().st_size)
    mpath.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    return {"files": len(manifest)}


def samples():
    """Sample table: gsm, title, tissue, disease, inflammation, panel (from the series matrix characteristics)."""
    import pandas as pd
    rows: dict[str, list[str]] = {}
    chars = []
    with gzip.open(DEST / f"{GSE}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key == "!Sample_characteristics_ch1":
                chars.append(vals)
            elif key in ("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1"):
                rows[key[8:]] = vals
    df = pd.DataFrame(rows)
    for col in chars:
        for i, v in enumerate(col):
            if ":" in v:
                k, val = v.split(":", 1)
                df.loc[i, k.strip()] = val.strip()
    return df


if __name__ == "__main__":
    print(run())
    print(samples().to_string(max_colwidth=40))
