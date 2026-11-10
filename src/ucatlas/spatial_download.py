"""M11 spatial: download only the expression files of the Xenium IBD series GSE312420 (Garrido-Trigo/consortia 2026).

The 34 GB archive is mostly morphology images; per sample, the gene x cell matrix, barcodes, features and the cell
table (centroids, areas) take about 160 MB in total. Sample metadata (tissue, disease, inflammation, panel) come from
the series matrix. Manifest with sha256 in data/raw/spatial/GSE312420/manifest.json.
"""
from __future__ import annotations

import gzip
import io
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


# ---------------------------------------------------------------------------------------------------------------
# Added 2026-11-10: two denser Xenium 5K sections and whole-transcriptome Visium sections.
# GSE339323: Xenium Prime 5K, one UC resection, most and least inflamed region (cell x gene h5 + cell table only).
# GSE189184: 10x Visium (whole transcriptome, 55 um spots), adult UC inflamed and healthy colon (checkpoint-inhibitor
# colitis sections are not used). Each sample is a 0.1-0.9 GB zip that is mostly images; only the filtered spot x
# gene matrix and the spot positions are read from it with HTTP range requests (about 1.5 MB per section).
XENIUM2 = {"GSM9893496": ("UC1_inf", "inflamed"), "GSM9893497": ("UC1_noninf", "noninflamed")}
VISIUM = {"GSM5695828": ("B5", "UC_inflamed"), "GSM5695829": ("B4", "UC_inflamed"), "GSM5695834": ("B8", "UC_inflamed"),
          "GSM5695835": ("B9", "UC_inflamed"), "GSM5695838": ("B12", "UC_inflamed"), "GSM5695839": ("B13", "UC_inflamed"),
          "GSM5695840": ("C2", "UC_inflamed"), "GSM5695837": ("B10", "healthy"), "GSM5695843": ("C5", "healthy")}
VISIUM_MEMBERS = ("filtered_feature_bc_matrix.h5", "spatial/tissue_positions_list.csv", "spatial/scalefactors_json.json")


class _HttpFile(io.RawIOBase):
    """Read-only seekable view of a remote file through HTTP range requests (enough for zipfile)."""

    def __init__(self, url: str):
        self.url, self.pos = url, 0
        self.size = int(requests.head(url, timeout=60, allow_redirects=True).headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        r = requests.get(self.url, headers={"Range": f"bytes={self.pos}-{end}"}, timeout=300)
        r.raise_for_status()
        self.pos += len(r.content)
        return r.content

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def _gsm_url(gsm: str, name: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"


def run_extra() -> dict:
    import zipfile
    out = {}
    # Xenium Prime 5K (GSE339323)
    d2 = RAW / "spatial" / "GSE339323"
    d2.mkdir(parents=True, exist_ok=True)
    m2p = d2 / "manifest.json"
    m2 = json.loads(m2p.read_text()) if m2p.exists() else {}
    for gsm, (stem, _) in XENIUM2.items():
        for suffix in ("cell_feature_matrix.h5", "cells.parquet.gz"):
            name = f"{gsm}_{stem}_{suffix}"
            dest = d2 / name
            if not dest.exists():
                fetch(_gsm_url(gsm, name), dest)
                m2[name] = {"url": _gsm_url(gsm, name), "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
            m2.setdefault(name, {}).update(sha256=sha256(dest), bytes=dest.stat().st_size)
    m2p.write_text(json.dumps(m2, indent=1, sort_keys=True))
    out["GSE339323"] = len(m2)
    # Visium (GSE189184), members read from the remote zips
    d3 = RAW / "spatial" / "GSE189184"
    d3.mkdir(parents=True, exist_ok=True)
    m3p = d3 / "manifest.json"
    m3 = json.loads(m3p.read_text()) if m3p.exists() else {}
    for gsm, (stem, _) in VISIUM.items():
        url = _gsm_url(gsm, f"{gsm}_{stem}.zip")
        targets = {mem: d3 / f"{gsm}_{mem.split('/')[-1]}" for mem in VISIUM_MEMBERS}
        if not all(p.exists() for p in targets.values()):
            z = zipfile.ZipFile(_HttpFile(url))
            names = {n.split("/", 1)[1]: n for n in z.namelist() if "/" in n}
            for mem, dest in targets.items():
                dest.write_bytes(z.read(names[mem]))
                m3[dest.name] = {"url": f"{url} (zip member {names[mem]})",
                                 "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        for dest in targets.values():
            m3.setdefault(dest.name, {}).update(sha256=sha256(dest), bytes=dest.stat().st_size)
        print("visium", gsm, flush=True)
    m3p.write_text(json.dumps(m3, indent=1, sort_keys=True))
    out["GSE189184"] = len(m3)
    return out
