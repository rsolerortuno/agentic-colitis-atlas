"""Download the CosMx sections of GSE312415 used by the atlas (added 2026-11-10, third round).

GSE312415 (IBD spatial platform comparison) ships every sample in one 6.5 GB tar.gz. The archive is streamed once and
only the members of the UC and control rectum sections are written: cell x gene count matrices, cell metadata (position, area) and field-of-view positions; the
transcript and polygon files are not kept. CosMx 1K panel: 4 inflamed and 2 non-inflamed UC rectum biopsies; CosMx 6K
panel: 2 inflamed, 1 non-inflamed UC and 1 healthy rectum biopsy. Every member name is listed in
data/raw/spatial/GSE312415/members.txt and every written file in the manifest with its SHA-256.
"""
from __future__ import annotations

import json
import re
import tarfile
from datetime import datetime, timezone

import requests

from .download import sha256
from .paths import RAW

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312415/suppl/GSE312415_files.tar.gz"
D = RAW / "spatial" / "GSE312415"
# Slides hold two biopsies each. Slide stem in the archive -> its two sections: (sample id, GSM, panel, state).
SLIDES = {
    "I0284RI0276R2119": [("I0284", "GSM9345415", "1k", "inflamed"), ("I0276", "GSM9345414", "1k", "noninflamed")],
    "I0286RI0277R1119": [("I0286", "GSM9345413", "1k", "inflamed"), ("I0277", "GSM9345412", "1k", "inflamed")],
    "310R321R228": [("I0310", "GSM9345410", "1k", "noninflamed"), ("I0321", "GSM9345411", "1k", "inflamed")],
    "Cos1380R387RRNA6K4u2716": [("I0380", "GSM9345418", "6k", "inflamed"), ("I0387", "GSM9345419", "6k", "inflamed")],
    "Cos1400R355RRNA6K4u1716": [("I0400", "GSM9345416", "6k", "healthy"), ("I0355", "GSM9345417", "6k", "noninflamed")],
}
KEEP = re.compile(r"^(" + "|".join(SLIDES) + r")_(exprMat|metadata|fov_positions)_file\.csv$")


def run() -> dict:
    D.mkdir(parents=True, exist_ok=True)
    mp = D / "manifest.json"
    manifest = json.loads(mp.read_text(encoding="utf-8")) if mp.exists() else {}
    names, kept = [], 0
    with requests.get(URL, stream=True, timeout=600) as r:
        r.raise_for_status()
        r.raw.decode_content = False
        with tarfile.open(fileobj=r.raw, mode="r|gz") as tar:
            for m in tar:
                names.append(f"{m.name}\t{m.size}")
                if not m.isfile() or not KEEP.match(m.name) or (D / m.name).exists():
                    continue
                dest = D / m.name
                with tar.extractfile(m) as src, open(dest, "wb") as out:
                    while chunk := src.read(1 << 20):
                        out.write(chunk)
                manifest[dest.name] = {"url": f"{URL} (tar member {m.name})", "bytes": dest.stat().st_size,
                                       "sha256": sha256(dest),
                                       "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
                kept += 1
                print("kept", m.name, m.size, flush=True)
    (D / "members.txt").write_text("\n".join(names), encoding="utf-8")
    mp.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    return {"members": len(names), "kept": kept}


if __name__ == "__main__":
    print(run())
