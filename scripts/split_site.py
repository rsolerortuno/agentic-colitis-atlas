"""Build the GitHub Pages site from the single-file report (standard library only; runs in the Pages workflow).

reports/uc_report.html is self-contained so it opens offline. Online, the four heaviest blocks (Plotly, Cytoscape and
3Dmol libraries, the 3D structures and the tissue-viewer images, about 11 MB) become separate files: the libraries are
cached by the browser, and the structure and tissue data are fetched only when their tab is opened. Blocks are marked
in the template with data-asset="<file name>".

Usage: python scripts/split_site.py [report.html] [site_dir]
"""
import re
import sys
from pathlib import Path

BLOCK = re.compile(r'<script([^>]*?)\sdata-asset="([^"]+)"([^>]*)>(.*?)</script>', re.S)


def split(report: Path, site: Path) -> dict:
    html = report.read_text(encoding="utf-8")
    (site / "assets").mkdir(parents=True, exist_ok=True)
    sizes = {}

    def repl(m):
        before, name, after, body = m.groups()
        (site / "assets" / name).write_text(body.strip(), encoding="utf-8")
        sizes[name] = len(body)
        attrs = (before + after).strip()
        if "application/gzip-base64" in attrs:          # data block: fetched on demand by the page script
            return f'<script {attrs} data-src="assets/{name}"></script>'
        return f'<script {attrs} src="assets/{name}"></script>' if attrs else f'<script src="assets/{name}"></script>'

    out = BLOCK.sub(repl, html)
    (site / "index.html").write_text(out, encoding="utf-8")
    (site / ".nojekyll").touch()
    return {"index_mb": round(len(out) / 1e6, 2), "assets_mb": {k: round(v / 1e6, 2) for k, v in sizes.items()}}


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("reports/uc_report.html")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("site")
    print(split(src, dst))
