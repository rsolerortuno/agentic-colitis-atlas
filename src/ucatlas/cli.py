"""`uv run ucatlas run-all`: deterministic rebuild from the raw downloads to the report (M9).

Steps: download (sha256 manifest) -> sample table -> expression and meta-analysis -> literature (frozen PMID lists,
re-validation of the saved extractions) -> external enrichment (disk cache) -> graph -> convergence of the latest
agent run (unchanged agent files) -> verify agents -> assemble -> verify final -> render -> validate report.
The agent step itself (/informe) is not re-run: its outputs are the frozen model responses of that run.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time

from .paths import GRAPH_DB, PROCESSED, ROOT


def _step(name, fn):
    t = time.time()
    out = fn()
    print(f"[{time.time() - t:6.1f} s] {name}: {out if not isinstance(out, dict) else {k: v for k, v in list(out.items())[:6]}}",
          flush=True)
    return out


def _script(*args) -> str:
    p = subprocess.run([sys.executable, *args], cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
    if p.returncode:
        raise SystemExit(f"{args[0]} failed:\n{p.stdout}\n{p.stderr}")
    return p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ""


def run_all(clean: bool = True) -> None:
    from . import cohorts, download, enrich, expression, kg_build, literature
    from .runs import current_run
    if clean:
        # the single-cell layer (sc_*) is rebuilt by its own milestone step (it streams 5 GB); everything else goes
        for f in PROCESSED.glob("*"):
            if not f.name.startswith("sc_"):
                f.unlink() if f.is_file() else shutil.rmtree(f)
        GRAPH_DB.unlink(missing_ok=True)
    run = current_run()
    _step("download", lambda: f"{len(download.run())} files checked")
    _step("samples", lambda: f"{len(cohorts.build())} samples")
    _step("expression", lambda: f"{len(expression.run()['meta'])} genes in the meta-analysis")
    _step("literature", literature.ingest)
    _step("relations", literature.validate_extractions)
    _step("enrichment", enrich.run)
    _step("graph", lambda: kg_build.build()["claims_without_resolvable_evidence"])
    _step("convergence", lambda: _rebuild_convergence(run))
    _step("verify agents", lambda: _script("scripts/verify.py", "agents", run.name))
    _step("assemble", lambda: _script("scripts/assemble.py", run.name))
    _step("verify final", lambda: _script("scripts/verify.py", "final", run.name))
    _step("render", lambda: _script("scripts/render.py", run.name))
    _step("validate report", lambda: _script("scripts/validate_report.py"))


def _rebuild_convergence(run) -> str:
    """Recompute convergence.json of an existing run without deleting its agent files (integrate.py would)."""
    from .integrate import write
    from .provenance import load_basis
    out = write(day=run.name)
    return f"top10 {out['top10']}; basis kept: {load_basis(run) is not None}"


def main() -> None:
    ap = argparse.ArgumentParser(prog="ucatlas")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run-all", help="rebuild everything from the raw downloads to the report")
    r.add_argument("--keep", action="store_true", help="do not delete processed data and the graph first")
    a = ap.parse_args()
    if a.cmd == "run-all":
        run_all(clean=not a.keep)


if __name__ == "__main__":
    main()
