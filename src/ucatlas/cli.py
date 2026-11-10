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


def run_all(clean: bool = True, run_name: str | None = None) -> None:
    from . import cohorts, download, enrich, expression, kg_build, literature
    from .runs import current_run
    if clean:
        # the single-cell pseudobulk files (sc_*, sc3_*) come from long streaming steps and are kept; everything
        # else is rebuilt. They are rebuilt here only when missing.
        keep = ("sc_pseudobulk", "sc_celltype_share", "sc_celltype_de", "sc_cells", "sc3_pseudobulk", "sc3_celltype_share")
        for f in PROCESSED.glob("*"):
            if not f.name.startswith(keep):
                f.unlink() if f.is_file() else shutil.rmtree(f)
        GRAPH_DB.unlink(missing_ok=True)
    from .runs import RUNS
    run = RUNS / run_name if run_name else current_run()
    _step("download", lambda: f"{len(download.run())} files checked")
    _step("samples", lambda: f"{len(cohorts.build())} samples")
    _step("expression", lambda: f"{len(expression.run()['meta'])} genes in the meta-analysis")
    _step("literature", literature.ingest)
    _step("relations", literature.validate_extractions)
    _step("enrichment", enrich.run)
    _step("graph", lambda: kg_build.build()["claims_without_resolvable_evidence"])
    from . import composition, kg_cells, sc_replication, sc_neutrophil, singlecell
    if not (PROCESSED / "sc_celltype_share.parquet").exists() and (singlecell.SC / "manifest.json").exists():
        _step("single cell (Smillie 2019)", lambda: singlecell.run(expression_top(None)))
    if not (PROCESSED / "sc3_celltype_share.parquet").exists() and (sc_neutrophil.D / "manifest.json").exists():
        _step("single cell with neutrophils (Gubatan 2026)", sc_neutrophil.build)
    if (PROCESSED / "sc_celltype_share.parquet").exists():   # single-cell layer present (M11)
        top = expression_top(500)
        _step("single-cell replication", lambda: sc_replication.run(top))
        _step("single-cell graph", kg_cells.load)
        if (PROCESSED / "sc3_celltype_share.parquet").exists():
            _step("composition", lambda: (composition.run(expression_top(None))["verdict"].value_counts().to_dict(),
                                          composition.load())[0])
    from . import activity, discovery, pathways
    if (PROCESSED / "sc3_celltype_share.parquet").exists():
        _step("discovery", lambda: (discovery.run()["status"].value_counts().to_dict(), discovery.load())[0])
        _step("pathways", lambda: (pathways.pathway_activity(), pathways.tf_activity(), pathways.communication(),
                                   pathways.load())[3])
    _step("disease activity", lambda: (activity.run(), activity.load())[1])
    from . import treatment
    _step("treatment response", lambda: (treatment.response_genes(),
                                         treatment.interaction(list(dict.fromkeys(expression_top(300) + TX_EXTRA))),
                                         treatment.response_pathways(), treatment.load())[3])
    from . import spatial
    if (spatial.D / "manifest.json").exists():
        extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
        _step("spatial", lambda: (spatial.run(list(dict.fromkeys(expression_top(300) + extra))), spatial.load())[1])
    _step("convergence", lambda: _rebuild_convergence(run))
    _step("verify agents", lambda: _script("scripts/verify.py", "agents", run.name))
    _step("assemble", lambda: _script("scripts/assemble.py", run.name))
    _step("verify final", lambda: _script("scripts/verify.py", "final", run.name))
    _step("render", lambda: _script("scripts/render.py", run.name))
    _step("validate report", lambda: _script("scripts/validate_report.py"))


TX_EXTRA = ["TNF", "IL23A", "IL12B", "ITGB7", "ITGA4", "MADCAM1", "OSM", "OSMR", "TREM1", "JAK1", "S1PR1",
            "TNFSF15", "IL11"]


def expression_top(n):
    import pandas as pd
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    m = m.reindex(m["meta_log2fc"].abs().sort_values(ascending=False).index)
    return (m["gene"].head(n) if n else m["gene"]).tolist()


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
    r.add_argument("--run", help="agent run folder under runs/ to re-verify and render (default: runs/.current)")
    a = ap.parse_args()
    if a.cmd == "run-all":
        run_all(clean=not a.keep, run_name=a.run)


if __name__ == "__main__":
    main()
