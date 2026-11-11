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
        # The graph is cleared first: on Windows an open handle (e.g. an MCP server of a running session) blocks
        # deleting the file, and failing after the processed data were deleted left a half-cleaned tree (2026-11-11).
        _clear_graph()
        # the single-cell pseudobulk files (sc_*, sc3_*) come from long streaming steps and are kept; everything
        # else is rebuilt. They are rebuilt here only when missing.
        keep = ("sc_pseudobulk", "sc_celltype_share", "sc_celltype_de", "sc_cells", "sc3_pseudobulk", "sc3_celltype_share",
                "sc4_")
        for f in PROCESSED.glob("*"):
            if not f.name.startswith(keep):
                f.unlink() if f.is_file() else shutil.rmtree(f)
    from .runs import RUNS
    run = RUNS / run_name if run_name else current_run()
    from . import reproducibility
    if reproducibility.MANIFEST.exists():   # frozen inputs: fail before rebuilding if any input changed
        _step("frozen inputs", lambda: f"{reproducibility.require_verified()['checked']} files verified")
    _step("download", lambda: f"{len(download.run())} files checked")
    _step("samples", lambda: f"{len(cohorts.build())} samples")
    _step("expression", lambda: f"{len(expression.run()['meta'])} genes in the meta-analysis")
    _step("literature", literature.ingest)
    _step("relations", literature.validate_extractions)
    _step("enrichment", enrich.run)
    _step("graph", lambda: kg_build.build()["claims_without_resolvable_evidence"])
    _step("PubTator3 relations", literature.load_pubtator_relations)
    from . import antitnf_atlas, expression_adjusted
    _step("covariate-adjusted expression", lambda: (expression_adjusted.run(), expression_adjusted.load())[0])
    if antitnf_atlas.SRC.exists():
        _step("anti-TNF single-cell atlas", lambda: (antitnf_atlas.run(), antitnf_atlas.load())[1])
    from . import proteomics
    if proteomics.SRC.exists():
        _step("protein-level replication", lambda: (proteomics.run(), proteomics.load())[0])
    from . import composition, kg_cells, sc_replication, sc_neutrophil, singlecell
    if not (PROCESSED / "sc_celltype_share.parquet").exists() and (singlecell.SC / "manifest.json").exists():
        _step("single cell (Smillie 2019)", lambda: singlecell.run(expression_top(None)))
    if not (PROCESSED / "sc3_celltype_share.parquet").exists() and (sc_neutrophil.D / "manifest.json").exists():
        _step("single cell with neutrophils (Gubatan 2026)", sc_neutrophil.build)
    from . import sc_gse282122
    if (sc_gse282122.D / "manifest.json").exists():   # fourth atlas: full GSE282122 matrices (added 2026-11-11)
        _step("download GSE282122", lambda: f"{len(download.run_gse282122())} files checked")
        if not (PROCESSED / "sc4_celltype_share.parquet").exists():
            _step("single cell, anti-TNF atlas (GSE282122)", sc_gse282122.run)
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
        if (PROCESSED / "sc4_celltype_share.parquet").exists():   # needs the discovery genes; before the shortlist
            from . import sc_gse282122
            _step("GSE282122 single-cell graph", sc_gse282122.load)
            _step("single-cell replication, four atlases",
                  lambda: (sc_replication.run4(sc_replication.genes4()), sc_replication.load4(), kg_cells.load())[0])
        _step("pathways", lambda: (pathways.pathway_activity(), pathways.tf_activity(), pathways.communication(),
                                   pathways.load())[3])
    _step("disease activity", lambda: (activity.run(), activity.load())[1])
    from . import treatment
    _step("treatment response", lambda: (treatment.response_genes(),
                                         treatment.interaction(list(dict.fromkeys(expression_top(300) + TX_EXTRA))),
                                         treatment.response_pathways(), treatment.calprotectin(), treatment.load())[4])
    from . import spatial
    extra = ["IL11", "OSMR", "TREM1", "TNF", "IL23A", "ITGB7", "CXCR1", "CXCR2", "IL1B", "OSM", "FCGR3B"]
    if (spatial.D / "manifest.json").exists():
        _step("spatial", lambda: (spatial.run(list(dict.fromkeys(expression_top(300) + extra))), spatial.load())[1])
    from . import spatial_visium
    if (spatial_visium.D / "manifest.json").exists():
        _step("spatial (Visium)", lambda: (spatial_visium.run(list(dict.fromkeys(expression_top(300) + extra))),
                                           spatial_visium.load_graph())[1])
        from . import spatial_deconv
        _step("spatial (Visium deconvolution)",
              lambda: (spatial_deconv.run(list(dict.fromkeys(expression_top(300) + extra + spatial_deconv.EXTRA))),
                       spatial_deconv.viewer(spatial_deconv.VIEWER_GENES), spatial_deconv.load_graph())[2])
    from . import spatial_cosmx
    if (spatial_cosmx.D / "manifest.json").exists():
        _step("spatial (CosMx)", lambda: (spatial_cosmx.run(list(dict.fromkeys(expression_top(300) + extra + spatial_deconv.EXTRA))),
                                          spatial_cosmx.load())[0])
    from . import mr
    if mr.EXTRACT.exists() and mr.GWAS.exists():
        _step("Mendelian randomisation (pQTL)", lambda: (mr.run(), mr.load())[0])
    from . import mr_colon
    if mr_colon.EXTRACT.exists() and mr.GWAS.exists():
        _step("Mendelian randomisation (gut eQTL)", lambda: (mr_colon.run(), mr_colon.load())[0])
        from . import mr_sensitivity
        _step("Mendelian randomisation sensitivity", lambda: (mr_sensitivity.run(), mr_sensitivity.load())[0])
    _step("convergence", lambda: _rebuild_convergence(run))   # also loads the citable ranking records
    from . import report_layers, structure
    _step("structure (AlphaFold pockets)", lambda: (structure.run(structure.targets(run)), structure.load_graph())[1])
    _step("report layers", lambda: list(report_layers.run()))
    from . import mouse_models
    _step("mouse models (MGI)", lambda: (mouse_models.run(), mouse_models.load())[1])
    from . import safety
    _step("mouse colitis literature", lambda: (mouse_models.literature_check(mouse_models.literature_genes(run)),
                                               mouse_models.load_literature())[1])
    _step("target safety", lambda: (safety.run(), safety.load())[0])
    from . import shortlist_literature as sl_lit, validation_shortlist
    if any(sl_lit.CLASSIFIED.glob("*.jsonl")):   # classifier output (agent step, kept in the repository)
        _step("shortlist literature", lambda: ((sl_lit.DIR / "records.parquet").exists() or sl_lit.retrieve(),
                                               sl_lit.validate(), sl_lit.summarize(), sl_lit.load())[2])
    _step("validation shortlist", validation_shortlist.run)
    _step("graph statistics", _write_graph_stats)   # after every load, ranking included (the report shows these counts)
    _step("verify agents", lambda: _script("scripts/verify.py", "agents", run.name))
    _step("assemble", lambda: _script("scripts/assemble.py", run.name))
    _step("verify final", lambda: _script("scripts/verify.py", "final", run.name))
    _step("render", lambda: _script("scripts/render.py", run.name))
    _step("validate report", lambda: _script("scripts/validate_report.py"))


def _clear_graph() -> None:
    """Delete the graph file; if another process holds it open (Windows), empty its tables instead."""
    try:
        GRAPH_DB.unlink(missing_ok=True)
    except PermissionError:
        import sqlite3
        con = sqlite3.connect(GRAPH_DB, timeout=60)
        tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type = 'table'")]
        with con:
            for t in tables:
                con.execute(f'DELETE FROM "{t}"')
        con.close()
        print("graph file in use by another process: emptied its tables instead of deleting it", flush=True)


TX_EXTRA = ["TNF", "IL23A", "IL12B", "ITGB7", "ITGA4", "MADCAM1", "OSM", "OSMR", "TREM1", "JAK1", "S1PR1",
            "TNFSF15", "IL11"]


def expression_top(n):
    import pandas as pd
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    m = m.reindex(m["meta_log2fc"].abs().sort_values(ascending=False).index)
    return (m["gene"].head(n) if n else m["gene"]).tolist()


def _write_graph_stats() -> str:
    import json
    from .graph import GraphStore
    st = GraphStore().stats()
    (PROCESSED / "graph_stats.json").write_text(json.dumps(st, indent=1))
    return f"{st['claims_total']} claims, {st['claims_without_resolvable_evidence']} without evidence"


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
