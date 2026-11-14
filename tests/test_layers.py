"""M14-M16 and review checks on the later layers (skip without processed data)."""
import json

import pandas as pd
import pytest

from ucatlas.paths import GRAPH_DB, PROCESSED

need = pytest.mark.skipif(not (PROCESSED / "response_genes.parquet").exists() or not GRAPH_DB.exists(),
                          reason="layers not built")


@need
def test_infliximab_positive_control_direction():
    """Arijs et al. 2009 infliximab non-response genes are higher in non-responders of GSE16879."""
    r = pd.read_parquet(PROCESSED / "response_genes.parquet")
    x = r[(r.arm == "infliximab") & r.gene.isin(["TNFRSF11B", "STC1", "PTGS2", "IL13RA2", "IL11"])]
    assert len(x) == 5 and (x.diff_nr_minus_r > 0).all() and (x.p < 0.01).all()


@need
def test_activity_tracks_inflammation_markers():
    a = pd.read_parquet(PROCESSED / "activity.parquet")
    a = a[(a.cohort == "GSE193677")].set_index("gene")
    assert a.loc["S100A8", "rho"] > 0.3 and a.loc["S100A8", "padj"] < 0.05


@need
def test_spatial_permutations_are_reproducible():
    """The neighbourhood permutation seed depends only on the gene (crc32), not on Python's hash randomisation."""
    import zlib
    import numpy as np
    assert np.random.default_rng(zlib.crc32(b"OSM")).random() == np.random.default_rng(zlib.crc32(b"OSM")).random()
    src = (PROCESSED.parents[1] / "src" / "ucatlas" / "spatial.py").read_text(encoding="utf-8")
    assert "hash(g)" not in src


@need
def test_graph_stats_are_current():
    from ucatlas.graph import GraphStore
    st = json.loads((PROCESSED / "graph_stats.json").read_text())
    assert st["claims_total"] == GraphStore().stats()["claims_total"]


@need
def test_tractability_covers_the_genome():
    from ucatlas.profiles import tractability_table
    t = tractability_table()
    assert len(t) > 15000 and t["predicted"].sum() > 3000


def test_pocket_control_ranks_in_the_texts_match_the_data():
    """The shortlist help says where the UC small-molecule targets rank by pocket score; it was typed and had drifted."""
    import json
    import re
    from ucatlas.paths import PROCESSED, ROOT
    p = PROCESSED / "structure_summary.json"
    if not p.exists():
        pytest.skip("structure layer not built")
    rows = [r for r in json.loads(p.read_text(encoding="utf-8"))["results"] if r.get("best_pocket_score") is not None]
    rows.sort(key=lambda r: -r["best_pocket_score"])
    ranks = [i + 1 for i, r in enumerate(rows) if r["control"]]
    text = (ROOT / "src/ucatlas/report_assets/shortlist_tab.html.j2").read_text(encoding="utf-8")
    m = re.search(r"small-molecule targets (\d+)th to (\d+)th of (\d+) proteins", text)
    assert m and (int(m.group(1)), int(m.group(2)), int(m.group(3))) == (min(ranks), max(ranks), len(rows))
