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


def test_neighbourhood_test_does_not_pass_clustered_patterns_placed_at_random():
    """Synthetic tissue with two large zones. Patterns made of patches placed at random are clustered but unrelated
    to the zones: few may pass (drawing random cells as the null passed a quarter to a half of them). A pattern
    placed inside one zone must pass. No data needed."""
    import numpy as np
    from scipy import sparse
    from ucatlas import spatial
    rng = np.random.default_rng(11)
    side = 300   # zones several tiles wide, so the association holds at every tile size the test uses
    gx, gy = np.meshgrid(np.arange(float(side)), np.arange(float(side)))
    xy = np.column_stack([gx.ravel(), gy.ravel()]) + rng.normal(0, 0.2, (side * side, 2))
    zone = (np.sin(xy[:, 0] / 60) + np.cos(xy[:, 1] / 75)) > 0
    lab = np.where(zone, "epithelial", "fibroblast").astype(object)

    def patches(inside=None, n=60, radius=4.0):
        v = np.zeros(len(xy), dtype=np.float32)
        done = 0
        while done < n:
            c = rng.random(2) * side
            near = np.flatnonzero(((xy - c) ** 2).sum(axis=1) < radius ** 2)
            if inside is not None and not zone[near].mean() > 0.9:
                continue
            v[near[rng.random(len(near)) < 0.5]] = 1
            done += 1
        return sparse.csr_matrix(v).T
    names = [f"NULL{i}" for i in range(40)] + ["REAL"]
    m = sparse.hstack([patches() for _ in range(40)] + [patches(inside=True)]).tocsr()
    nb = spatial.neighbourhood(m, names, lab, xy, names)
    null = nb[nb.gene != "REAL"]
    assert (null.p_enriched < 0.05).mean() < 0.15
    real = nb[(nb.gene == "REAL") & (nb.neighbour_lineage == "epithelial")].iloc[0]
    assert real.p_enriched < 0.05 and real.enrichment > 1.5


def test_viewer_rings_the_spots_the_test_used():
    """The tissue viewer ringed the top tenth of the expressing spots while the test used the top tenth of all spots."""
    vp, npq = PROCESSED / "spatial_viewer.json", PROCESSED / "spatial_visium_neighbors.parquet"
    if not (vp.exists() and npq.exists()):
        pytest.skip("Visium layer not built")
    v = json.loads(vp.read_text(encoding="utf-8"))
    nb = pd.read_parquet(npq).drop_duplicates("gene").set_index("gene")
    for g in ("S100A8", "LCN2"):
        ringed = sum(len(s["high"].get(g, [])) for s in v["sections"] if s["region"] == "UC inflamed")
        tested = int(nb.at[g, "n_high_spots"])        # the test leaves out spots without neighbours
        assert tested <= ringed <= tested * 1.05, (g, ringed, tested)
    assert all(r["padj"] is None or 0 <= r["padj"] <= 1 for rows in v["stats"].values() for r in rows)
