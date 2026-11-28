"""Citable ranks (rank_context.py): ranks within each named set are consistent with the values."""
import pytest

from ucatlas.paths import PROCESSED


@pytest.mark.skipif(not (PROCESSED / "discovery.parquet").exists(), reason="processed data not built")
def test_ranks_follow_the_values():
    from ucatlas import rank_context
    recs = rank_context.records()
    down = sorted((r for r in recs if r["set"] == "meta_down"), key=lambda r: r["rank"])
    assert down[0]["rank"] == 1 and down[0]["meta_log2fc"] == min(r["meta_log2fc"] for r in down)
    assert all(a["meta_log2fc"] <= b["meta_log2fc"] for a, b in zip(down, down[1:]))
    fic = [r for r in recs if r["set"] == "first_in_class"]
    assert {r["of"] for r in fic} == {len(fic)}
    assert sorted(r["rank_by_log2fc"] for r in fic)[0] == 1
    for name in ("mr_plasma", "mr_gut"):
        rows = [r for r in recs if r["set"] == name]
        assert len({r["gene"] for r in rows}) == len(rows)          # one record per gene
        best = min(rows, key=lambda r: r["p"])
        assert best["rank"] == 1


@pytest.mark.skipif(not (PROCESSED / "drug_indications.parquet").exists(), reason="processed data not built")
def test_drug_tags_say_where_the_drug_is():
    """A drug for another disease is never shown bare (fedratinib on JAK2 was read as a UC proposal), and a UC drug
    listed for several genes names them (IL-23 blockers under IL12B and IL23A)."""
    import pandas as pd
    from ucatlas.render import drug_tags
    tags = drug_tags()
    d = pd.read_parquet(PROCESSED / "drug_indications.parquet").set_index("gene")
    for g, tag in tags.items():
        cat = d.at[g, "category"]
        if cat == "UC":
            assert tag.startswith("UC ")
        elif cat == "other diseases":
            assert "(" in tag and tag.endswith(")")
    if "IL12B" in tags:
        assert "IL23A" in tags["IL12B"]
    if "JAK2" in tags:
        assert "Fedratinib" not in tags["JAK2"] and "JAK1" in tags["JAK2"]
