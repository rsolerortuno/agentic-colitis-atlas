"""M1: the sample table matches the counts checked against GEO (data/README.md)."""
import pandas as pd
import pytest

from ucatlas.paths import PROCESSED

pytestmark = pytest.mark.skipif(not (PROCESSED / "samples.parquet").exists(), reason="run ucatlas.cohorts first")


@pytest.fixture(scope="module")
def s():
    d = pd.read_parquet(PROCESSED / "samples.parquet")
    return d[d["timepoint"] == "baseline"]


def n(s, cohort, group, **kw):
    x = s[(s.cohort == cohort) & (s.group == group)]
    for k, v in kw.items():
        x = x[x[k] == v]
    return len(x)


def test_counts_match_geo(s):
    assert (n(s, "GSE193677", "UC_inflamed"), n(s, "GSE193677", "UC_noninflamed"), n(s, "GSE193677", "healthy")) == (293, 391, 340)
    assert (n(s, "GSE206285", "UC"), n(s, "GSE206285", "healthy")) == (550, 18)
    assert (n(s, "GSE92415", "UC"), n(s, "GSE92415", "healthy")) == (87, 21)
    assert (n(s, "GSE16879", "UC"), n(s, "GSE16879", "healthy")) == (24, 6)      # colon controls only
    assert (n(s, "GSE73661", "UC"), n(s, "GSE73661", "healthy")) == (67, 12)


def test_leuven_people_are_counted_once(s):
    """GSE73661 re-measures the infliximab patients and six controls of GSE16879; they stay out of its contrast."""
    b = s[(s.cohort == "GSE73661") & (s.timepoint == "baseline")]
    dup = b[b.in_gse16879]
    assert (dup.group == "UC").sum() == 23 and set(dup[dup.group == "UC"].treatment) == {"infliximab"}
    assert (dup.group == "healthy").sum() == 6
    assert not s[s.cohort != "GSE73661"].in_gse16879.any()
    import pandas as pd
    from ucatlas.paths import PROCESSED
    de = pd.read_parquet(PROCESSED / "de_results.parquet")
    n73 = de[de.cohort == "GSE73661"].iloc[0]
    n16 = de[de.cohort == "GSE16879"].iloc[0]
    assert (int(n73.n1), int(n73.n2)) == (44, 6) and (int(n16.n1), int(n16.n2)) == (24, 6)


def test_response_labels(s):
    uc = s[s.group == "UC"]
    ifx = uc[(uc.cohort == "GSE16879")]["response"].value_counts()
    assert ifx[True] == 8 and ifx[False] == 16
    ust = uc[(uc.cohort == "GSE206285") & (uc.treatment == "ustekinumab")]
    assert len(ust) == 364
