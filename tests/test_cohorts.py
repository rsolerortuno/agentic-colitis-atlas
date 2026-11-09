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
    assert (n(s, "GSE16879", "UC"), n(s, "GSE16879", "healthy")) == (24, 12)
    assert (n(s, "GSE73661", "UC"), n(s, "GSE73661", "healthy")) == (67, 12)


def test_response_labels(s):
    uc = s[s.group == "UC"]
    ifx = uc[(uc.cohort == "GSE16879")]["response"].value_counts()
    assert ifx[True] == 8 and ifx[False] == 16
    ust = uc[(uc.cohort == "GSE206285") & (uc.treatment == "ustekinumab")]
    assert len(ust) == 364
