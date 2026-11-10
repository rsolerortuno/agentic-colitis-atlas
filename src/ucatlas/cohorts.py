"""M1: one tidy sample table for the eight bulk cohorts (data/processed/samples.parquet).

Columns: sample, cohort, patient, group (UC_inflamed | UC_noninflamed | UC | healthy | other), region, timepoint,
treatment, response (True/False/None, the cohort's own endpoint), mayo (endoscopic or total, see column note),
plus the cohort's raw fields in `extra` (JSON). Only colon / rectum biopsies of UC patients and controls are used.
"""
from __future__ import annotations

import gzip
import json
import re

import pandas as pd

from .paths import PROCESSED, RAW

COLON = {"Rectum", "Sigmoid", "LeftColon", "RightColon", "Cecum", "Transverse", "Colon", "colon"}


def series_meta(gse: str) -> pd.DataFrame:
    """Sample rows of a GEO series matrix: title, source, geo id and every 'key: value' characteristic."""
    rows: dict[str, list[str]] = {}
    chars: list[list[str]] = []
    with gzip.open(RAW / gse / f"{gse}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key == "!Sample_characteristics_ch1":
                chars.append(vals)
            elif key in ("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1"):
                rows[key[8:]] = vals
    df = pd.DataFrame(rows)
    for col in chars:
        for i, v in enumerate(col):
            if ":" in v:
                k, val = v.split(":", 1)
                df.loc[i, k.strip()] = val.strip()
    return df


def _mssccr() -> pd.DataFrame:
    m = series_meta("GSE193677")
    m["patient"] = m["title"].str.extract(r"(MSCCR_reGRID_\d+)")[0]
    grp = m["diseasetypere"].map({"UC.I": "UC_inflamed", "UC.NonI": "UC_noninflamed", "Control.NonI": "healthy"})
    out = pd.DataFrame({"sample": m["title"].str.split(",").str[0], "geo": m["geo_accession"], "cohort": "GSE193677",
                        "patient": m["patient"], "group": grp.fillna("other"), "region": m["regionre"],
                        "timepoint": "baseline", "treatment": None, "response": None,
                        "mayo": pd.to_numeric(m["ibdmesuc_mayo_score"], errors="coerce"),
                        "age": pd.to_numeric(m["study_eligibility_age_at_endo"], errors="coerce"),
                        "sex": m["demographics_gender"]})
    return out[out["region"].isin(COLON)]


def _unifi() -> pd.DataFrame:
    m = series_meta("GSE206285")
    yn = {"Y": True, "N": False}
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE206285",
                        "patient": m["donor id"],
                        "group": m["diagnosis"].map({"ulcerative colitis": "UC", "healthy control": "healthy"}),
                        "region": "Sigmoid", "timepoint": "baseline",
                        "treatment": m["treatment"].replace("NA", None).map(
                            lambda t: None if not isinstance(t, str) else ("placebo" if "Placebo" in t else "ustekinumab")),
                        "response": m["mucosal healing at week 8"].map(yn), "mayo": None, "age": None, "sex": None})
    out["response_endpoint"] = "mucosal healing at week 8"
    out["remission"] = m["clinical remission at week 8"].map(yn)
    return out


def _pursuit() -> pd.DataFrame:
    m = series_meta("GSE92415")
    healthy = m["treatment"].eq("Healthy") | m.get("disease").isna()
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE92415",
                        "patient": m["subject"], "group": healthy.map({True: "healthy", False: "UC"}),
                        "region": "Colon",
                        "timepoint": m["visit"].map({"Week 0": "baseline", "Week 6": "week 6"}),
                        "treatment": m["treatment"].map({"golimumab": "golimumab", "Placebo": "placebo"}),
                        "response": m["wk6response"].map({"Yes": True, "No": False}),
                        "mayo": pd.to_numeric(m["mayo score"], errors="coerce"), "age": None, "sex": None})
    out.loc[out["group"] == "healthy", ["timepoint", "treatment", "response", "patient"]] = ["baseline", None, None, None]
    out.loc[out["group"] == "healthy", "patient"] = out.loc[out["group"] == "healthy", "sample"]
    out["response_endpoint"] = "clinical response at week 6"
    return out


def _gse16879() -> pd.DataFrame:
    m = series_meta("GSE16879")
    keep = m["disease"].isin(["UC", "Control"])
    m = m[keep].reset_index(drop=True)
    when = m["before or after first infliximab treatment"].map(
        {"Before first infliximab treatment": "baseline", "After first infliximab treatment": "after infliximab",
         "Not applicable": "baseline"})
    pid = m["title"].str.extract(r"^(\w+?)(?:_|$)")[0]
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE16879",
                        "patient": m["title"].str.replace(r"(before|after|_\w+)$", "", regex=True),
                        "group": m["disease"].map({"UC": "UC", "Control": "healthy"}), "region": "Colon",
                        "timepoint": when,
                        "treatment": m["disease"].map({"UC": "infliximab"}),
                        "response": m["response to infliximab"].map({"Yes": True, "No": False}),
                        "mayo": None, "age": None, "sex": None, "title": m["title"], "_pid": pid})
    out["response_endpoint"] = "endoscopic and histological response (authors)"
    return out


def _gse73661() -> pd.DataFrame:
    m = series_meta("GSE73661")
    ind = m["study individual number"]
    week = m["week (w)"]
    tx = m["induction therapy_maintenance therapy"]
    group = week.eq("CO").map({True: "healthy", False: "UC"})
    drug = tx.map(lambda t: None if t == "CO" else ("infliximab" if "IFX" in t else "vedolizumab"))
    # response is written on the post-treatment samples ("UC R ..." / "UC NR ..."); carry it to the patient's week 0
    # The induction response is the label of the EARLIEST post-baseline sample (W4-6 infliximab, W6 or W12 vedolizumab);
    # later labels (W52) describe maintenance. "UC other" (early stop, placebo switch) carries no label.
    lab = m["title"].str.extract(r"_UC (NR|R) ", expand=False).map({"R": True, "NR": False})
    wk = week.str.extract(r"W(\d+)")[0].astype(float)
    wk = wk.where(~week.eq("W4_W6"), 5.0)
    lab_df = pd.DataFrame({"ind": ind, "wk": wk, "lab": lab}).dropna(subset=["lab"]).sort_values("wk")
    resp_by_patient = lab_df.groupby("ind")["lab"].first()
    induction = tx.str.split("_").str[0]
    drug = induction.map(lambda t: None if t == "CO" else ("infliximab" if t == "IFX" else
                                                           "placebo" if t == "plac" else "vedolizumab"))
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE73661",
                        "patient": "GSE73661_" + ind, "group": group, "region": "Colon",
                        "timepoint": week.map(lambda w: "baseline" if w in ("W0", "CO") else w),
                        "treatment": drug, "response": ind.map(resp_by_patient),
                        "mayo": pd.to_numeric(m["mayo endoscopic subscore"], errors="coerce"),
                        "age": None, "sex": None, "title": m["title"]})
    out.loc[out["group"] == "healthy", "response"] = None
    out["response_endpoint"] = "endoscopic response (authors, Mayo endoscopic subscore)"
    return out


def _gse109142() -> pd.DataFrame:
    """PROTECT: treatment-naive paediatric UC at diagnosis and non-IBD controls, rectal biopsies (RNA-seq)."""
    m = series_meta("GSE109142")
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE109142",
                        "patient": m["geo_accession"],
                        "group": m["diagnosis"].map({"Ulcerative Colitis": "UC", "Control": "healthy"}).fillna("other"),
                        "region": "Rectum", "timepoint": "baseline", "treatment": None, "response": None,
                        "mayo": pd.to_numeric(m["total mayo score"], errors="coerce"),
                        "age": pd.to_numeric(m["age at diagnosis"], errors="coerce"), "sex": m["Sex"]})
    out["response_endpoint"] = None
    return out


def _gse117993() -> pd.DataFrame:
    """RISK: paediatric UC at diagnosis and non-IBD controls, rectal biopsies (RNA-seq); Crohn's disease left out."""
    m = series_meta("GSE117993")
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE117993",
                        "patient": m["geo_accession"],
                        "group": m["diagnosis"].map({"UC": "UC", "Control": "healthy"}).fillna("other"),
                        "region": "Rectum", "timepoint": "baseline", "treatment": None, "response": None,
                        "mayo": None, "age": None, "sex": None})
    out["response_endpoint"] = None
    return out


def _gse165512() -> pd.DataFrame:
    """Adult UC and healthy colon biopsies (RNA-seq); ileum and Crohn's disease left out."""
    m = series_meta("GSE165512")
    colon = m["tissue"].eq("colon")
    grp = m["disease state"].map({"Ulcerative colitis": "UC", "Healthy": "healthy"}).where(colon).fillna("other")
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE165512",
                        "patient": "GSE165512_" + m["title"].str.extract(r"^(R\d+)")[0].fillna(m["geo_accession"]),
                        "group": grp, "region": "Colon", "timepoint": "baseline", "treatment": None, "response": None,
                        "mayo": None, "age": None, "sex": None})
    out["response_endpoint"] = None
    return out


def _gse72819() -> pd.DataFrame:
    """EUCALYPTUS (etrolizumab phase 2): baseline colonic biopsies, clinical remission at week 10 (RNA-seq). Treatment-
    response cohort only (no controls). The series holds the etrolizumab-treated patients; 3 without a label."""
    m = series_meta("GSE72819")
    rem = m["remission at week 10"].map({"Remitter": True, "Non-remitter": False})
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE72819",
                        "patient": "GSE72819_" + m["title"], "group": "UC", "region": "Colon", "timepoint": "baseline",
                        "treatment": "etrolizumab", "response": rem, "mayo": None, "age": None, "sex": None,
                        "title": m["title"], "prior_anti_tnf": m["prior anti.tnf"]})
    out["response_endpoint"] = "clinical remission at week 10 (EUCALYPTUS)"
    return out


def _gse23597() -> pd.DataFrame:
    """ACT1 (infliximab 5 or 10 mg/kg vs placebo, Toedter et al. 2011): colonic biopsies at weeks 0, 8 and 30; response
    at week 8 (clinical response, as given by the authors). Treatment-response cohort only (no controls)."""
    m = series_meta("GSE23597")
    drug = m["dose"].map(lambda d: "placebo" if d == "placebo" else "infliximab")
    out = pd.DataFrame({"sample": m["geo_accession"], "geo": m["geo_accession"], "cohort": "GSE23597",
                        "patient": "GSE23597_" + m["subject"], "group": "UC", "region": "Colon",
                        "timepoint": m["time"].map(lambda t: "baseline" if t == "W0" else t), "treatment": drug,
                        "response": m["wk8 response"].map({"Yes": True, "No": False}), "mayo": None, "age": None,
                        "sex": None, "title": m["title"]})
    out["response_endpoint"] = "clinical response at week 8 (ACT1)"
    return out


def build() -> pd.DataFrame:
    parts = [_mssccr(), _unifi(), _pursuit(), _gse16879(), _gse73661(), _gse109142(), _gse117993(), _gse165512(),
             _gse72819(), _gse23597()]
    df = pd.concat(parts, ignore_index=True)
    df["response"] = df["response"].astype("object").where(df["response"].notna(), None)
    keep = [c for c in df.columns if not c.startswith("_")]
    df = df[keep]
    df.to_parquet(PROCESSED / "samples.parquet", index=False)
    return df


def summary(df: pd.DataFrame) -> pd.DataFrame:
    base = df[df["timepoint"] == "baseline"]
    t = base.groupby(["cohort", "group"]).size().unstack(fill_value=0)
    return t


if __name__ == "__main__":
    d = build()
    print(summary(d).to_string())
