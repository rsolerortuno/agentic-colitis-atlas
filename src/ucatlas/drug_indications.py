"""Which diseases the drugs of each gene are for (added 2026-11-12).

The discovery screen sets aside every gene with a clinical drug ("known druggable target"), and the convergence table
lists drugs without their indication. A drug that reached the clinic for another disease and was never tried in UC is
a repurposing lead: the target is reachable, its safety in people is partly known, and UC is open. This layer splits
the clinical record of each gene by indication, from Open Targets 26.09 clinical precedence (the same table as the
UC clinical-target set of evaluate.py):
- UC: a trial or approval whose disease is ulcerative colitis;
- Crohn's / IBD: Crohn's disease or inflammatory bowel disease, but not UC;
- other diseases: clinical drugs, none in UC or IBD (repurposing);
- broad drugs only: every drug hits more than five targets (kinase inhibitors and similar), so the gene is not what the
  drug was made for (MAX_TARGETS_PER_DRUG, as in evaluate.py); reported apart.
Stage = most advanced clinical stage (approval > phase 4 > 3 > 2/3 > 2 > 1/2 > 1). Trials stopped or withdrawn keep
their stage. An approval whose only record is an ATC code is shown as "ATC-coded, approval not confirmed": the ATC
classification is assigned to some drugs before or without approval (epacadostat, IDO1, is ATC-coded L01XX58 and
listed as approved, but was never approved).
Drug action (Open Targets directionOnTarget): LoF = the drug inhibits or blocks the target, GoF = it activates it,
missing = not stated (often antibody-drug conjugates, cell therapies or vaccines that kill cells carrying the protein).
The genes of the screen rise in UC, so a repurposing lead needs an inhibitor: other_action says whether one exists
(nitroglycerin activates GUCY1A1; sacituzumab govitecan carries a toxin to TACSTD2-positive cells: neither fits).
Trials outside Open Targets (added 2026-11-12): for the genes the target lists and the Top 100 can hold whose drugs
are only for other diseases, two more sources are searched for each drug (up to six per gene): ClinicalTrials.gov (the
drug as intervention, a condition naming colitis, Crohn's disease or IBD) and Europe PMC trial papers and protocols
(the drug in the title or abstract, colitis, Crohn's or IBD in the title, a trial publication type or "trial" or
"randomised" in the title). Calibrated: anakinra finds the IASO protocol of IL-1 blockade in acute severe colitis
(PMID 30772849, not in Open Targets or ClinicalTrials.gov); crizanlizumab, evolocumab, avacopan, iptacopan,
tilarginine and sotatercept find nothing. First run's false positives, now filtered: ipilimumab in colitis caused by
cancer immunotherapy (CTLA4), methylene blue as a chromoendoscopy dye (GUCY1A1), a suPAR biomarker paper (PLAU), and
other IL-18-pathway antibodies matched by a fuzzy registry search (IL18R1). A gene with such a trial moves to "tried in UC/IBD (outside Open Targets)"
and leaves the repurposing list.
Outputs: data/processed/drug_indications.parquet; graph evidence type clinical_indications (new type, so no record that
the agents cite changes).
"""
from __future__ import annotations

import json
import re

import pandas as pd
import requests

from .paths import PROCESSED, RAW

UC = re.compile(r"ulcerative colitis", re.I)
IBD = re.compile(r"crohn|inflammatory bowel", re.I)
MAX_TARGETS_PER_DRUG = 5
SALTS = (r"(?i)\s+(hydrochloride|dihydrochloride|sodium|potassium|calcium|mesylate|dimesylate|sulfate|sulphate|acetate|"
         r"maleate|tosylate|besylate|citrate|tartrate|fumarate|succinate|phosphate|bromide|anhydrous|hydrate|monohydrate)\b")
STAGES = ["APPROVAL", "ATC", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1", "EARLY_PHASE_1"]
RANK = {s: i for i, s in enumerate(STAGES)}
LABEL = {"APPROVAL": "approved", "ATC": "ATC-coded, approval not confirmed", "PHASE_4": "phase 4", "PHASE_3": "phase 3", "PHASE_2_3": "phase 2/3",
         "PHASE_2": "phase 2", "PHASE_1_2": "phase 1/2", "PHASE_1": "phase 1", "EARLY_PHASE_1": "early phase 1"}
CATEGORIES = ["UC", "Crohn's / IBD", "tried in UC/IBD (outside Open Targets)", "other diseases", "broad drugs only"]
IBD_COND = re.compile(r"ulcerative colitis|crohn|inflammatory bowel", re.I)
CANCER_COND = re.compile(r"neoplasm|cancer|carcinoma|tumou?r|malignan|melanoma|lymphoma", re.I)
TITLE_IBD = re.compile(r"ulcerative colitis|crohn|inflammatory bowel|severe colitis", re.I)
# colitis caused by cancer immunotherapy, and biomarker papers, are not trials of the drug in UC or IBD
TITLE_NOT = re.compile(r"checkpoint|immune[- ]related|inhibitor[- ]induced|ipilimumab|nivolumab|marker", re.I)


def _ctgov(drug: str) -> list[dict]:
    """ClinicalTrials.gov studies with the drug as intervention and a colitis, Crohn's or IBD condition (cached)."""
    from .external import _cached, _request
    params = {"query.intr": drug, "query.cond": "ulcerative colitis OR crohn disease OR inflammatory bowel disease",
              "pageSize": 20, "fields": "NCTId,BriefTitle,Phase,OverallStatus,Condition,InterventionName,DesignPrimaryPurpose"}
    d = _cached("clinical_trials", json.dumps(params, sort_keys=True),
                lambda: _request("GET", "https://clinicaltrials.gov/api/v2/studies", params=params).json())["data"]
    out = []
    for s in d.get("studies", []):
        p = s.get("protocolSection", {})
        conds = p.get("conditionsModule", {}).get("conditions", [])
        names = " | ".join(i.get("name", "") for i in p.get("armsInterventionsModule", {}).get("interventions", []))
        purpose = p.get("designModule", {}).get("designInfo", {}).get("primaryPurpose", "")
        # the search is fuzzy: keep studies that give this drug, in UC/Crohn's/IBD (not cancer), not as a diagnostic
        # (methylene blue is a chromoendoscopy dye in IBD surveillance)
        if (drug.lower() in names.lower() and any(IBD_COND.search(c) for c in conds)
                and not any(CANCER_COND.search(c) for c in conds) and purpose not in ("DIAGNOSTIC", "SCREENING")):
            out.append({"source": "ClinicalTrials.gov", "id": p["identificationModule"]["nctId"],
                        "title": p["identificationModule"].get("briefTitle", ""),
                        "phase": ", ".join(p.get("designModule", {}).get("phases", []) or [])})
    return out


def _epmc_trials(drug: str) -> list[dict]:
    """Europe PMC trial papers and protocols of the drug in colitis, Crohn's disease or IBD (cached)."""
    from .external import _cached, _request
    q = (f'(ABSTRACT:"{drug}" OR TITLE:"{drug}") AND (TITLE:colitis OR TITLE:crohn* OR TITLE:"inflammatory bowel") AND '
         '(PUB_TYPE:"Clinical Trial" OR PUB_TYPE:"Randomized Controlled Trial" OR PUB_TYPE:"Clinical Trial Protocol" OR '
         'TITLE:trial OR TITLE:randomised OR TITLE:randomized)')
    params = {"query": q, "format": "json", "pageSize": 10, "resultType": "lite"}
    d = _cached("europepmc_trials", json.dumps(params, sort_keys=True),
                lambda: _request("GET", "https://www.ebi.ac.uk/europepmc/webservices/rest/search", params=params).json())["data"]
    return [{"source": "Europe PMC", "id": f"PMID {r.get('pmid')}", "title": r.get("title", "")[:200],
             "phase": ""} for r in d.get("resultList", {}).get("result", [])
            if r.get("pmid") and TITLE_IBD.search(r.get("title", "")) and not TITLE_NOT.search(r.get("title", ""))]


def check_trials(d: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    """Move genes whose 'other disease' drugs were tried in UC/IBD outside Open Targets."""
    found = {}
    for r in d[(d.category == "other diseases") & d.gene.isin(genes)].itertuples():
        drugs = [x for x in dict.fromkeys(list(r.inhibitor_drugs or []) + list(r.other_drugs or [])) if len(x) <= 60][:6]
        hits = []
        for drug in drugs:   # long names are descriptions of cell or gene therapies, not searchable drug names
            for fn in (_ctgov, _epmc_trials):
                try:
                    hits += [{**h, "drug": drug} for h in fn(drug)]
                except requests.HTTPError:   # a name the service cannot parse: no result
                    pass
        if hits:
            found[r.gene] = hits[:5]
    d = d.copy()
    d["ibd_trials_elsewhere"] = d.gene.map(lambda g: json.dumps(found.get(g, []), ensure_ascii=False))
    d.loc[d.gene.isin(found), "category"] = "tried in UC/IBD (outside Open Targets)"
    return d


def _best(x: pd.DataFrame) -> tuple[str | None, list[str], list[str]]:
    """Most advanced stage, its drugs (up to three) and the diseases at that stage (up to three)."""
    x = x[x.clinicalStage.isin(RANK)]
    if x.empty:
        return None, [], []
    s = min(x.clinicalStage, key=RANK.get)
    top = x[x.clinicalStage == s]
    names = top.drugFromSource.fillna(top.drugId).str.replace(r"\s*\[.*?\]", "", regex=True).str.strip().str.title()
    names = names.str.replace(SALTS, "", regex=True).str.strip()   # one name per drug, not one per salt
    drugs = list(dict.fromkeys(names))[:3]
    dis = list(dict.fromkeys(top.diseaseFromSource.dropna().str.lower()))[:3]
    return s, drugs, dis


def run(genes: list[str] | None = None) -> pd.DataFrame:
    cp = pd.read_parquet(RAW / "extra" / "ot_clinical_precedence.parquet",
                         columns=["targetId", "drugId", "drugFromSource", "clinicalStage", "diseaseFromSource",
                                  "clinicalReportId", "directionOnTarget"])
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    cp["gene"] = cp.targetId.map(dict(zip(hg.ensembl_gene_id, hg.symbol)))
    cp = cp.dropna(subset=["gene"])
    atc = cp.clinicalReportId.fillna("").str.fullmatch(r"[a-z]\d\d[a-z]{2}\d\d")
    cp.loc[(cp.clinicalStage == "APPROVAL") & atc, "clinicalStage"] = "ATC"
    cp["broad"] = cp.drugId.map(cp.groupby("drugId").targetId.nunique()) > MAX_TARGETS_PER_DRUG
    if genes is not None:
        cp = cp[cp.gene.isin(set(genes))]
    dis = cp.diseaseFromSource.fillna("")
    cp["is_uc"], cp["is_ibd"] = dis.str.contains(UC), dis.str.contains(IBD) & ~dis.str.contains(UC)
    rows = []
    for g, x in cp.groupby("gene"):
        sel = x[~x.broad]
        uc_s, uc_d, _ = _best(sel[sel.is_uc])
        ibd_s, ibd_d, _ = _best(sel[sel.is_ibd])
        ot_s, ot_d, ot_dis = _best(sel[~sel.is_uc & ~sel.is_ibd])
        br_s, br_d, br_dis = _best(x[x.broad])
        oth = sel[~sel.is_uc & ~sel.is_ibd]
        inh_s, inh_d, inh_dis = _best(oth[oth.directionOnTarget == "LoF"])
        has_act = bool((oth.directionOnTarget == "GoF").any())
        action = ("inhibitor" if inh_s else "activator only" if has_act else "action not stated") if len(oth) else None
        if uc_s:
            cat = "UC"
        elif ibd_s:
            cat = "Crohn's / IBD"
        elif ot_s:
            cat = "other diseases"
        elif br_s:
            cat = "broad drugs only"
        else:
            continue
        rows.append({"gene": g, "category": cat,
                     "uc_stage": uc_s, "uc_drugs": uc_d, "ibd_stage": ibd_s, "ibd_drugs": ibd_d,
                     "other_stage": ot_s, "other_drugs": ot_d, "other_diseases": ot_dis,
                     "broad_stage": br_s, "broad_drugs": br_d, "other_action": action,
                     "inhibitor_stage": inh_s, "inhibitor_drugs": inh_d, "inhibitor_diseases": inh_dis,
                     "uc_trials": int(sel[sel.is_uc].shape[0])})
    d = pd.DataFrame(rows)
    d = check_trials(d, trial_check_genes())
    d.to_parquet(PROCESSED / "drug_indications.parquet", index=False)
    return d


def trial_check_genes() -> set[str]:
    """Genes the target lists and the Top 100 can hold (the trial search is per drug, so it is kept to these)."""
    out = set()
    p = PROCESSED / "discovery.parquet"
    if p.exists():
        disc = pd.read_parquet(p)
        out |= set(disc.gene[disc.status.isin(["known druggable target", "well studied in UC", "first-in-class opportunity"])])
    try:
        from .runs import current_run
        conv = json.loads((current_run() / "convergence.json").read_text(encoding="utf-8"))["ranking"][:100]
        out |= {r["gene"] for r in conv}
    except (FileNotFoundError, KeyError):
        pass
    return out


def inhibitor_text(r) -> str:
    """For a repurposing lead: is there an inhibitor in the clinic for another disease?"""
    if r.get("other_action") == "inhibitor":
        return (f"yes, {LABEL[r['inhibitor_stage']]} ({', '.join(r['inhibitor_drugs'])}; "
                f"{'; '.join(r['inhibitor_diseases'])})")
    if r.get("other_action") == "activator only":
        return "no: the clinical drugs activate it"
    return "no: action not stated (often cell-killing antibody conjugates, cell therapies or vaccines)"


def describe(r) -> str:
    """One line for the report: where the gene's drugs are, most advanced stage first."""
    if r["category"] == "UC":
        return f"UC: {LABEL[r['uc_stage']]} ({', '.join(r['uc_drugs'])})"
    if r["category"] == "tried in UC/IBD (outside Open Targets)":
        t = json.loads(r.get("ibd_trials_elsewhere") or "[]")
        first = t[0] if t else {}
        return (f"Tried in UC/IBD, not in Open Targets: {first.get('drug', '')} ({first.get('source', '')} "
                f"{first.get('id', '')}){'; ' + str(len(t) - 1) + ' more' if len(t) > 1 else ''}")
    if r["category"] == "Crohn's / IBD":
        return f"Crohn's / IBD only: {LABEL[r['ibd_stage']]} ({', '.join(r['ibd_drugs'])})"
    if r["category"] == "other diseases":
        return f"Not in UC; {LABEL[r['other_stage']]} for {'; '.join(r['other_diseases'])} ({', '.join(r['other_drugs'])})"
    return f"Not in UC; only broad drugs ({LABEL[r['broad_stage']]}: {', '.join(r['broad_drugs'])})"


def load(store=None) -> dict:
    """One clinical_indications record per gene of the meta-analysis with a clinical drug."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "drug_indications.parquet")
    genes = set(pd.read_parquet(PROCESSED / "meta.parquet")["gene"])
    n = 0
    with store.batch():
        if store.get_node("DATASET:OT_CLINICAL_PRECEDENCE") is None:
            store.add_node("DATASET:OT_CLINICAL_PRECEDENCE", "Dataset",
                           "Open Targets 26.09 clinical precedence (trials and approvals, all diseases)")
        for r in d[d.gene.isin(genes)].to_dict("records"):
            node = store.resolve(r["gene"])
            if not node or node["type"] != "Gene":
                continue
            rec = {k: (list(v) if hasattr(v, "__len__") and not isinstance(v, str) else v) for k, v in r.items()}
            rec["summary"] = describe(r)
            ev = store.add_evidence(evidence_id("clinical_indications", r["gene"]), "clinical_indications",
                                    "DATASET:OT_CLINICAL_PRECEDENCE", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "has_clinical_drug_in", "DATASET:OT_CLINICAL_PRECEDENCE", ev)
            n += 1
    return {"clinical_indications_records": n}


if __name__ == "__main__":
    d = run()
    print(d.category.value_counts().to_dict())
    print(json.dumps({g: describe(r) for g, r in d.set_index("gene").head(5).iterrows()}, indent=1))
