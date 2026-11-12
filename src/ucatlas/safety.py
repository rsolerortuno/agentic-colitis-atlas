"""Target safety (added 2026-11-10, fourth round). Before, safety was only "low RNA in vital organs" (discovery).

Three independent views, for the report candidates, the convergence top 100, every gene of the discovery screen that
the target lists can hold (first-in-class, known druggable, well studied in UC; widened 2026-11-12) and the UC
clinical targets (as reference):
1. Loss-of-function constraint (gnomAD v4.1, MANE Select transcript): LOEUF (upper bound of the observed / expected
   loss-of-function ratio). LOEUF < 0.35 = strongly constrained (people rarely carry a broken copy): life-long
   inhibition may be poorly tolerated; it does not forbid a drug (many drug targets are constrained).
2. Known safety liabilities (Open Targets target safety: curated adverse events and toxicity reports, with the
   direction of the intervention and the tissue); API answers cached in data/raw/extra/ot_safety_cache.json.
3. Genetic pleiotropy ("PheWAS" of the instruments): for genes with a significant Mendelian-randomisation estimate,
   other GWAS traits that colocalise (H4 >= 0.8) with the gene's instrument credible sets in Open Targets, IBD traits
   excluded, molecular traits (protein or metabolite levels) excluded; these are the other conditions a drug acting on the same gene could change. Uses the cached
   colocalisation answers of mr.py.
4. Mendelian disease (added 2026-11-11): the human diseases UniProt curates for the gene (DISEASE comments, from OMIM
   and the literature), kept only when UniProt states that the disease is caused by variants of the gene (not
   susceptibility loci, not somatic involvement), with the inheritance and whether the variants are gain of function. A gene whose loss causes
   a disease may give that disease when it is inhibited for years (CCNO and primary ciliary dyskinesia is the case
   that prompted it); a gain-of-function disease does not argue against inhibition. Cached through external.py.
Outputs: data/processed/safety.parquet, safety_summary.json.
"""
from __future__ import annotations

import json
import re
import time

import numpy as np
import pandas as pd
import requests

from .paths import PROCESSED, RAW

API = "https://api.platform.opentargets.org/api/v4/graphql"
CACHE = RAW / "extra" / "ot_safety_cache.json"
Q = """query($id:String!){ target(ensemblId:$id){ approvedSymbol safetyLiabilities{ event eventId
  effects{ direction dosing } biosamples{ tissueLabel cellLabel } datasource literature } } }"""
LOEUF_STRONG = 0.35
IBD = re.compile(r"ulcerative colitis|inflammatory bowel|crohn|colitis", re.I)
# molecular traits (protein or metabolite levels) are not conditions a drug would change
MOLECULAR = re.compile(r"level|circulating|protein|cerebrospinal|abundance|concentration", re.I)


def _ot(ens: str) -> list[dict]:
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    if ens not in cache:
        for attempt in range(5):
            try:
                r = requests.post(API, json={"query": Q, "variables": {"id": ens}}, timeout=120)
                if r.status_code == 200 and r.text.strip():
                    cache[ens] = ((r.json().get("data") or {}).get("target") or {}).get("safetyLiabilities") or []
                    break
            except (requests.RequestException, ValueError):
                pass
            time.sleep(4 * (attempt + 1))
        else:
            raise RuntimeError(f"Open Targets safety query failed for {ens}")
        CACHE.write_text(json.dumps(cache), encoding="utf-8")
    return cache[ens]


GOF = re.compile(r"gain[- ]of[- ]function|activating mutation|constitutive(ly)? activ", re.I)


def mendelian(symbol: str) -> list[dict]:
    from .external import UNIPROT_URL, _cached, _request
    params = {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true", "fields": "accession,cc_disease",
              "format": "json", "size": 1}
    data = _cached("uniprot", json.dumps(params, sort_keys=True),
                   lambda: _request("GET", UNIPROT_URL, params=params).json())["data"]
    out = []
    for r in data.get("results", [])[:1]:
        for c in r.get("comments", []):
            if c.get("commentType") != "DISEASE" or "disease" not in c:
                continue
            dz = c["disease"]
            note = " ".join(t.get("value", "") for t in (c.get("note") or {}).get("texts", []))
            if not re.search(r"disease is caused by", note, re.I):
                continue   # susceptibility loci and somatic involvement are not Mendelian diseases
            text = " ".join([dz.get("description", "")] + [t.get("value", "") for t in (c.get("note") or {}).get("texts", [])])
            inh = ("recessive" if re.search(r"recessive", text, re.I) else "dominant" if re.search(r"dominant", text, re.I)
                   else "X-linked" if re.search(r"X-linked", text, re.I) else "not stated")
            out.append({"disease": dz.get("diseaseId", ""), "inheritance": inh, "gain_of_function": bool(GOF.search(text))})
    return out


def genes() -> list[str]:
    from .evaluate import clinical_targets
    from .runs import current_run
    out = []
    f = current_run() / "findings.json"
    if f.exists():
        fj = json.loads(f.read_text(encoding="utf-8"))
        out += [t["gene"] for t in fj.get("targets", [])] + [r["gene"] for r in fj["convergence"]["ranking"][:100]]
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    # every gene the three target lists can hold (2026-11-12): before, only the first-in-class genes and the convergence
    # top 40 were checked, and list genes outside them showed "no Mendelian disease" without having been checked
    out += d[d.status.isin(["first-in-class opportunity", "known druggable target", "well studied in UC"])].gene.tolist()
    out += sorted(set(clinical_targets()["gene"]))
    return list(dict.fromkeys(out))


def phewas(gene: str) -> list[dict]:
    from . import mr
    traits = {}
    for f, layer in (("mr_instruments.parquet", "plasma"), ("mr_colon_instruments.parquet", "gut")):
        p = PROCESSED / f
        if not p.exists():
            continue
        h = pd.read_parquet(p, columns=["gene", "studyLocusId"])
        for sl in h[h.gene == gene].studyLocusId.unique():
            if not mr.COLOC_CACHE.exists() or sl not in json.loads(mr.COLOC_CACHE.read_text(encoding="utf-8")):
                continue
            data = mr._coloc_query(sl)
            rows = (((data.get("data") or {}).get("credibleSet") or {}).get("colocalisation") or {}).get("rows") or []
            for x in rows:
                st = (x.get("otherStudyLocus") or {}).get("study") or {}
                name = st.get("traitFromSource") or ""
                if (x.get("h4") or 0) >= 0.8 and name and not IBD.search(name) and not MOLECULAR.search(name):
                    traits[name] = max(traits.get(name, 0), round(float(x["h4"]), 3))
    return [{"trait": t, "h4": h} for t, h in sorted(traits.items(), key=lambda kv: -kv[1])][:12]


def run() -> dict:
    con = pd.read_csv(RAW / "extra" / "gnomad_v4.1_constraint.tsv", sep="\t", low_memory=False,
                      usecols=["gene", "gene_id", "mane_select", "lof.oe_ci.upper", "lof.pLI", "mis.z_score"])
    con = con[con.mane_select.astype(str).str.lower().isin(["true", "1"])].drop_duplicates("gene").set_index("gene")
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    ens = dict(zip(hg.symbol, hg.ensembl_gene_id))
    from .evaluate import clinical_targets
    clin = set(clinical_targets()["gene"])
    mr_sig = set()
    for f in ("mr_pqtl.parquet", "mr_colon.parquet"):
        r = pd.read_parquet(PROCESSED / f)
        mr_sig |= set(r[r.q < 0.05].gene)
    rows = []
    for g in genes():
        loeuf = float(con.at[g, "lof.oe_ci.upper"]) if g in con.index else np.nan
        liab = _ot(ens[g]) if g in ens else []
        events = sorted({x.get("event") for x in liab if x.get("event")})
        pw = phewas(g) if g in mr_sig else []
        md = mendelian(g)
        rows.append({"gene": g, "loeuf": loeuf, "pli": float(con.at[g, "lof.pLI"]) if g in con.index else np.nan,
                     "constrained": bool(loeuf < LOEUF_STRONG) if not np.isnan(loeuf) else None,
                     "safety_events": events, "n_safety_events": len(events),
                     "phewas_traits": [x["trait"] for x in pw], "n_phewas_traits": len(pw), "clinical_target": g in clin,
                     "mendelian_diseases": [f"{x['disease']} ({x['inheritance']}{', gain of function' if x['gain_of_function'] else ''})"
                                            for x in md],
                     "mendelian_loss_disease": any(not x["gain_of_function"] for x in md)})
    d = pd.DataFrame(rows)
    d.to_parquet(PROCESSED / "safety.parquet", index=False)
    cand = d[~d.clinical_target]
    ct = d[d.clinical_target]
    s = {"genes": int(len(d)), "constrained_candidates": int(cand.constrained.fillna(False).sum()),
         "constrained_clinical_targets": int(ct.constrained.fillna(False).sum()), "clinical_targets": int(len(ct)),
         "candidates_with_safety_events": int((cand.n_safety_events > 0).sum()),
         "clinical_targets_with_safety_events": int((ct.n_safety_events > 0).sum()),
         "genes_with_phewas": int((d.n_phewas_traits > 0).sum()),
         "candidates_with_mendelian_disease": int(cand.mendelian_loss_disease.sum()),
         "clinical_targets_with_mendelian_disease": int(ct.mendelian_loss_disease.sum()),
         "median_loeuf_candidates": round(float(cand.loeuf.median()), 3),
         "median_loeuf_clinical_targets": round(float(ct.loeuf.median()), 3)}
    (PROCESSED / "safety_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "safety.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:SAFETY", "Dataset", "Target safety: gnomAD v4.1 constraint, Open Targets safety liabilities, MR PheWAS")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "loeuf": None if np.isnan(r.loeuf) else round(float(r.loeuf), 3),
                   "strongly_constrained": r.constrained, "known_safety_events": list(r.safety_events)[:10],
                   "pleiotropic_traits_sharing_the_genetic_signal": list(r.phewas_traits)[:10],
                   "interpretation": "constraint and pleiotropy flag possible on-target risks; they do not forbid a drug"}
            ev = store.add_evidence(evidence_id("target_safety", r.gene), "target_safety", "DATASET:SAFETY",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "has_safety_profile_in", "DATASET:SAFETY", ev)
            n += 1
            if "mendelian_diseases" in d.columns and len(r.mendelian_diseases):
                # a separate record, so the target_safety records cited by earlier agent runs stay unchanged
                ev2 = store.add_evidence(evidence_id("mendelian_disease", r.gene), "mendelian_disease", "DATASET:SAFETY",
                                         date.today().isoformat(), gene=r.gene,
                                         diseases_caused_by_variants=list(r.mendelian_diseases)[:10],
                                         loss_of_function_disease=bool(r.mendelian_loss_disease), database="UniProt DISEASE",
                                         interpretation="inhibiting the gene for years may reproduce a disease caused by its loss")
                store.add_claim(node["id"], "has_mendelian_disease_in", "DATASET:SAFETY", ev2)
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene")
    print(d.loc[[g for g in ["CXCR2", "IL12B", "JAK2", "CXCL8", "MMP1", "CHI3L1", "DUOX2", "IDO1", "CD274", "CHRDL2", "LAX1",
                             "TNFRSF6B", "IL10", "TNF", "JAK1", "ITGB7"] if g in d.index],
                ["loeuf", "constrained", "n_safety_events", "phewas_traits"]].to_string())
