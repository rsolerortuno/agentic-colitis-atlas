"""Target safety (added 2026-11-10, fourth round). Before, safety was only "low RNA in vital organs" (discovery).

Three independent views, for the report candidates, the convergence top 40, the first-in-class genes and the UC
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


def genes() -> list[str]:
    from .evaluate import clinical_targets
    from .runs import current_run
    out = []
    f = current_run() / "findings.json"
    if f.exists():
        fj = json.loads(f.read_text(encoding="utf-8"))
        out += [t["gene"] for t in fj.get("targets", [])] + [r["gene"] for r in fj["convergence"]["ranking"][:40]]
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    out += d[d.status == "first-in-class opportunity"].gene.tolist()
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
        rows.append({"gene": g, "loeuf": loeuf, "pli": float(con.at[g, "lof.pLI"]) if g in con.index else np.nan,
                     "constrained": bool(loeuf < LOEUF_STRONG) if not np.isnan(loeuf) else None,
                     "safety_events": events, "n_safety_events": len(events),
                     "phewas_traits": [x["trait"] for x in pw], "n_phewas_traits": len(pw), "clinical_target": g in clin})
    d = pd.DataFrame(rows)
    d.to_parquet(PROCESSED / "safety.parquet", index=False)
    cand = d[~d.clinical_target]
    ct = d[d.clinical_target]
    s = {"genes": int(len(d)), "constrained_candidates": int(cand.constrained.fillna(False).sum()),
         "constrained_clinical_targets": int(ct.constrained.fillna(False).sum()), "clinical_targets": int(len(ct)),
         "candidates_with_safety_events": int((cand.n_safety_events > 0).sum()),
         "clinical_targets_with_safety_events": int((ct.n_safety_events > 0).sum()),
         "genes_with_phewas": int((d.n_phewas_traits > 0).sum()),
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
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene")
    print(d.loc[[g for g in ["CXCR2", "IL12B", "JAK2", "CXCL8", "MMP1", "CHI3L1", "DUOX2", "IDO1", "CD274", "CHRDL2", "LAX1",
                             "TNFRSF6B", "IL10", "TNF", "JAK1", "ITGB7"] if g in d.index],
                ["loeuf", "constrained", "n_safety_events", "phewas_traits"]].to_string())
