"""External validation in mouse models (added 2026-11-10, third round).

MGI genotype-phenotype annotations (MGI_PhenoGenoMP.rpt; ontology names from VOC_MammalianPhenotype.rpt; human
orthologs from HMD_HumanPhenotype.rpt). Only genotypes that alter one gene are used (one marker with a human ortholog; Cre and reporter
transgenes have none, so conditional knockouts count), so a phenotype can be attributed to that gene. Terms:
- spontaneous: colitis (MP:0002816), large intestinal inflammation (MP:0003304), crypt abscesses (MP:0000491);
- induced, more severe: increased susceptibility to induced colitis (MP:0008537) or to colitis-induced
  morbidity/mortality (MP:0013298);
- induced, protected: decreased susceptibility to induced colitis (MP:0008539) or to colitis-induced
  morbidity/mortality (MP:0013299).
Most alleles are loss of function, so "loss gives colitis" suggests the gene protects the mucosa (inhibiting it is a
risk), and "loss protects" suggests the gene drives colitis (inhibiting it is the therapeutic direction). Transgenic
gain-of-function alleles are rare in these calls and are not separated. Check: how often UC clinical targets have any
colitis phenotype compared with all genes with an MGI ortholog.
Outputs: data/processed/mouse_colitis.parquet, mouse_colitis_summary.json.
"""
from __future__ import annotations

import gzip
import json

import pandas as pd

from .paths import PROCESSED, RAW

SPONT = {"MP:0002816", "MP:0003304", "MP:0000491"}
WORSE = {"MP:0008537", "MP:0013298"}
PROTECT = {"MP:0008539", "MP:0013299"}


def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else open(path, encoding="utf-8", errors="replace")


def run() -> dict:
    cols = ["genotype", "alleles", "background", "mp", "pmids", "markers", "allele_ids"]
    with _open(RAW / "extra" / "mgi_phenogenomp.rpt") as f:
        g = pd.read_csv(f, sep="\t", header=None, names=cols, usecols=range(7), dtype=str)
    with _open(RAW / "extra" / "mgi_hmd_human.rpt") as f:
        h = pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["human", "mouse", "mgi"], dtype=str)
    human = h.dropna().drop_duplicates("mgi").set_index("mgi")["human"]
    # one gene altered: exactly one marker with a human ortholog (Cre and reporter transgenes have none, so a
    # conditional knockout of one gene is kept; double mutants are not)
    g = g[g.markers.notna()]
    mapped = g.markers.str.split(r"[|,]").map(lambda ms: sorted({human[m] for m in ms if m in human.index}))
    g = g[mapped.map(len) == 1].copy()
    g["gene"] = mapped[mapped.map(len) == 1].str[0]
    universe = set(h.human.dropna())
    rel = g[g.mp.isin(SPONT | WORSE | PROTECT)]
    rows = []
    for gene, d in rel.groupby("gene"):
        sp, wo, pr = d.mp.isin(SPONT), d.mp.isin(WORSE), d.mp.isin(PROTECT)
        loss_gives = sp.any() or wo.any()
        direction = ("mixed" if loss_gives and pr.any() else
                     "loss gives or worsens colitis" if loss_gives else "loss protects from colitis")
        rows.append({"gene": gene, "spontaneous_genotypes": int(d[sp].genotype.nunique()),
                     "worse_induced_genotypes": int(d[wo].genotype.nunique()),
                     "protected_genotypes": int(d[pr].genotype.nunique()), "direction": direction,
                     "pmids": ";".join(sorted({p for x in d.pmids.dropna() for p in str(x).split("|") if p})[:10]),
                     "alleles": ";".join(sorted(set(d.alleles.dropna()))[:6])})
    out = pd.DataFrame(rows).sort_values("gene")
    out.to_parquet(PROCESSED / "mouse_colitis.parquet", index=False)
    from .evaluate import clinical_targets
    clin = set(clinical_targets()["gene"]) & universe
    hit = set(out.gene)
    conv = pd.read_parquet(PROCESSED / "convergence_table.parquet")
    top100 = set(conv.sort_values("score", ascending=False).gene.head(100))
    tprof = set(conv.sort_values("score_target", ascending=False).gene.head(100))
    s = {"genes_with_colitis_phenotype": len(hit), "genes_with_ortholog": len(universe),
         "background_rate": round(len(hit) / len(universe), 4),
         "clinical_targets_rate": round(len(clin & hit) / len(clin), 3), "clinical_targets": len(clin),
         "clinical_targets_with_phenotype": sorted(clin & hit),
         "convergence_top100_rate": round(len(top100 & hit) / 100, 3),
         "target_profile_top100_rate": round(len(tprof & hit) / 100, 3),
         "direction_counts": out.direction.value_counts().to_dict()}
    (PROCESSED / "mouse_colitis_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "mouse_colitis.parquet")
    n = 0
    with store.batch():
        store.add_node("DATASET:MGI", "Dataset", "MGI mouse genotype-phenotype annotations (colitis terms)")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "database": "MGI", "spontaneous_colitis_genotypes": r.spontaneous_genotypes,
                   "more_severe_induced_colitis_genotypes": r.worse_induced_genotypes,
                   "protected_from_induced_colitis_genotypes": r.protected_genotypes, "direction": r.direction,
                   "pmids": r.pmids.split(";") if r.pmids else [], "alleles": r.alleles.split(";") if r.alleles else []}
            ev = store.add_evidence(evidence_id("mouse_colitis", r.gene), "mouse_colitis", "DATASET:MGI",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "mouse_phenotype_in", "DATASET:MGI", ev)
            n += 1
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "mouse_colitis.parquet").set_index("gene")
    for g in ["IL10", "IL12B", "IL23A", "TNF", "CXCR2", "OSMR", "IDO1", "CHI3L1", "MMP1", "DUOX2", "CD274", "TNFRSF6B",
              "S100A8", "LCN2", "IRF5", "CHRDL2", "XKR9", "CPXM1", "ITGB7", "JAK1", "STAT3", "NOD2"]:
        if g in d.index:
            print(g, d.loc[g, ["spontaneous_genotypes", "worse_induced_genotypes", "protected_genotypes", "direction"]].to_dict())
        else:
            print(g, "no colitis phenotype")


# ------------------------------------------------------------------ literature (added 2026-11-10, fourth round)
MOUSE_QUERY = ('(TITLE:"{g}" OR ABSTRACT:"{g}") AND (TITLE:colitis OR ABSTRACT:colitis) AND '
               '(ABSTRACT:mice OR ABSTRACT:mouse OR ABSTRACT:murine) AND (ABSTRACT:knockout OR ABSTRACT:deficient OR '
               'ABSTRACT:deficiency OR ABSTRACT:"-/-" OR ABSTRACT:conditional OR ABSTRACT:antibody OR ABSTRACT:inhibitor '
               'OR ABSTRACT:blockade OR ABSTRACT:antagonist)')


def literature_genes(run=None) -> list[str]:
    """Report candidates, convergence top 30 and every first-in-class gene."""
    from .runs import current_run
    run = run or current_run()
    f = json.loads((run / "findings.json").read_text(encoding="utf-8"))
    g = [t["gene"] for t in f["targets"]] + [r["gene"] for r in f["convergence"]["ranking"][:30]]
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    return list(dict.fromkeys(g + d[d.status == "first-in-class opportunity"].gene.tolist()))


def literature_check(genes: list[str]) -> pd.DataFrame:
    """Europe PMC: papers on colitis in mice that perturb the gene (knockout, deficiency, blockade), title and abstract
    search (a gene symbol in the full text is not enough); the five most cited are listed. It says whether a gene was tested in mouse colitis, not with which result
    (MGI calls give the direction where curated). Answers cached (data/cache/europepmc)."""
    import json as _json
    from .external import EPMC_URL, _cached, _request
    rows = []
    for g in genes:
        q = MOUSE_QUERY.format(g=g)
        params = {"query": q, "format": "json", "pageSize": 5, "resultType": "lite", "sort": "CITED desc"}
        data = _cached("europepmc_mouse", _json.dumps(params, sort_keys=True),
                       lambda: _request("GET", EPMC_URL, params=params).json())["data"]
        hits = data.get("resultList", {}).get("result", [])
        rows.append({"gene": g, "papers": int(data.get("hitCount", 0)),
                     "top_pmids": [str(h.get("pmid")) for h in hits if h.get("pmid")][:5],
                     "top_titles": [h.get("title", "")[:160] for h in hits][:3]})
    d = pd.DataFrame(rows)
    d.to_parquet(PROCESSED / "mouse_literature.parquet", index=False)
    return d


def load_literature(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "mouse_literature.parquet")
    mgi = pd.read_parquet(PROCESSED / "mouse_colitis.parquet").set_index("gene")
    n = 0
    with store.batch():
        store.add_node("DATASET:EUROPEPMC_MOUSE", "Dataset", "Europe PMC search: gene perturbation in mouse colitis")
        for r in d.itertuples():
            node = store.resolve(r.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": r.gene, "mouse_colitis_papers": int(r.papers), "top_pmids": list(r.top_pmids),
                   "top_titles": list(r.top_titles), "mgi_direction": mgi.at[r.gene, "direction"] if r.gene in mgi.index else None,
                   "status": "tested in mouse colitis" if r.papers > 0 else "never tested in mouse colitis (by this search)"}
            ev = store.add_evidence(evidence_id("mouse_colitis_literature", r.gene), "mouse_colitis_literature",
                                    "DATASET:EUROPEPMC_MOUSE", date.today().isoformat(), **rec)
            store.add_claim(node["id"], "mouse_phenotype_in", "DATASET:EUROPEPMC_MOUSE", ev)
            n += 1
    return n
