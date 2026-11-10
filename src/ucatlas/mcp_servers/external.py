"""external-mcp: Open Targets, STRING and UniProt (disk-cached), with graph evidence ids when available."""
from __future__ import annotations

import pandas as pd
from mcp.server.fastmcp import FastMCP

from .. import external as ext
from ..config import load
from ..graph import evidence_id
from ..paths import PROCESSED, RAW
from .common import canonical_gene, clean, existing

mcp = FastMCP("external-mcp")


def _ensembl(symbol: str) -> str | None:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str)
    m = hg[hg["symbol"] == symbol]["ensembl_gene_id"].dropna()
    return m.iloc[0] if len(m) else None


@mcp.tool()
def open_targets_association(gene: str) -> dict:
    """Open Targets Platform overall association score between a gene and ulcerative colitis
    (MONDO_0005101), its rank among associated targets, and the evidence_id."""
    g = canonical_gene(gene)
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").reset_index(drop=True)
    sub = ot[ot["symbol"] == g]
    if sub.empty:
        return {"gene": g, "message": f"not among the top {len(ot)} Open Targets associations for ulcerative colitis"}
    r = sub.iloc[0]
    return clean({"gene": g, "score": round(r["score"], 4), "rank": int(sub.index[0]) + 1, "of": len(ot),
                  "evidence_id": existing(evidence_id("open_targets_assoc", r["ensembl"]))})


@mcp.tool()
def open_targets_target(gene: str) -> dict:
    """Open Targets target profile: tractability (small molecule, antibody, PROTAC), known drugs with
    maximum clinical stage, and Reactome pathways. Includes evidence ids for each drug/tractability record."""
    g = canonical_gene(gene)
    ens = _ensembl(g)
    if not ens:
        return {"gene": g, "message": "no Ensembl id"}
    t = ext.ot_target(ens)
    if not t:
        return {"gene": g, "message": "not in Open Targets"}
    t["tractability_evidence_id"] = existing(evidence_id("open_targets_tract", g))
    for d in t["drugs"]:
        d["evidence_id"] = existing(evidence_id("open_targets_drug", d["chembl"], g))
    t["pathways"] = t["pathways"][:15]
    return clean(t)


@mcp.tool()
def string_interactions(gene: str, limit: int = 10) -> list[dict]:
    """High-confidence (score >= 0.7) STRING protein interaction partners for a gene."""
    g = canonical_gene(gene)
    out = ext.string_partners(g, limit=limit)
    for p in out:
        p["evidence_id"] = existing(evidence_id("string", *sorted([g, p["partner"]])))
    return clean(out)


@mcp.tool()
def uniprot_entry(gene: str) -> dict:
    """UniProtKB/Swiss-Prot entry for the human gene: accession, protein name, function, subcellular location."""
    g = canonical_gene(gene)
    u = ext.uniprot(g)
    if not u:
        return {"gene": g, "message": "no reviewed UniProt entry"}
    u["evidence_id"] = existing(evidence_id("uniprot", u["accession"]))
    return clean(u)


@mcp.tool()
def disease_top_targets(n: int = 25) -> list[dict]:
    """Top Open Targets associated targets for ulcerative colitis (MONDO_0005101), with scores."""
    ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet").head(min(n, 200))
    return clean([{"gene": r.symbol, "score": round(r.score, 4),
                   "evidence_id": existing(evidence_id("open_targets_assoc", r.ensembl))} for r in ot.itertuples()])



@mcp.tool()
def mendelian_randomization(gene: str) -> dict:
    """Mendelian randomisation of plasma protein level on ulcerative colitis risk: UKB-PPP cis-pQTL instruments (Open
    Targets 26.09 fine-mapped credible sets) on the de Lange 2017 UC GWAS. Returns the Wald ratio or IVW estimate
    (log odds of UC per SD of protein), odds ratio, p, BH q over all tested proteins, number of instruments, and, for
    significant proteins, whether the pQTL colocalises with an IBD/UC GWAS signal (coloc H4 >= 0.8 or CLPP >= 0.01).
    beta > 0: genetically higher protein raises UC risk (inhibition would be protective); beta < 0: it lowers risk.
    A significant estimate without colocalisation may come from a neighbouring causal variant. Cite the evidence_id.
    The same design is repeated with gut eQTL instruments (gene expression in the colon wall rather than blood):
    colon_tissue_expression (GTEx v10 transverse/sigmoid colon, CEDAR colon and rectum; per SD of expression) and
    gut_cell_type_expression (IBDverse single-cell eQTL of gut biopsies; the strongest cell type, p Bonferroni-corrected
    for the number of cell types tested). Each block has its own evidence_id; a gene may be tested in one and not
    in another."""
    from ..graph import evidence_id as _eid
    from .common import store as _store
    g = canonical_gene(gene)
    eid = _eid("mendelian_randomization", g)
    rec = _store().get_evidence([eid])
    if not rec:
        plasma = {"message": "no cis-pQTL instrument for this protein in UKB-PPP (not tested)"}
    else:
        plasma = {**{k: v for k, v in rec[0].items() if k not in ("source", "retrieved", "date", "type")}, "evidence_id": eid}
    gut = {}
    # each block also carries the FinnGen R12 replication (independent UC GWAS) and the combined estimate
    for layer in ("colon_tissue", "gut_cell"):
        e2 = _eid("mendelian_randomization_gut", f"{g}|{layer}")
        r2 = _store().get_evidence([e2])
        gut[layer] = ({**{k: v for k, v in r2[0].items() if k not in ("source", "retrieved", "date", "type")},
                       "evidence_id": e2} if r2 else {"message": "no gut eQTL instrument for this gene (not tested)"})
    return clean({"gene": g, "plasma_protein": plasma, "colon_tissue_expression": gut["colon_tissue"],
                  "gut_cell_type_expression": gut["gut_cell"]})


@mcp.tool()
def mr_hits(n: int = 25) -> list[dict]:
    """Proteins with the strongest Mendelian-randomisation evidence for ulcerative colitis (q < 0.05 first, ordered by
    p): gene, odds ratio per SD, q, direction, instruments and colocalisation. Each row carries its evidence_id."""
    from ..graph import evidence_id as _eid
    r = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").sort_values("p").head(max(1, min(n, 100)))
    return clean([{"gene": x.gene, "odds_ratio_per_sd": round(x.or_per_sd, 3), "p": float(f"{x.p:.3g}"),
                   "q": float(f"{x.q:.3g}"), "n_instruments": int(x.n_instruments),
                   "direction": "higher protein raises UC risk" if x.beta > 0 else "higher protein lowers UC risk",
                   "colocalises_with_ibd_gwas": None if x.colocalises is None or pd.isna(x.colocalises) else bool(x.colocalises),
                   "evidence_id": existing(_eid("mendelian_randomization", x.gene))} for x in r.itertuples()])


@mcp.tool()
def mouse_phenotype(gene: str) -> dict:
    """Mouse models (MGI): single-gene genotypes (conditional knockouts included) with spontaneous colitis, more
    severe or less severe induced colitis; direction 'loss gives or worsens colitis' (the gene protects the mucosa;
    inhibiting it is a risk) or 'loss protects from colitis' (the gene drives colitis; inhibition is the therapeutic
    direction). Cite the evidence_id. No record means no colitis phenotype recorded, which includes genes never
    tested in colitis models."""
    from ..graph import evidence_id as _eid
    from .common import store as _store
    g = canonical_gene(gene)
    eid = _eid("mouse_colitis", g)
    r = _store().get_evidence([eid])
    if not r:
        return {"gene": g, "message": "no colitis phenotype recorded in MGI for a single-gene mouse model"}
    return clean({**{k: v for k, v in r[0].items() if k not in ("source", "date", "type")}, "evidence_id": eid})


@mcp.tool()
def mr_gut_hits(layer: str = "colon_tissue", n: int = 25) -> list[dict]:
    """Genes with the strongest Mendelian-randomisation evidence for UC from gut eQTL instruments (layer
    'colon_tissue': GTEx v10 colon and CEDAR; 'gut_cell': IBDverse single-cell eQTL), ordered by p: gene, tissue or
    cell type, odds ratio per SD of expression, q, direction, colocalisation. Each row carries its evidence_id."""
    from ..graph import evidence_id as _eid
    r = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    r = r[r["layer"] == layer].sort_values("p").head(max(1, min(n, 100)))
    return clean([{"gene": x.gene, "tissue_or_cell": x.label, "odds_ratio_per_sd": round(x.or_per_sd, 3),
                   "p": float(f"{x.p:.3g}"), "q": float(f"{x.q:.3g}"), "n_instruments": int(x.n_instruments),
                   "direction": "higher expression raises UC risk" if x.beta > 0 else "higher expression lowers UC risk",
                   "colocalises_with_ibd_gwas": None if x.colocalises is None or pd.isna(x.colocalises) else bool(x.colocalises),
                   "evidence_id": existing(_eid("mendelian_randomization_gut", f"{x.gene}|{layer}"))} for x in r.itertuples()])


if __name__ == "__main__":
    load()
    mcp.run()
