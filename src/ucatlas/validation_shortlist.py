"""Shortlist for experimental validation (added 2026-11-10, revised 2026-11-11): one row per first-in-class gene, with
a traffic light per criterion and a proposed first experiment. Everything comes from tables already in data/processed;
nothing is re-estimated.

Criteria (green / amber / red; grey = not measured):
- reproducible: cohorts with a significant increase (green >= 6 of 8, amber 4-5, red < 4), and the covariate-adjusted
  MSCCR change holds (if it does not, at most amber);
- protein: the protein also rises in the independent proteomics study (green: replicates; red: measured, does not
  replicate; grey: not quantified);
- tissue change: two views of "is it more RNA per cell, or more cells of the type that makes it?". (a) the bulk
  composition verdict of the discovery screen (intrinsic or partly composition, lineage level); (b) the change of the
  gene inside its own single-cell clusters (Smillie 2019, inflamed vs healthy, mean expression per cell in the three
  clusters that hold most of its RNA, weighted by that share; clusters with >= 50 cells in both states). (b) sees
  shifts inside a lineage that (a) cannot (crypt expansion, loss of mature enterocytes). Green: per-cell change
  >= 1 log2 with an intrinsic verdict, or >= 1.5 with a partly-composition verdict; red: per-cell change < 0.3 (the
  single cells explain the rise by cell numbers); amber otherwise. Calibrated on markers: DUOX2 +7.1, MMP3 +4.9,
  CHI3L1 +4.1 (state changes), TNF +0.2, IL1B -0.5, OSM -0.6 (rise by immune influx), AQP8 -2.7;
- cell source: share of the gene's transcripts in its dominant lineage in inflamed single cells (Smillie 2019, the
  lineage the discovery screen uses; green >= 0.5, amber 0.3-0.5, red < 0.3). Green also needs the lineage to
  replicate: since 2026-11-11 it must be the dominant lineage in at least 2 of the 3 other atlases (Garrido-Trigo
  2023, Gubatan 2026, GSE282122; sc_replication.run4), otherwise amber; before, a lineage not replicated in the second
  atlas (Garrido-Trigo 2023) capped it. A regeneration flag also caps it at amber: an epithelial gene with >= 76 % of its
  epithelial RNA in the crypt zone (stem, TA, cycling TA, secretory TA, enterocyte progenitors) that rises per cell
  there. 76 % is halfway between EPCAM (58 %, any epithelial gene) and the crypt markers LGR5, OLFM4, MKI67 and REG1A
  (92-97 %). Such a gene may belong to the repair programme of the crypt, and blocking it could slow healing.
  A crypt-niche flag (added 2026-11-11) does the same for stromal genes: the fibroblasts that sit under the crypt
  base and feed the stem cells (Smillie clusters WNT2B+ Fos-hi, WNT2B+ Fos-lo 1, WNT2B+ Fos-lo 2, RSPO3+; the
  Inflammatory Fibroblasts are the inflammation-associated state and are not niche) hold about half of the fibroblast
  cells (3,594 of 7,239 in inflamed mucosa), so a generic fibroblast gene has about half of its fibroblast RNA there:
  COL1A2 52 % (COL3A1 50 %, VIM 50 %, DCN 59 %). Known niche markers sit higher: WNT2B 62 %, GREM1 70 %, RSPO3 83 %
  (mean 72 %). The threshold, 62 %, is halfway between COL1A2 and the niche-marker mean, as the crypt threshold was
  set between EPCAM and the crypt markers. Flag: a stromal gene with >= 62 % of its fibroblast RNA (inflamed) in the
  niche clusters that rises per cell there (the per-cell measure of the tissue light restricted to the niche
  clusters, >= 1 log2). It caps the cell light at amber ("possible repair programme of the crypt niche") and adds an
  organoid-fibroblast co-culture read-out to the experiment. None of the 36 genes is flagged (2026-11-11): the group A
  stromal genes (CHRDL2 48 %, CPXM1 52 %, SGIP1 54 %, OLFML2B 50 %, PCDH17 13 %, an endothelial gene) have the niche
  share of a generic fibroblast gene, although all but PCDH17 rise per cell there (+2.0 to +3.4 log2);
- neutrophil-made genes (added 2026-11-12): when neutrophils make more than half of a gene's transcripts in inflamed
  mucosa (Gubatan 2026, the only atlas that captures them), the Smillie-based views cannot see the cell that makes it:
  the per-cell change is measured in other cells (ADGRE3 rose +2.4 log2 "in macrophages") and the lineage agreement
  across atlases folds neutrophils into myeloid cells, where three atlases have none. Such a gene is reported as
  neutrophil-directed (cell family "neutrophil"): the tissue light uses the per-cell change inside the neutrophil
  states of Gubatan 2026 instead (inflamed vs uninflamed UC mucosa, states with >= 50 cells in both, weighted by the
  gene's inflamed share; same thresholds; the baseline is uninflamed UC, not healthy colon, which Gubatan lacks).
  Checked: S100A8 +1.1, CXCR2 +0.7, FCGR3B +0.7 log2 per neutrophil; ADGRE3 -0.3 (its rise is all neutrophil
  arrival, red). The cell light is amber (source known from one atlas only, so it cannot replicate) and the first
  experiment is in blood neutrophils. This is a different kind of target (it acts on cell recruitment, like CXCR2),
  not a weaker one. Shortlist genes: ADGRE3, CEACAM4, CLEC1A, CSGALNACT2, RILPL2; clinical targets: CXCR2, MMP9,
  ALOX5, PDE4B, ICAM1, PTGS2. The cell family column also shows "endothelial" when vessels hold most of a gene's RNA
  (PCDH17), which the lineage grouping folds into stroma;
- causal genetics: a significant Mendelian-randomisation estimate in plasma or colon in the direction that supports
  inhibition (genetically higher level raises risk) is green, the opposite direction red, tested and null amber (it
  scores 0: a null test is not support), no instrument grey;
- druggable: an antibody route or strong small-molecule evidence in
  Open Targets (high-quality ligand, high-quality pocket, druggable family, clinical precedent) is green; a pocket
  score >= 0.6 alone, a weak bucket (structure with ligand) or a predicted membrane or signal peptide only is amber
  (the pocket score ranked the UC small-molecule targets 4th to 25th of 45 proteins: it does not discriminate alone).
  Antibody route = an extracellular part, from UniProt topology: an extracellular topological domain, a GPI anchor or
  secretion. "Membrane" or "Cell membrane" as a location is not enough (OSBPL3, EHD3 and RFTN1 sit on the cytoplasmic
  face); a protein marked secreted but also resident in the ER membrane (ERO1A) counts as amber. Checked on the UC
  clinical targets: ITGA2, S1PR1, CXCR2, ICAM1, IL2RA, ITGB7 (extracellular domain), IL23A, IL12B, CCL11, CXCL10
  (secreted) pass; JAK3 and PDE4B (small-molecule targets) do not;
- safety: highest vital-organ RNA (> 50 nTPM), loss-of-function constraint (LOEUF < 0.35), known safety events, and a
  Mendelian disease caused by loss of the gene (UniProt); green with none, amber with one flag, red with two or more;
- UC literature: what the papers about the gene in UC or IBD show (shortlist_literature.py: Europe PMC records read by
  a classifier subagent, every call backed by a verbatim sentence). Green: functional evidence of a pathogenic role
  (inhibition should help) and none protective; red: functional evidence of a protective role only (inhibition may
  harm); amber: a specific expression, biomarker or genetic result, or functional evidence both ways; grey: no
  paper about the gene, or only named in gene lists or reviews. Clinical targets in the calibration are not reviewed (grey). Only papers in a colitis or IBD
  setting count; functional results in other inflammatory settings are shown as context;
- mouse: whether the gene has ever been perturbed in mouse colitis (information, not a score).
Points: green 2, amber 1 (0 for a null genetic test), grey 0, red -1.
Groups, not positions (2026-11-11): many genes tie on points and a single weighting decides the order inside a tie, so
the table reports how stable each gene is across six ways of scoring (the base one; without the protein light;
without the literature light; without both; red lights costing 2; an unmeasured protein counted as amber). Stability
is the number of the six in which the gene is among the first ten (ties at the border included): group A, in all six;
group B, in one to five; group C, in none. The rule was fixed before it was applied. Inside a group genes are ordered by
stability, base points, greens and meta log2FC, and that order carries no extra information.
Protein bias: a measured protein can earn up to 2 points that an unmeasured one cannot, and mass spectrometry
quantifies abundant proteins more often, so the protein light favours abundant genes (quantified genes sit at the 78th
percentile of RNA abundance, the rest near the 40th, summary json). The "without protein" schemes show how much each gene
depends on it.
Calibration: the same rules applied to the UC clinical targets that are in the discovery table
(validation_shortlist_calibration.parquet). They are not expected to score high on every light (several rise by
immune influx, which is how an anti-inflammatory works); the comparison shows what each light separates.
Proposed experiment: by the concrete cell systems that hold the gene's RNA in inflamed mucosa (Smillie 2019 clusters
grouped into organoid epithelium, fibroblasts, endothelium, pericytes and glia, macrophages and DC, regulatory T cells,
T and NK cells, B and plasma cells, mast cells): the largest, and a second one when it holds at least 20 %; by route
(antibody, small molecule or genetic); with the published functional lead when a colitis study exists; followed by a
mouse colitis model when none exists.
Outputs: data/processed/validation_shortlist.parquet, validation_shortlist_calibration.parquet,
validation_shortlist_summary.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .paths import PROCESSED, RAW

ORDER = {"green": 2, "amber": 1, "grey": 0, "red": -1}
CRYPT = {"Stem", "TA 1", "TA 2", "Cycling TA", "Secretory TA", "Enterocyte Progenitors"}
CRYPT_SHARE = 0.76
NEUTROPHIL_MADE = 0.5   # share of transcripts made by neutrophils (Gubatan 2026, inflamed) above which a gene is neutrophil-directed
FIBROBLASTS = {"Inflammatory Fibroblasts", "WNT2B+ Fos-hi", "WNT2B+ Fos-lo 1", "WNT2B+ Fos-lo 2", "WNT5B+ 1",
               "WNT5B+ 2", "RSPO3+", "Myofibroblasts"}
NICHE = {"WNT2B+ Fos-hi", "WNT2B+ Fos-lo 1", "WNT2B+ Fos-lo 2", "RSPO3+"}   # crypt-niche fibroblasts (not inflammatory)
NICHE_SHARE = 0.62   # halfway between COL1A2 (0.52) and the mean of WNT2B, GREM1, RSPO3 (0.72); see the docstring
NICHE_MARKERS = {"niche": ("WNT2B", "GREM1", "RSPO3"), "generic": ("COL1A2", "COL3A1", "VIM", "DCN")}
STRONG_SM = ("High-Quality Ligand", "High-Quality Pocket", "Druggable Family", "Approved Drug", "Advanced Clinical",
             "Phase 1 Clinical")
CRITERIA = ["reproducible", "protein", "tissue", "cell", "genetics", "druggable", "safety", "literature", "mouse"]


def _lineage_share() -> dict[str, tuple[str, float]]:
    # the same source as the lineage of the discovery screen (Smillie 2019, inflamed), so the two always agree
    from .sc_replication import shares_2019
    lin = shares_2019("inflamed").sort_values("share", ascending=False).drop_duplicates("gene").set_index("gene")
    return {g: (r.lineage, float(r.share)) for g, r in lin.iterrows()}


class _SingleCell:
    """Per-cell change inside the gene's own clusters, the crypt-zone share and the crypt-niche share of the
    fibroblast RNA (Smillie 2019)."""

    def __init__(self):
        s = pd.read_parquet(PROCESSED / "sc_celltype_share.parquet")
        w = s.pivot_table(index=["gene", "cluster"], columns="health", values=["mean_cp10k", "n_cells", "share"])
        w.columns = [f"{a}_{b}" for a, b in w.columns]
        self.w = w.reset_index()
        x = s[(s.health == "inflamed") & (s.compartment == "epithelial")]
        self.crypt = (x[x.cluster.isin(CRYPT)].groupby("gene").cp10k_sum.sum() / x.groupby("gene").cp10k_sum.sum())
        f = s[(s.health == "inflamed") & s.cluster.isin(FIBROBLASTS)]
        self.niche = (f[f.cluster.isin(NICHE)].groupby("gene").cp10k_sum.sum() / f.groupby("gene").cp10k_sum.sum())

    def systems(self, g: str) -> list[tuple[str, float, str]]:
        """Cell systems holding the gene's RNA in inflamed mucosa: the largest, and a second one with >= 20 %."""
        d = self.w[self.w.gene == g].dropna(subset=["share_inflamed"]).copy()
        d["system"] = d.cluster.map(CLUSTER_SYSTEM)
        d = d.dropna(subset=["system"])
        if d.empty:
            return []
        tot = d.groupby("system").share_inflamed.sum().sort_values(ascending=False)
        out = []
        for name, share in tot.items():
            if out and (share < SYSTEM_MIN_SHARE or len(out) == 2):
                break
            top = d[d.system == name].nlargest(2, "share_inflamed").cluster.tolist()
            out.append((name, float(share), " and ".join(top)))
        return out

    def per_cell(self, g: str, k: int = 3, clusters: set | None = None) -> tuple[float, str]:
        d = self.w[self.w.gene == g].dropna(subset=["mean_cp10k_healthy", "mean_cp10k_inflamed"])
        if clusters is not None:
            d = d[d.cluster.isin(clusters)]
        d = d[(d.n_cells_healthy >= 50) & (d.n_cells_inflamed >= 50)].nlargest(k, "share_inflamed")
        if d.empty:
            return np.nan, ""
        wt = d.share_inflamed / d.share_inflamed.sum()
        lfc = np.log2((d.mean_cp10k_inflamed + 0.01) / (d.mean_cp10k_healthy + 0.01))
        return float((wt * lfc).sum()), ", ".join(f"{c} {v:+.1f}" for c, v in zip(d.cluster, lfc))


def _neutrophil_per_cell():
    """Per-cell change inside the neutrophil states of Gubatan 2026 (inflamed vs uninflamed UC), weighted by the
    gene's inflamed share; states with >= 50 cells in both conditions."""
    s = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet",
                        columns=["gene", "health", "cluster", "compartment", "mean_cp10k", "n_cells", "share"])
    s = s[s.compartment == "Neutrophil"]
    w = s.pivot_table(index=["gene", "cluster"], columns="health", values=["mean_cp10k", "n_cells", "share"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index().dropna(subset=["mean_cp10k_inflamed", "mean_cp10k_noninflamed"])
    w = w[(w.n_cells_inflamed >= 50) & (w.n_cells_noninflamed >= 50)]
    by_gene = dict(tuple(w.groupby("gene")))

    def get(g: str) -> tuple[float, str]:
        d = by_gene.get(g)
        if d is None or d.empty:
            return np.nan, ""
        d = d.sort_values("share_inflamed", ascending=False)
        lfc = np.log2((d.mean_cp10k_inflamed + 0.01) / (d.mean_cp10k_noninflamed + 0.01))
        wt = d.share_inflamed / d.share_inflamed.sum()
        return float((wt * lfc).sum()), ", ".join(f"{c} {v:+.1f}" for c, v in zip(d.cluster.head(3), lfc.head(3)))
    return get


def _topology(symbol: str) -> dict:
    from .external import UNIPROT_URL, _cached, _request
    params = {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
              "fields": "accession,ft_topo_dom,ft_transmem,ft_lipid,cc_subcellular_location", "format": "json", "size": 1}
    d = _cached("uniprot", json.dumps(params, sort_keys=True),
                lambda: _request("GET", UNIPROT_URL, params=params).json())["data"]
    if not d.get("results"):
        return {}
    r = d["results"][0]
    f = r.get("features", [])
    locs = [x["location"]["value"] for c in r.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION"
            for x in c.get("subcellularLocations", []) if "location" in x]
    return {"extracellular_domain": any(x["type"] == "Topological domain" and "xtracellular" in x.get("description", "")
                                        for x in f),
            "gpi": any(x["type"] == "Lipidation" and "GPI" in x.get("description", "") for x in f),
            "secreted": any(x.startswith("Secreted") for x in locs),
            "er_membrane": any(x.startswith("Endoplasmic reticulum membrane") for x in locs)}


def antibody_route(symbol: str) -> tuple[str, str]:
    """('full' | 'partial' | 'none', reason) from UniProt topology."""
    t = _topology(symbol)
    if t.get("extracellular_domain"):
        return "full", "extracellular domain"
    if t.get("gpi"):
        return "full", "GPI-anchored at the surface"
    if t.get("secreted"):
        return ("partial", "secreted, but also ER membrane") if t.get("er_membrane") else ("full", "secreted")
    return "none", ""


def _sm_labels():
    from . import external
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    ens = dict(zip(hg.symbol, hg.ensembl_gene_id))

    def get(g):
        t = external.ot_target(ens[g]) if g in ens else None
        return [x.split(": ", 1)[1] for x in (t or {}).get("tractability_labels", []) if x.startswith("SM: ")]
    return get


CELL_SYSTEMS = {   # cell system -> (first experiment, Smillie 2019 clusters that belong to it)
    "epithelium": ("knock-down or CRISPR knockout in patient-derived colon organoids stimulated with TNF and IFN-γ "
                   "(read-outs: barrier integrity, chemokine release, epithelial death)",
                   {"Stem", "TA 1", "TA 2", "Cycling TA", "Secretory TA", "Enterocyte Progenitors", "Enterocytes",
                    "Immature Enterocytes 1", "Immature Enterocytes 2", "Best4+ Enterocytes", "Goblet",
                    "Immature Goblet", "Tuft", "Enteroendocrine", "M cells"}),
    "fibroblasts": ("silencing in primary colonic fibroblasts from inflamed UC stimulated with OSM and IL-1β "
                    "(read-outs: IL-6, CXCL1, MMP release; neutrophil attraction in co-culture)",
                    {"Inflammatory Fibroblasts", "WNT2B+ Fos-hi", "WNT2B+ Fos-lo 1", "WNT2B+ Fos-lo 2", "WNT5B+ 1",
                     "WNT5B+ 2", "RSPO3+", "Myofibroblasts"}),
    "endothelium": ("silencing in primary human intestinal microvascular endothelial cells stimulated with TNF and "
                    "IL-1β (read-outs: ICAM1 and VCAM1, leukocyte adhesion and transmigration, permeability)",
                    {"Endothelial", "Microvascular", "Post-capillary Venules"}),
    "pericytes and glia": ("silencing in primary pericytes or enteric glial cells stimulated with IL-1β (read-outs: "
                           "IL-6 and chemokine release, support of endothelial tubes)", {"Pericytes", "Glia"}),
    "macrophages and DC": ("knock-down in monocyte-derived macrophages stimulated with LPS and IFN-γ (read-outs: IL-1β, "
                           "TNF, IL-23 and CXCL8 release, M1 polarisation)",
                           {"Macrophages", "Inflammatory Monocytes", "DC1", "DC2", "Cycling Monocytes"}),
    "regulatory T cells": ("knock-down or blockade in regulatory T cells (read-outs: suppression assay, FOXP3 stability, "
                           "IL-10)", {"Tregs"}),
    "T and NK cells": ("knock-down in blood or lamina propria CD4 and CD8 T cells (read-outs: activation, proliferation, "
                       "IFN-γ and IL-17 production, Th17 differentiation)",
                       {"CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "CD4+ Memory", "CD4+ PD1+", "CD8+ IELs",
                        "CD8+ IL17+", "CD8+ LP", "Cycling T", "ILCs", "NKs"}),
    "B and plasma cells": ("knock-down in activated B cells and plasmablasts (read-outs: activation, class switching, "
                           "IgA and IgG secretion)", {"Follicular", "GC", "Cycling B", "Plasma"}),
    "mast cells": ("knock-down in mast cells (read-outs: degranulation, cytokine release)", {"CD69+ Mast", "CD69- Mast"}),
    "neutrophils": ("blocking antibody or knock-down in blood neutrophils (read-outs: migration towards CXCL8, activation "
                    "markers, degranulation and NET formation; co-culture with colon organoids for epithelial damage)",
                    set()),   # not in Smillie 2019: chosen from the Gubatan 2026 neutrophil share
}
CLUSTER_SYSTEM = {c: name for name, (_, cl) in CELL_SYSTEMS.items() for c in cl}
LINEAGE_SYSTEM = {"epithelial": "epithelium", "stromal": "fibroblasts", "endothelial": "endothelium",
                  "myeloid": "macrophages and DC", "T/NK": "T and NK cells", "B": "B and plasma cells",
                  "plasma": "B and plasma cells", "mast": "mast cells"}
SYSTEM_MIN_SHARE = 0.2   # a second cell system is proposed when it holds at least 20 % of the gene's RNA


COCULTURE = (", with an organoid-fibroblast co-culture read-out (crypt organoid growth and stem-cell markers with "
             "silenced niche fibroblasts), since the gene rises in crypt-niche fibroblasts")


def _experiment(systems: list, lineage: str, surface: bool, sm: bool, mouse_papers: int | None, regen: bool,
                lead: str, niche: bool = False) -> str:
    """systems: [(name, share of the gene's RNA, its two main clusters)], from the gene's own single-cell clusters."""
    if not systems and LINEAGE_SYSTEM.get(lineage):
        systems = [(LINEAGE_SYSTEM[lineage], np.nan, "")]
    if systems:
        name, share, clus = systems[0]
        first = CELL_SYSTEMS[name][0]
        if name == "epithelium" and regen:
            first += ", with a wound-healing or organoid regrowth assay, since the gene rises in regenerating crypt cells"
        if name == "fibroblasts" and niche:
            first += COCULTURE
        where = f"{name}, {share:.0%} of its RNA, mainly {clus}" if not np.isnan(share) else name
        system = f"{first} [{where}]"
        if len(systems) > 1 and name == "neutrophils":
            system += (f". In the atlas without neutrophils (Smillie 2019) it is mostly in {systems[1][0]} "
                       f"({systems[1][1]:.0%} of the RNA there, mainly {systems[1][2]}): check that system too")
        elif len(systems) > 1:
            system += (f". Also in {systems[1][0]} ({systems[1][1]:.0%} of its RNA, mainly {systems[1][2]}) "
                       "with the matching system" + (COCULTURE if systems[1][0] == "fibroblasts" and niche else ""))
    else:
        system = "genetic tools first, in the cell type that expresses it"
    route = ("a blocking antibody can be tested directly" if surface else
             "small-molecule tool compounds may exist (Open Targets ligand or pocket evidence)" if sm else "genetic tools first")
    mouse = (" Then a knockout or blockade in DSS or T-cell transfer colitis (never tested in mouse colitis)."
             if mouse_papers == 0 else f" Mouse colitis data exist ({mouse_papers} paper{'s' if mouse_papers > 1 else ''}): check them first."
             if mouse_papers else "")
    return f"{system}; {route}.{(' Published lead: ' + lead) if lead else ''}{mouse}"


def _load():
    t = {"adj": pd.read_parquet(PROCESSED / "de_adjusted.parquet").set_index("gene"),
         "prot": pd.read_parquet(PROCESSED / "proteomics.parquet").set_index("gene"),
         "lin": _lineage_share(), "sc": _SingleCell(), "sm": _sm_labels(),
         "saf": pd.read_parquet(PROCESSED / "safety.parquet").set_index("gene"),
         "ml": pd.read_parquet(PROCESSED / "mouse_literature.parquet").set_index("gene")["papers"].to_dict()}
    rep_p = PROCESSED / "sc_replication.parquet"
    t["screp"] = pd.read_parquet(rep_p).set_index("gene")["same_lineage"].to_dict() if rep_p.exists() else {}
    rep4_p = PROCESSED / "sc_replication4.parquet"
    t["rep4"] = pd.read_parquet(rep4_p).set_index("gene") if rep4_p.exists() else pd.DataFrame()
    vis = json.loads((PROCESSED / "spatial_visium_summary.json").read_text())["targets"]
    t["vis"] = {v["gene"]: v for v in vis}
    mr_rows = {}
    for f_, lab in (("mr_pqtl.parquet", "plasma"), ("mr_colon.parquet", "colon")):
        r = pd.read_parquet(PROCESSED / f_)
        if "layer" in r:
            r = r[r.layer == "colon_tissue"]
        for x in r.itertuples():
            mr_rows.setdefault(x.gene, []).append((lab, float(x.beta), float(x.q)))
    t["mr"] = mr_rows
    s3 = pd.read_parquet(PROCESSED / "sc3_celltype_share.parquet", columns=["gene", "health", "cluster", "compartment", "share"])
    s3 = s3[(s3.health == "inflamed") & (s3.compartment == "Neutrophil")].sort_values("share", ascending=False)
    t["neu_states"] = {g: " and ".join(d.cluster.head(2)) for g, d in s3.groupby("gene", sort=False)}
    t["neu_pc"] = _neutrophil_per_cell()
    lit_p = PROCESSED / "shortlist_literature.parquet"
    t["lit"] = pd.read_parquet(lit_p).set_index("gene") if lit_p.exists() else pd.DataFrame()
    pk = pd.read_parquet(PROCESSED / "structure_pockets.parquet")
    t["pocket"] = pk[pk["rank"] == 1].set_index("gene")["pocket_score"].to_dict()
    return t


def _row(x, t) -> dict:
    g = x.gene
    c = {}
    rep = "green" if x.cohorts_up >= 6 else "amber" if x.cohorts_up >= 4 else "red"
    holds = bool(t["adj"].at[g, "holds"]) if g in t["adj"].index else None
    if holds is False and rep == "green":
        rep = "amber"
    c["reproducible"] = (rep, f"up in {int(x.cohorts_up)}/8 cohorts" +
                         ("" if holds is None else f"; adjusted: {'holds' if holds else 'does not hold'}"))
    prot = t["prot"]
    if g in prot.index and pd.notna(prot.at[g, "padj"]):
        r_ = prot.at[g, "replicates"]
        c["protein"] = ("green" if r_ == "yes" else "red" if r_ == "no" else "amber",
                        f"log2 ratio {prot.at[g, 'log2_ratio']:.2f}, padj {prot.at[g, 'padj']:.1e}")
    else:
        c["protein"] = ("grey", "not quantified")
    pc, pc_txt = t["sc"].per_cell(g)
    comp = x.composition or "not assessed"
    if np.isnan(pc):
        tis = "amber" if comp in ("intrinsic", "partly_composition") else "red"
    elif pc < 0.3:
        tis = "red"
    elif (comp == "intrinsic" and pc >= 1) or (comp == "partly_composition" and pc >= 1.5):
        tis = "green"
    else:
        tis = "amber"
    ret = f"{float(x.retained):.0%} kept" if pd.notna(x.retained) else ""
    neu = float(x.neutrophil_share) if pd.notna(getattr(x, "neutrophil_share", np.nan)) else 0.0
    neutro = neu > NEUTROPHIL_MADE
    if neutro:   # measure the per-cell change in the cell that makes it: neutrophil states of Gubatan 2026
        npc_, npc_txt = t["neu_pc"](g)
        if np.isnan(npc_):
            tis = "amber" if tis == "green" else tis
        else:
            tis = "red" if npc_ < 0.3 else "green" if ((comp == "intrinsic" and npc_ >= 1) or
                                                       (comp == "partly_composition" and npc_ >= 1.5)) else "amber"
        c["tissue"] = (tis, f"bulk: {comp.replace('_', ' ')}{', ' + ret if ret else ''}; neutrophils make {neu:.0%} of "
                            "its transcripts (Gubatan 2026); per neutrophil, inflamed vs uninflamed UC: " +
                       (f"{npc_:+.1f} log2 ({npc_txt})" if not np.isnan(npc_) else "not enough neutrophils") +
                       (f"; in the atlas without neutrophils (Smillie 2019): {pc:+.1f} log2 ({pc_txt})" if not np.isnan(pc) else ""))
    else:
        c["tissue"] = (tis, f"bulk: {comp.replace('_', ' ')}{', ' + ret if ret else ''}; per cell: " +
                       (f"{pc:+.1f} log2 ({pc_txt})" if not np.isnan(pc) else "no matched clusters"))
    lname, lshare = t["lin"].get(g, (x.lineage, np.nan))
    v = t["vis"].get(g, {})
    vtxt = f"; in situ {v['top_in_situ_expresser'].replace('_', ' ')}" if v.get("top_in_situ_expresser") else ""
    col = "grey" if np.isnan(lshare) else "green" if lshare >= 0.5 else "amber" if lshare >= 0.3 else "red"
    r4 = t["rep4"]
    atlases = {}
    if g in r4.index:   # four atlases: green needs the lineage dominant in >= 2 of the 3 others
        from .sc_replication import ATLASES
        n_agree = int(r4.at[g, "n_agree"])
        atlases = {ATLASES[k]: (r4.at[g, f"lineage_{k}"] if pd.notna(r4.at[g, f"lineage_{k}"]) else None,
                                None if pd.isna(r4.at[g, f"share_{k}"]) else float(r4.at[g, f"share_{k}"]))
                   for k in ATLASES}
        agree = [ATLASES[k].split(" (")[0] for k in ("garrido", "gubatan", "gse282122") if r4.at[g, f"agree_{k}"]]
        rtxt = f"; same lineage in {n_agree} of 3 other atlases" + (f" ({', '.join(agree)})" if agree else "")
        if col == "green" and n_agree < 2:
            col = "amber"
    else:
        n_agree = None
        sr = t["screp"]
        rtxt = "" if g not in sr else "; same lineage in the second atlas" if sr[g] else "; not replicated in the second atlas"
        if col == "green" and sr.get(g) is False:
            col = "amber"
    crypt = float(t["sc"].crypt.get(g, np.nan))
    regen = bool(lname == "epithelial" and crypt >= CRYPT_SHARE and (pc if not np.isnan(pc) else 0) >= 1)
    if regen and col == "green":
        col = "amber"
    ctxt = f"; {crypt:.0%} of epithelial RNA in the crypt zone, rising there: possible repair programme" if regen else ""
    nshare = float(t["sc"].niche.get(g, np.nan))
    npc, _ = t["sc"].per_cell(g, clusters=NICHE) if lname == "stromal" else (np.nan, "")
    niche = bool(lname == "stromal" and nshare >= NICHE_SHARE and (npc if not np.isnan(npc) else 0) >= 1)
    if niche and col == "green":
        col = "amber"
    if niche:
        ctxt += (f"; {nshare:.0%} of fibroblast RNA in crypt-niche fibroblasts, rising there ({npc:+.1f} log2): "
                 "possible repair programme of the crypt niche")
    if neutro:   # the source is known from the one atlas with neutrophils: it cannot replicate
        col = "amber" if col in ("green", "grey") else col
        c["cell"] = (col, f"neutrophil {neu:.0%} of transcripts (Gubatan 2026, the only atlas with neutrophils; "
                          f"cannot replicate); without neutrophils: {lname} {lshare:.0%} in Smillie 2019{rtxt}{vtxt}")
    else:
        c["cell"] = (col, (f"{lname} {lshare:.0%} of transcripts{rtxt}{ctxt}{vtxt}" if not np.isnan(lshare)
                           else "not in the single-cell atlas"))
    mrs = t["mr"].get(g, [])
    sig = [m for m in mrs if m[2] < 0.05]
    if sig:
        c["genetics"] = ("green" if any(m[1] > 0 for m in sig) else "red",
                         "; ".join(f"{m[0]} OR {np.exp(m[1]):.2f} (q {m[2]:.1e})" for m in sig))
    elif mrs:
        c["genetics"] = ("amber", "tested, no causal effect (scores 0)")
    else:
        c["genetics"] = ("grey", "no genetic instrument")
    ab, ab_txt = antibody_route(g)
    surface = ab == "full"
    sm = t["sm"](g)
    strong = [s for s in sm if s in STRONG_SM]
    weak = [s for s in sm if s not in STRONG_SM]
    pscore = t["pocket"].get(g)
    mods = list(x.tractable) if x.tractable is not None else []
    drug = "green" if surface or strong else "amber" if ab == "partial" or (pscore or 0) >= 0.6 or sm or mods else "red"
    parts = ([f"antibody route: {ab_txt}"] if ab != "none" else ["no extracellular part"]) + (["Open Targets: " + ", ".join(strong)] if strong else []) + \
            ([f"pocket {pscore:.2f}"] if pscore is not None else []) + ([f"weak: {', '.join(weak)}"] if weak else [])
    c["druggable"] = (drug, "; ".join(parts) or ("predicted " + "/".join(mods) if mods else "no route"))
    flags = []
    if pd.notna(x.vital_max_ntpm) and x.vital_max_ntpm > 50:
        flags.append(f"vital organs {x.vital_max_ntpm:.0f} nTPM")
    saf = t["saf"]
    if g in saf.index:
        if saf.at[g, "constrained"] is True:
            flags.append(f"constrained (LOEUF {saf.at[g, 'loeuf']:.2f})")
        if saf.at[g, "n_safety_events"] > 0:
            flags.append("known safety events")
        if "mendelian_loss_disease" in saf.columns and saf.at[g, "mendelian_loss_disease"]:
            flags.append("loss causes " + "; ".join(saf.at[g, "mendelian_diseases"]))
    vital = f"vital organs ≤ {x.vital_max_ntpm:.0f} nTPM" if pd.notna(x.vital_max_ntpm) else "vital organs not measured"
    c["safety"] = ("green" if not flags else "amber" if len(flags) == 1 else "red",
                   "; ".join(flags) if flags else f"{vital}, not constrained, no Mendelian disease")
    lit = t["lit"]
    c["literature"] = ((lit.at[g, "light"], lit.at[g, "verdict"])
                       if g in lit.index else ("grey", "not reviewed"))
    lead = ""
    if g in lit.index and "findings" in lit.columns:
        fl = [f for f in json.loads(lit.at[g, "findings"]) if f["kind"] == "functional" and f["ibd_setting"]
              and f["direction"] in ("pathogenic", "protective")]
        if fl:
            lead = f"{fl[0]['summary'].rstrip('.')} (PMID {fl[0]['pmid']})."
    papers = t["ml"].get(g)
    papers = int(papers) if papers is not None and pd.notna(papers) else None
    mouse_lit = lit.at[g, "mouse_colitis"] if g in lit.index and "mouse_colitis" in lit.columns else ""
    if mouse_lit:   # the literature review found a perturbation in mouse colitis that the title search missed
        papers = max(papers or 0, mouse_lit.count("PMID"))
    c["mouse"] = ("grey", "not searched" if papers is None else "never tested in mouse colitis" if papers == 0
                  else f"tested in mouse colitis: {mouse_lit}" if mouse_lit else f"{papers} mouse colitis papers")
    pts = {k: (0 if k == "genetics" and v[0] == "amber" else ORDER[v[0]]) for k, v in c.items() if k != "mouse"}
    systems = t["sc"].systems(g)
    if neutro:   # first system: neutrophils (Gubatan share), then the main Smillie system if it holds >= 20 %
        states = t["neu_states"].get(g, "neutrophil states")
        systems = [("neutrophils", neu, states)] + [s for s in systems[:1] if s[1] >= SYSTEM_MIN_SHARE]
    family = ("neutrophil" if neutro else "endothelial" if systems and systems[0][0] == "endothelium"
              and systems[0][1] >= 0.5 else x.lineage)
    return {"gene": g, "lineage": x.lineage, "cell_family": family, "neutrophil_directed": neutro,
            "neutrophil_share": round(neu, 3), "track": x.track, "status": x.status,
            "meta_log2fc": round(float(x.meta_log2fc), 2), "uc_papers": int(x.uc_papers),
            "criteria": json.dumps(c, ensure_ascii=False), "per_cell_log2fc": None if np.isnan(pc) else round(pc, 2),
            "crypt_share": None if np.isnan(crypt) else round(crypt, 3), "regeneration_flag": regen,
            "niche_share": None if np.isnan(nshare) else round(nshare, 3),
            "niche_per_cell_log2fc": None if np.isnan(npc) else round(npc, 2), "niche_flag": niche,
            "atlases_agree": n_agree, "atlas_lineages": json.dumps(atlases, ensure_ascii=False),
            "greens": sum(v[0] == "green" for k, v in c.items() if k != "mouse"),
            "reds": sum(v[0] == "red" for k, v in c.items() if k != "mouse"), "score": int(sum(pts.values())),
            "systems": json.dumps(systems, ensure_ascii=False),
            "experiment": _experiment(systems, x.lineage, surface, bool(strong) or (pscore or 0) >= 0.6, papers, regen, lead,
                                      niche)}


SCHEMES = {   # name -> arguments of _points
    "base": {},
    "without protein": {"drop": ("protein",)},
    "without literature": {"drop": ("literature",)},
    "without protein and literature": {"drop": ("protein", "literature")},
    "red lights cost 2": {"red": -2},
    "unmeasured protein as amber": {"grey_protein": 1},
}
TOP = 10


def _points(c: dict, drop: tuple = (), red: int = -1, grey_protein: int = 0) -> int:
    pts = 0
    for k in CRITERIA[:-1]:   # the mouse light is not scored
        if k in drop:
            continue
        v = c[k][0]
        if k == "genetics" and v == "amber":
            continue   # a null genetic test scores 0
        if k == "protein" and v == "grey":
            pts += grey_protein
            continue
        pts += red if v == "red" else ORDER[v]
    return pts


def _rank(rows: list[dict], stability: bool = False) -> pd.DataFrame:
    out = pd.DataFrame(rows)
    if stability:
        crit = [json.loads(c) for c in out.criteria]
        ranks = {}
        for name, kw in SCHEMES.items():
            pts = pd.Series([_points(c, **kw) for c in crit])
            ranks[name] = pts.rank(ascending=False, method="min")   # ties share the better position
        R = pd.DataFrame(ranks)
        out["stability"] = (R <= TOP).sum(axis=1).astype(int)
        out["best_rank"], out["worst_rank"] = R.min(axis=1).astype(int), R.max(axis=1).astype(int)
        out["group"] = np.where(out.stability == len(SCHEMES), "A", np.where(out.stability > 0, "B", "C"))
        out = out.sort_values(["stability", "score", "greens", "meta_log2fc"], ascending=False)
    else:
        out = out.sort_values(["score", "greens", "meta_log2fc"], ascending=False)
    out = out.reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def protein_bias() -> dict:
    """Do the genes that mass spectrometry quantified have more RNA than the rest? (percentile of RNA abundance)"""
    p = pd.read_parquet(PROCESSED / "proteomics.parquet")
    s = pd.read_parquet(PROCESSED / "sc_celltype_share.parquet", columns=["gene", "health", "cluster", "n_cells", "cp10k_sum"])
    s = s[s.health == "inflamed"]
    expr = s.groupby("gene").cp10k_sum.sum() / s.drop_duplicates("cluster").n_cells.sum()
    pct = expr.rank(pct=True)
    q = pct.index.isin(p.gene[p.padj.notna()])
    return {"genes_compared": int(len(pct)), "quantified": int(q.sum()),
            "median_rna_percentile_quantified": round(float(pct[q].median()), 2),
            "median_rna_percentile_not_quantified": round(float(pct[~q].median()), 2)}


def niche_markers(sc: "_SingleCell") -> dict:
    """Calibration of the niche threshold: share of fibroblast RNA in the crypt-niche clusters for the markers."""
    return {k: {g: round(float(sc.niche.get(g, np.nan)), 3) for g in v} for k, v in NICHE_MARKERS.items()}


def run() -> dict:
    from .evaluate import clinical_targets
    d = pd.read_parquet(PROCESSED / "discovery_scored.parquet")
    t = _load()
    out = _rank([_row(x, t) for x in d[d.status == "first-in-class opportunity"].itertuples()], stability=True)
    out.to_parquet(PROCESSED / "validation_shortlist.parquet", index=False)
    cal = _rank([_row(x, t) for x in d[d.gene.isin(set(clinical_targets()["gene"]))].itertuples()])
    cal.to_parquet(PROCESSED / "validation_shortlist_calibration.parquet", index=False)

    def greens(df):
        return {k: round(float(df.criteria.map(lambda s: json.loads(s)[k][0] == "green").mean()), 2) for k in CRITERIA[:-1]}
    s = {"genes": len(out), "groups": {g: out.gene[out.group == g].tolist() for g in "ABC"},
         "schemes": list(SCHEMES), "stability": dict(zip(out.gene, out.stability.astype(int))),
         "protein_quantified_in_shortlist": int((out.criteria.map(lambda c: json.loads(c)["protein"][0]) != "grey").sum()),
         "protein_bias": protein_bias(), "no_red": int((out.reds == 0).sum()),
         "regeneration_flag": out.gene[out.regeneration_flag].tolist(),
         "niche_flag": out.gene[out.niche_flag].tolist(), "niche_threshold": NICHE_SHARE,
         "neutrophil_directed": out.gene[out.neutrophil_directed].tolist(),
         "neutrophil_directed_clinical": cal.gene[cal.neutrophil_directed].tolist(), "neutrophil_made": NEUTROPHIL_MADE,
         "atlases_agree": {str(k): int(v) for k, v in out.atlases_agree.dropna().astype(int).value_counts().sort_index().items()},
         "niche_markers": niche_markers(t["sc"]),
         "niche_stromal": {r.gene: {"group": r.group, "niche_share": r.niche_share,
                                    "niche_per_cell_log2fc": r.niche_per_cell_log2fc}
                           for r in out[out.cell_family == "stromal"].itertuples()},
         "mendelian_flag": out.gene[out.criteria.str.contains("loss causes")].tolist(),
         "green_share_first_in_class": greens(out), "green_share_clinical_targets": greens(cal),
         "clinical_targets": len(cal), "median_points_first_in_class": float(out.score.median()),
         "median_points_clinical_targets": float(cal.score.median())}
    (PROCESSED / "validation_shortlist_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, ensure_ascii=False))
    for f in ("validation_shortlist.parquet", "validation_shortlist_calibration.parquet"):
        d = pd.read_parquet(PROCESSED / f)
        print("\n", f)
        for r in d.itertuples():
            c = json.loads(r.criteria)
            print(f"{r.rank:>2} {getattr(r, 'group', '-')} {getattr(r, 'stability', '-')}/6 {r.gene:<11}{r.lineage:<11}{r.score:>3} " +
                  " ".join(f"{k[:4]}={v[0][:2]}" for k, v in c.items() if k != "mouse"))
