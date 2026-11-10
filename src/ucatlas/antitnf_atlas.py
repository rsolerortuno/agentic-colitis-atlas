"""Anti-TNF single-cell atlas (Thomas et al. 2024, Nat Immunol 25:2152; GEO GSE282122; added 2026-11-10, fourth round).

About one million cells from 216 gut biopsies of 41 people (Crohn's disease, UC, healthy), before and after
adalimumab, annotated into 109 cell states by the authors. GEO holds the count matrices without cell annotations, so
re-annotating them would not reproduce the published cell states; the analysis therefore uses the authors' published
results (supplementary tables, Europe PMC package of PMC11519010):
- Supplementary Table 1A iii: UC patients and their remission status after adalimumab;
- Supplementary Table 8 (Pre_tx_UC_cell_state, Pre_tx_UC_compartment): abundance of each cell state at baseline,
  non-remission vs remission (MASC odds ratios, BH);
- Supplementary Table 9 (Baseline_UC): genes differing at baseline between non-remission and remission within each
  cell state (MAST with age, sex, site, disease duration and nested random effects; positive coef = higher in
  non-remission);
- Supplementary Table 10 (Longitudinal_UC_cellstate): genes changing after treatment within each cell state, in
  paired samples (positive coef = higher after treatment).
Outputs: data/processed/antitnf_baseline_genes.parquet, antitnf_abundance.parquet, antitnf_longitudinal.parquet,
antitnf_summary.json.
"""
from __future__ import annotations

import io
import json
import warnings
import zipfile

import pandas as pd

from .paths import PROCESSED, RAW

SRC = RAW / "extra" / "thomas2024_supplement.zip"
STUDY = "Thomas et al. 2024, Nat Immunol (GSE282122)"


def _sheet(name: str, sheet: str) -> pd.DataFrame:
    with zipfile.ZipFile(SRC) as z, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pd.read_excel(io.BytesIO(z.read(name)), sheet)


def run() -> dict:
    uc = _sheet("41590_2024_1994_MOESM3_ESM.xlsx", "Supp_Table_1A_iii_UCspecific")
    rem = uc.drop_duplicates("Patient").Remission.value_counts().to_dict()
    base = _sheet("41590_2024_1994_MOESM9_ESM.xlsx", "Baseline_UC")
    base = base.rename(columns={"adjusted_pval": "padj", "coef": "coef_nonremission"})
    # Excel turned a few gene names (MARCH*, SEPT*) into dates in the published table; they are dropped
    base = base[base.gene.map(lambda v: isinstance(v, str))]
    base[["cell_state", "gene", "coef_nonremission", "pvalue", "padj"]].to_parquet(PROCESSED / "antitnf_baseline_genes.parquet", index=False)
    ab = _sheet("41590_2024_1994_MOESM8_ESM.xlsx", "Pre_tx_UC_cell_state")
    ab = ab.rename(columns={"model.pvalue": "p", "RemissionNon_Remission.OR": "or_nonremission", "BH": "padj"})
    ab[["cell_state", "size", "or_nonremission", "p", "padj"] + (["category"] if "category" in ab else [])].to_parquet(
        PROCESSED / "antitnf_abundance.parquet", index=False)
    lo = _sheet("41590_2024_1994_MOESM10_ESM.xlsx", "Longitudinal_UC_cellstate")
    padj_col = "fdr" if "fdr" in lo else "adjusted_pval"
    lo = lo.rename(columns={padj_col: "padj", "coef": "coef_after"})
    lo = lo[lo.gene.map(lambda v: isinstance(v, str))]
    lo[["cell_state", "gene", "coef_after", "pvalue", "padj"]].to_parquet(PROCESSED / "antitnf_longitudinal.parquet", index=False)
    s = {"study": STUDY, "uc_patients": int(sum(rem.values())), "uc_remission": int(rem.get("Remission", 0)),
         "uc_non_remission": int(rem.get("Non_Remission", 0)),
         "baseline_genes": int(base.gene.nunique()), "baseline_cell_states": int(base.cell_state.nunique()),
         "baseline_pairs_padj05": int((base.padj < 0.05).sum()),
         "abundance_cell_states_padj05": ab[ab.padj < 0.05].sort_values("p")[["cell_state", "or_nonremission"]]
         .round(3).to_dict("records"),
         "longitudinal_pairs_padj05": int((lo.padj < 0.05).sum())}
    (PROCESSED / "antitnf_summary.json").write_text(json.dumps(s, indent=1), encoding="utf-8")
    return s


def gene_view(gene: str) -> dict:
    b = pd.read_parquet(PROCESSED / "antitnf_baseline_genes.parquet")
    lo = pd.read_parquet(PROCESSED / "antitnf_longitudinal.parquet")
    bb = b[(b.gene == gene) & (b.padj < 0.05)].sort_values("padj")
    ll = lo[(lo.gene == gene) & (lo.padj < 0.05)].sort_values("padj")
    return {"baseline_nonremission_vs_remission": [{"cell_state": r.cell_state, "coef": round(float(r.coef_nonremission), 3),
                                                     "padj": float(f"{r.padj:.3g}")} for r in bb.head(8).itertuples()],
            "after_vs_before_treatment": [{"cell_state": r.cell_state, "coef": round(float(r.coef_after), 3),
                                           "padj": float(f"{r.padj:.3g}")} for r in ll.head(8).itertuples()]}


def load(store=None) -> int:
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    b = pd.read_parquet(PROCESSED / "antitnf_baseline_genes.parquet")
    lo = pd.read_parquet(PROCESSED / "antitnf_longitudinal.parquet")
    s = json.loads((PROCESSED / "antitnf_summary.json").read_text())
    genes = set(b[b.padj < 0.05].gene) | set(lo[lo.padj < 0.05].gene)
    n = 0
    with store.batch():
        store.add_node("DATASET:GSE282122", "Dataset", f"Anti-TNF (adalimumab) longitudinal single-cell atlas of IBD ({STUDY})")
        for g in sorted(genes):
            node = store.resolve(g)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": g, "study": STUDY, "disease": "ulcerative colitis", "treatment": "adalimumab",
                   "patients": f"{s['uc_remission']} remission, {s['uc_non_remission']} non-remission",
                   "direction_note": "baseline coef > 0: higher in non-remission; after-treatment coef > 0: higher after",
                   **gene_view(g)}
            ev = store.add_evidence(evidence_id("antitnf_singlecell", g), "antitnf_singlecell", "DATASET:GSE282122",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "associated_with_response_in", "DATASET:GSE282122", ev)
            n += 1
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    for g in ["DUOX2", "LCN2", "S100A8", "CXCL8", "OSM", "OSMR", "IL11", "CHI3L1", "MMP1", "IDO1", "TREM1", "CD274", "CXCR2"]:
        print(g, gene_view(g))
