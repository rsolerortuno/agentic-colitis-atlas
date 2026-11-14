"""UC GWAS against the three QTL layers: how much genetic support each layer gives (added 2026-11-14).

The Mendelian randomisation modules (mr.py, mr_colon.py) report genes tested and genes significant. This module adds
what lies before and beside that, for each layer (plasma pQTL from UKB-PPP, colon-tissue eQTL from GTEx v10 and CEDAR,
gut cell-type sc-eQTL from IBDverse), with the same instruments those modules use (Open Targets 26.09 fine-mapped cis
credible sets, lead variant p < 5e-8, MHC out):

1. QTL signals: credible sets that qualify as instruments, and their genes. SuSiE credible sets of one gene in one
   study are independent signals by construction; the same signal seen in several tissues or cell types is not.
2. Matched to the UC GWAS (de Lange 2017): signals whose lead variant, or another variant of the same 95 % credible set,
   is a non-palindromic SNV present in the GWAS with matching alleles. These are the MR instruments.
3. Instruments at a genome-wide UC variant: the matched variant itself has UC p < 5e-8 in de Lange 2017.
4. MR: genes significant (q < 0.05, mr.py and mr_colon.py) and, of those, the ones that colocalise.
5. From the GWAS side, for the 241 UC risk signals with an Open Targets credible set (genetic_targets.signals(); the 22
   preprint-only Lyu 2025 loci have none): does a QTL instrument of the layer share a 95 % credible-set variant with
   the signal, and does it colocalise with it (Open Targets coloc H4 >= 0.8 or eCAVIAR CLPP >= 0.01, the rule of
   mr.py)? Sharing a variant is not colocalisation: two signals can share a variant by chance in a dense region, or
   share the causal variant with different lead variants because of linkage disequilibrium. Colocalisation is
   queried for the credible set that represents each signal (others at the same locus are not queried).

The overlap between layers is counted by gene (MR significant; MR significant and colocalised) and by UC risk signal
(colocalised with at least one gene of the layer).

Outputs: data/processed/qtl_overlap.json (summary), qtl_overlap_signals.parquet (UC signal x layer x QTL gene, every
shared or colocalised pair). Colocalisation answers are cached in data/raw/extra/ot_coloc_uc_signals_cache.json.
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
import requests

from . import mr, mr_colon
from .paths import PROCESSED, RAW

LAYERS = {"plasma_protein": "Plasma protein (UKB-PPP pQTL)", "colon_tissue": "Colon tissue (GTEx v10, CEDAR eQTL)",
          "gut_cell": "Gut cell types (IBDverse sc-eQTL)"}
H4, CLPP, P_GW = 0.8, 0.01, 5e-8
CACHE = RAW / "extra" / "ot_coloc_uc_signals_cache.json"
QUERY = """query($id:String!,$i:Int!){ credibleSet(studyLocusId:$id){ colocalisation(studyTypes:[pqtl,eqtl,sceqtl],
  page:{index:$i,size:500}){ count rows{ h4 clpp colocalisationMethod otherStudyLocus{ studyLocusId } } } } }"""


# ------------------------------------------------------------------ instruments of each layer
def layer_instruments() -> dict[str, pd.DataFrame]:
    """Instrument credible sets per layer (before matching to the GWAS), with gene and layer."""
    p = mr.instruments().assign(layer="plasma_protein")
    g = mr_colon.instruments()
    g["layer"] = g["studyId"].map(mr_colon.studies().set_index("studyId")["layer"])
    return {"plasma_protein": p, "colon_tissue": g[g.layer == "colon_tissue"], "gut_cell": g[g.layer == "gut_cell"]}


def matched_instruments() -> dict[str, pd.DataFrame]:
    """The harmonised instruments that mr.py and mr_colon.py saved (one per credible set)."""
    h = pd.read_parquet(PROCESSED / "mr_colon_instruments.parquet")
    lay = mr_colon.studies().set_index("studyId")["layer"]
    h["layer"] = h["studyId"].map(lay)
    return {"plasma_protein": pd.read_parquet(PROCESSED / "mr_instruments.parquet"),
            "colon_tissue": h[h.layer == "colon_tissue"], "gut_cell": h[h.layer == "gut_cell"]}


def _lead_usable(v: str) -> bool:
    c, p, ref, alt = v.split("_")
    return len(ref) == 1 and len(alt) == 1 and not mr._palindromic(ref, alt)


# ------------------------------------------------------------------ colocalisation from the GWAS side
def _query(sl: str) -> list[dict]:
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    if sl in cache:
        return cache[sl]
    rows, i = [], 0
    while True:
        for attempt in range(5):
            try:
                r = requests.post(mr.API, json={"query": QUERY, "variables": {"id": sl, "i": i}}, timeout=120)
                if r.status_code == 200 and r.text.strip():
                    data = r.json()
                    break
            except (requests.RequestException, ValueError):
                pass
            time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f"Open Targets colocalisation query failed for {sl}")
        col = ((data.get("data") or {}).get("credibleSet") or {}).get("colocalisation") or {}
        page = col.get("rows") or []
        rows += [{"other": x["otherStudyLocus"]["studyLocusId"], "h4": x.get("h4"), "clpp": x.get("clpp"),
                  "method": x.get("colocalisationMethod")} for x in page if x.get("otherStudyLocus")]
        i += 1
        if len(rows) >= (col.get("count") or 0) or not page:
            break
    cache[sl] = rows
    CACHE.write_text(json.dumps(cache), encoding="utf-8")
    return rows


def signal_table(inst: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, int]:
    """UC signal x layer x QTL gene: shared 95 % credible-set variant, colocalisation (best H4 and CLPP)."""
    sig = pd.read_parquet(PROCESSED / "genetic_signals.parquet")
    sig = sig[~sig.preprint_only]
    sig_vars = sig.groupby(["signal", "cs"])["variant"].apply(set)
    by_cs = {}
    for layer, d in inst.items():
        for r in d.itertuples():
            by_cs[r.studyLocusId] = (layer, r.gene, r.studyId)
    var_index: dict[str, list[str]] = {}
    for layer, d in inst.items():
        for r in d.itertuples():
            for v in {r.variantId, *[x[0] for x in json.loads(r.locus)]}:
                var_index.setdefault(v, []).append(r.studyLocusId)
    rows = {}
    for (signal, cs), vs in sig_vars.items():
        for v in vs:
            for sl in var_index.get(v, []):
                layer, gene, study = by_cs[sl]
                k = (signal, layer, gene)
                rows.setdefault(k, {"signal": signal, "uc_credible_set": cs, "layer": layer, "gene": gene,
                                    "shares_variant": True, "h4": np.nan, "clpp": np.nan, "coloc_study": None})
        for x in _query(cs):
            if x["other"] not in by_cs:
                continue
            layer, gene, study = by_cs[x["other"]]
            k = (signal, layer, gene)
            r = rows.setdefault(k, {"signal": signal, "uc_credible_set": cs, "layer": layer, "gene": gene,
                                    "shares_variant": False, "h4": np.nan, "clpp": np.nan, "coloc_study": None})
            better = (x["h4"] or 0) > (0 if np.isnan(r["h4"]) else r["h4"]) or \
                     (x["clpp"] or 0) > (0 if np.isnan(r["clpp"]) else r["clpp"])
            if x["h4"] is not None:
                r["h4"] = max(x["h4"], 0 if np.isnan(r["h4"]) else r["h4"])
            if x["clpp"] is not None:
                r["clpp"] = max(x["clpp"], 0 if np.isnan(r["clpp"]) else r["clpp"])
            if better:
                r["coloc_study"] = study
    t = pd.DataFrame(list(rows.values()), columns=["signal", "uc_credible_set", "layer", "gene", "shares_variant",
                                                   "h4", "clpp", "coloc_study"])
    t["colocalises"] = (t.h4.fillna(0) >= H4) | (t.clpp.fillna(0) >= CLPP)
    return t.sort_values(["signal", "layer", "gene"]).reset_index(drop=True), sig.signal.nunique()


# ------------------------------------------------------------------ summary
def _venn(sets: dict[str, set]) -> dict:
    a, b, c = (sets[k] for k in LAYERS)
    return {"plasma_only": len(a - b - c), "colon_only": len(b - a - c), "cell_only": len(c - a - b),
            "plasma_and_colon_only": len((a & b) - c), "plasma_and_cell_only": len((a & c) - b),
            "colon_and_cell_only": len((b & c) - a), "all_three": len(a & b & c), "any": len(a | b | c),
            "all_three_genes": sorted(a & b & c)}


def run() -> dict:
    inst = layer_instruments()
    matched = matched_instruments()
    mrp = pd.read_parquet(PROCESSED / "mr_pqtl.parquet").assign(layer="plasma_protein")
    mrc = pd.read_parquet(PROCESSED / "mr_colon.parquet")
    mrr = pd.concat([mrp, mrc], ignore_index=True)
    sig_t, n_signals = signal_table(inst)
    layers, mr_sig, mr_coloc, sig_coloc, sig_share = {}, {}, {}, {}, {}
    for layer in LAYERS:
        d, h = inst[layer], matched[layer]
        r = mrr[mrr.layer == layer]
        s = r[r.q < 0.05]
        sc = s[s.colocalises.fillna(False).astype(bool)]
        st = sig_t[sig_t.layer == layer]
        mr_sig[layer], mr_coloc[layer] = set(s.gene), set(sc.gene)
        sig_coloc[layer], sig_share[layer] = set(st[st.colocalises].signal), set(st[st.shares_variant].signal)
        coloc_genes = set(st[st.colocalises].gene)
        layers[layer] = {
            "label": LAYERS[layer],
            "qtl_signals": len(d), "qtl_genes": int(d.gene.nunique()),
            "lead_variants": int(d.variantId.nunique()),
            "lead_variants_snv_not_palindromic": int(pd.Series(d.variantId.unique()).map(_lead_usable).sum()),
            "matched_instruments": len(h), "matched_genes": int(h.gene.nunique()),
            "matched_via_tag_variant": int(h.substitute_pp.notna().sum()),
            "matched_gwas_variants": int(h.gwas_variant.nunique()),
            "matched_share": round(len(h) / len(d), 3) if len(d) else None,
            "genes_with_2plus_instruments": int((r.n_instruments >= 2).sum()),
            "instruments_at_gw_uc_variant": int((h.p_outcome < P_GW).sum()),
            "genes_with_gw_uc_instrument": int(h.loc[h.p_outcome < P_GW, "gene"].nunique()),
            "mr_genes_tested": len(r), "mr_significant": len(s), "mr_significant_colocalised": len(sc),
            "uc_signals_sharing_a_variant": len(sig_share[layer]),
            "uc_signals_colocalised": len(sig_coloc[layer]),
            "colocalised_genes": len(coloc_genes),
            "colocalised_genes_mr_significant": len(coloc_genes & mr_sig[layer]),
        }
    shared_any = set().union(*sig_share.values())
    coloc_any = set().union(*sig_coloc.values())
    out = {"uc_signals": int(n_signals), "gwas": "de Lange 2017 (GCST004133) for matching and MR",
           "rule": f"colocalised when coloc H4 >= {H4} or eCAVIAR CLPP >= {CLPP}", "layers": layers,
           "overlap_mr_significant": _venn(mr_sig), "overlap_mr_significant_colocalised": _venn(mr_coloc),
           "overlap_uc_signals_colocalised": _venn(sig_coloc),
           "uc_signals_sharing_a_variant_any_layer": len(shared_any),
           "uc_signals_colocalised_any_layer": len(coloc_any),
           "uc_signals_sharing_but_not_colocalised": len(shared_any - coloc_any)}
    sig_t.to_parquet(PROCESSED / "qtl_overlap_signals.parquet", index=False)
    (PROCESSED / "qtl_overlap.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
