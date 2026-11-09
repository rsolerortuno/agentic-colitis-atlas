"""External enrichment: Open Targets associations for the disease; target details (tractability, drugs, Reactome),
STRING and UniProt for a candidate set. Cached under data/cache/, summarised in data/processed/external_*."""
from __future__ import annotations

import json

import pandas as pd

from .external import ot_disease_associations, ot_target, string_partners, uniprot
from .paths import PROCESSED, RAW

N_EXPRESSION = 300   # strongest replicated meta-analysis changes
N_OT_ALWAYS = 100    # top Open Targets targets of the disease (druggability of known disease genes)


def candidate_genes(ot: pd.DataFrame) -> list[str]:
    m = pd.read_parquet(PROCESSED / "meta.parquet")
    m = m[(m["meta_padj"] < 0.05) & ((m["n_sig_up"] >= 2) | (m["n_sig_down"] >= 2))]
    m = m.reindex(m["meta_log2fc"].abs().sort_values(ascending=False).index)
    return list(dict.fromkeys(m["gene"].head(N_EXPRESSION).tolist() + ot["symbol"].head(N_OT_ALWAYS).tolist()))


def run() -> dict:
    ot = pd.DataFrame(ot_disease_associations())
    ot = ot.sort_values("score", ascending=False).drop_duplicates("symbol")
    ot.to_parquet(PROCESSED / "external_ot_associations.parquet", index=False)
    chosen = candidate_genes(ot)
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str)
    ens = hg.dropna().set_index("symbol")["ensembl_gene_id"].to_dict()
    ens.update(ot.set_index("symbol")["ensembl"].to_dict())
    details = {}
    for i, g in enumerate(chosen):
        details[g] = {"gene": g, "ensembl": ens.get(g),
                      "open_targets": ot_target(ens[g]) if ens.get(g) else None,
                      "string": string_partners(g), "uniprot": uniprot(g)}
        if i % 50 == 0:
            print(f"{i}/{len(chosen)}", flush=True)
    (PROCESSED / "external_targets.json").write_text(json.dumps(details, indent=1))
    return {"ot_associations": len(ot), "enriched_genes": len(details),
            "with_drugs": sum(1 for d in details.values() if d["open_targets"] and d["open_targets"]["drugs"])}


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
