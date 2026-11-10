"""M4: build the knowledge graph from expression results (M2), literature (M3) and the external cache."""
from __future__ import annotations

import json
import re
from datetime import date

import pandas as pd

from .config import load
from .graph import GraphStore, evidence_id
from .literature_entities import is_disease, resolve_gene
from .paths import LITERATURE, PROCESSED, RAW

DISEASE = "MONDO:0005101"
TODAY = date.today().isoformat()
DATASETS = {
    "GSE193677": "Mount Sinai Crohn's and Colitis Registry (MSCCR) biopsies, RNA-seq (discovery)",
    "GSE206285": "UNIFI ustekinumab trial, baseline colon biopsies, Affymetrix",
    "GSE92415": "PURSUIT-SC golimumab trial, colon biopsies, Affymetrix",
    "GSE16879": "Infliximab cohort (Leuven), colon biopsies, Affymetrix",
    "GSE73661": "Infliximab and vedolizumab cohort (Leuven), colon biopsies, Affymetrix",
    "GSE109142": "PROTECT, treatment-naive paediatric UC, rectal biopsies, RNA-seq",
    "GSE117993": "RISK, paediatric UC at diagnosis, rectal biopsies, RNA-seq",
    "GSE165512": "Adult UC and healthy colon biopsies, RNA-seq",
    "META": "Random-effects meta-analysis of UC vs control across the eight cohorts",
}
CONTRAST_TEXT = {"inflamed_vs_control": "UC inflamed vs control", "noninflamed_vs_control": "UC non-inflamed vs control",
                 "inflamed_vs_noninflamed_paired": "UC inflamed vs non-inflamed, same patient",
                 "uc_vs_control": "UC vs healthy control"}


def gene_node_id(hgnc_id: str) -> str:
    return hgnc_id if hgnc_id.startswith("HGNC:") else f"HGNC:{hgnc_id}"


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")[:60]


def build(store: GraphStore | None = None) -> dict:
    cfg = load()
    store = store or GraphStore()
    store.reset()
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", dtype=str,
                     usecols=["hgnc_id", "symbol", "name", "location", "ensembl_gene_id", "uniprot_ids", "entrez_id"])
    by_symbol = hg.set_index("symbol")
    by_entrez = hg.dropna(subset=["entrez_id"]).set_index("entrez_id")["symbol"].to_dict()

    def gene(symbol: str) -> str | None:
        if symbol not in by_symbol.index:
            return None
        r = by_symbol.loc[symbol]
        nid = gene_node_id(r["hgnc_id"])
        store.add_node(nid, "Gene", symbol, hgnc_id=r["hgnc_id"], full_name=r["name"], location=r["location"],
                       ensembl=r["ensembl_gene_id"], uniprot=r["uniprot_ids"])
        return nid

    with store.batch():
        store.add_node(DISEASE, "Disease", cfg["disease"]["name"], mondo=DISEASE,
                       open_targets=cfg["disease"]["open_targets_id"])
        for code, dname in DATASETS.items():
            store.add_node(f"DATASET:{code}", "Dataset", dname, geo=code if code.startswith("GSE") else None)

        # ---------------------------------------------------------- expression (M2)
        thr = cfg["expression"]
        de = pd.read_parquet(PROCESSED / "de_results.parquet").dropna(subset=["padj"])
        for r in de.itertuples():
            gid = gene(r.gene)
            if gid is None:
                continue
            sig = r.padj < thr["padj"] and abs(r.log2fc) >= thr["min_abs_log2fc"]
            did = f"DATASET:{r.cohort}"
            ev = store.add_evidence(evidence_id("omics_de", r.cohort, r.contrast, r.gene), "omics_de", did, TODAY,
                                    layer="rna", gene=r.gene, cohort=r.cohort, contrast=r.contrast,
                                    contrast_text=CONTRAST_TEXT[r.contrast], log2fc=round(r.log2fc, 4),
                                    se=round(r.se, 4), padj=float(f"{r.padj:.3g}"), pvalue=float(f"{r.p:.3g}"),
                                    n_disease=int(r.n1), n_control=int(r.n2), significant=bool(sig))
            changed = sig
            pred = ("upregulated_in" if r.log2fc > 0 else "downregulated_in") if changed else "not_changed_in"
            store.add_claim(gid, pred, DISEASE, ev)
            store.add_claim(gid, "measured_in", did, ev)
        meta = pd.read_parquet(PROCESSED / "meta.parquet").dropna(subset=["meta_padj"])
        for r in meta.itertuples():
            gid = gene(r.gene)
            if gid is None:
                continue
            sig = r.meta_padj < thr["padj"] and abs(r.meta_log2fc) >= thr["min_abs_log2fc"]
            ev = store.add_evidence(evidence_id("expression_meta", r.gene), "expression_meta", "DATASET:META", TODAY,
                                    gene=r.gene, layer="rna", log2fc=round(r.meta_log2fc, 4), se=round(r.meta_se, 4),
                                    padj=float(f"{r.meta_padj:.3g}"), i2=round(r.i2, 3), cohorts=int(r.k),
                                    cohorts_significant_up=int(r.n_sig_up), cohorts_significant_down=int(r.n_sig_down),
                                    significant=bool(sig), method="DerSimonian-Laird random effects")
            pred = ("upregulated_in" if r.meta_log2fc > 0 else "downregulated_in") if sig else "not_changed_in"
            store.add_claim(gid, pred, DISEASE, ev)

        # ---------------------------------------------------------- literature (M3)
        pubs = pd.read_parquet(LITERATURE / "publications.parquet").set_index("pmid")

        def publication(pmid: int) -> str:
            nid = f"PMID:{pmid}"
            p = pubs.loc[pmid]
            store.add_node(nid, "Publication", p["title"][:300] or nid, pmid=int(pmid), year=p["year"],
                           journal=p["journal"][:200], url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            return nid

        def entity(text: str) -> str:
            g = resolve_gene(text)
            if g and (gid := gene(g)):
                return gid
            if is_disease(text):
                return DISEASE
            return store.add_node(f"CONCEPT:{_slug(text)}", "Concept", text.strip()[:120])

        rel = pd.read_parquet(LITERATURE / "relations.parquet")
        for r in rel.itertuples():
            pid = publication(int(r.pmid))
            s, o = entity(r.subject), entity(r.object)
            ev = store.add_evidence(evidence_id("literature", r.relation_id), "literature_relation", pid, TODAY,
                                    relation_id=r.relation_id, pmid=int(r.pmid), sentence=r.literal_sentence,
                                    confidence=float(r.confidence), subject_text=r.subject, object_text=r.object,
                                    extractor="extractor-relaciones (Claude Sonnet 5.5), verbatim-validated")
            store.add_claim(s, r.relation, o, ev)
            store.add_claim(pid, "reports", s, ev)

        ann = pd.read_parquet(LITERATURE / "annotations.parquet")
        genes_ann = ann[(ann["type"] == "Gene") & ann["identifier"].notna()].drop_duplicates(["pmid", "identifier"])
        for r in genes_ann.itertuples():
            sym = None
            for ident in str(r.identifier).split(";"):
                sym = by_entrez.get(ident.strip())
                if sym:
                    break
            if not sym or (gid := gene(sym)) is None:
                continue
            pid = publication(int(r.pmid))
            ev = store.add_evidence(evidence_id("pubtator", r.pmid, r.identifier), "pubtator_annotation", "PubTator3",
                                    TODAY, pmid=int(r.pmid), ncbi_gene=str(r.identifier), mention=r.text)
            store.add_claim(pid, "mentions", gid, ev)

        # gene mentions per publication (PubTator-normalized): mention text -> gene node id
        mentions_by_pmid: dict[int, dict[str, str]] = {}
        for r in genes_ann.itertuples():
            sym = next((by_entrez.get(i.strip()) for i in str(r.identifier).split(";") if by_entrez.get(i.strip())), None)
            if sym and sym in by_symbol.index and isinstance(r.text, str) and len(r.text) >= 2:
                mentions_by_pmid.setdefault(int(r.pmid), {})[r.text] = gene_node_id(by_symbol.loc[sym, "hgnc_id"])

        chunks = pd.read_parquet(LITERATURE / "chunks.parquet")
        for r in chunks.itertuples():
            pid = publication(int(r.pmid))
            cid = f"CHUNK:{r.chunk_id}"
            store.add_node(cid, "Chunk", r.text[:120], text=r.text, pmid=int(r.pmid))
            ev = store.add_evidence(evidence_id("chunk", r.chunk_id), "provenance", "PubTator3 abstract", TODAY,
                                    pmid=int(r.pmid))
            store.add_claim(cid, "part_of", pid, ev)
            for mention, gid in mentions_by_pmid.get(int(r.pmid), {}).items():
                if re.search(rf"(?<![A-Za-z0-9]){re.escape(mention)}(?![A-Za-z0-9])", r.text):
                    store.add_claim(cid, "mentions", gid, ev)

        # ---------------------------------------------------------- external (M5 cache)
        ot = pd.read_parquet(PROCESSED / "external_ot_associations.parquet")
        for r in ot.itertuples():
            gid = gene(r.symbol)
            if gid is None:
                continue
            ev = store.add_evidence(evidence_id("open_targets_assoc", r.ensembl), "open_targets_association",
                                    "Open Targets Platform", TODAY, score=round(r.score, 4), ensembl=r.ensembl,
                                    gene=r.symbol, genetic_association=round(r.dt_genetic_association, 4),
                                    literature=round(r.dt_literature, 4), known_drug=round(r.dt_known_drug, 4),
                                    disease=cfg["disease"]["open_targets_id"],
                                    url=f"https://platform.opentargets.org/evidence/{r.ensembl}/{cfg['disease']['open_targets_id']}")
            store.add_claim(gid, "associated_with", DISEASE, ev)

        details = json.loads((PROCESSED / "external_targets.json").read_text())
        for sym, d in details.items():
            gid = gene(sym)
            if gid is None:
                continue
            t = d.get("open_targets")
            if t:
                ev_t = store.add_evidence(evidence_id("open_targets_tract", sym), "open_targets_tractability",
                                          "Open Targets Platform", TODAY, gene=sym,
                                          tractable_modalities=t["tractable_modalities"],
                                          has_approved_or_clinical_drug=t["has_approved_or_clinical_drug"],
                                          labels=t["tractability_labels"][:20])
                store.add_claim(gid, "tractability_assessed_in", DISEASE, ev_t)
                for drug in t["drugs"]:
                    did = f"CHEMBL:{drug['chembl']}"
                    store.add_node(did, "Drug", drug["name"], chembl=drug["chembl"])
                    ev = store.add_evidence(evidence_id("open_targets_drug", drug["chembl"], sym), "open_targets_drug",
                                            "Open Targets Platform", TODAY, gene=sym, drug=drug["name"],
                                            max_clinical_stage=drug["max_stage"],
                                            url=f"https://platform.opentargets.org/drug/{drug['chembl']}")
                    store.add_claim(did, "targets", gid, ev)
                for pw in t["pathways"]:
                    pid = f"REACTOME:{pw['reactome']}"
                    store.add_node(pid, "Pathway", pw["name"], reactome=pw["reactome"], top_level=pw["top_level"])
                    ev = store.add_evidence(evidence_id("reactome", pw["reactome"], sym), "pathway_membership",
                                            "Reactome via Open Targets", TODAY, gene=sym, reactome=pw["reactome"])
                    store.add_claim(gid, "in_pathway", pid, ev)
            for p in d.get("string") or []:
                pg = gene(p["partner"])
                if pg is None:
                    continue
                ev = store.add_evidence(evidence_id("string", *sorted([sym, p["partner"]])), "string_interaction",
                                        "STRING", TODAY, a=sym, b=p["partner"], score=p["score"],
                                        experimental=p["experimental"])
                store.add_claim(gid, "interacts_with", pg, ev)
            u = d.get("uniprot")
            if u:
                prot = f"UNIPROT:{u['accession']}"
                store.add_node(prot, "Protein", u["name"] or u["accession"], accession=u["accession"],
                               function=u["function"], subcellular_location=u["subcellular_location"])
                ev = store.add_evidence(evidence_id("uniprot", u["accession"]), "uniprot_entry", "UniProtKB/Swiss-Prot",
                                        TODAY, accession=u["accession"],
                                        url=f"https://www.uniprot.org/uniprotkb/{u['accession']}")
                store.add_claim(gid, "encodes", prot, ev)

    stats = store.stats()
    (PROCESSED / "graph_stats.json").write_text(json.dumps(stats, indent=1))
    return stats


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
