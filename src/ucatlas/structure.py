"""Protein structure and pockets of the candidates (added 2026-11-10; ported from the kidney atlas).

For every report candidate, the top of the convergence ranking, the leading first-in-class genes and a set of UC
small-molecule targets used as controls: the AlphaFold DB model, residues below pLDDT MIN_PLDDT removed, cavities
detected with pyKVFinder, and a pocket score that orders the cavities of one protein:

    pocket score = 0.5 x confidence + 0.4 x geometry + 0.1 x hydrophobicity
    confidence    = clip((mean pLDDT of the wall residues - 50) / 40, 0, 1)
    geometry      = fit of the volume to 150-1,500 A^3 (1 inside, log decay outside) x clip(max depth / 10 A, 0, 1)
    hydrophobicity = clip((average Eisenberg hydropathy of the cavity + 1) / 2, 0, 1)

Weights fixed before looking at results (evals/history.md). The kidney atlas also scored lysines and exposed loops
for molecular-glue degraders; that question is not asked here, so those components are not used. The score is a
heuristic about the shape of a predicted monomer; it does not predict binding.
Check: residues in contact with ligands in experimental structures (PDBe-KB) - are they covered by the detected
cavities (detector), and by the top three pockets (score)? Controls: UC small-molecule targets with approved or late
drugs (JAK1, JAK2, JAK3, TYK2, S1PR1, S1PR5, PDE4B, PTGS2, PPARG, DHODH), whose pockets are real.
Outputs: data/processed/structure_pockets.parquet, structure_summary.json; models in data/cache/alphafold.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from .paths import CACHE, PROCESSED, RAW

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/"
PDBE_LIGAND_SITES = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/"
AF_DIR = CACHE / "alphafold"
UA = {"User-Agent": "uc-graphrag/0.1 (research; multi-omics GraphRAG case study)"}
MIN_PLDDT, CONFIDENT_PLDDT, POCKET_VOLUME, TOP_POCKETS = 50, 70, (150, 1500), 3
WEIGHTS = {"confidence": 0.5, "geometry": 0.4, "hydrophobicity": 0.1}
CONTROLS = ["JAK1", "JAK2", "JAK3", "TYK2", "S1PR1", "S1PR5", "PDE4B", "PTGS2", "PPARG", "DHODH"]


def _get_json(url: str, cache_name: str):
    AF_DIR.mkdir(parents=True, exist_ok=True)
    path = AF_DIR / cache_name
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    import time
    for attempt in range(4):   # PDBe answers 504 under load; a server error is retried, not cached
        try:
            r = requests.get(url, timeout=120, headers=UA)
        except requests.RequestException:
            time.sleep(5 * (attempt + 1))
            continue
        if r.status_code == 404:
            path.write_text("null", encoding="utf-8")
            return None
        if r.status_code >= 500:
            time.sleep(5 * (attempt + 1))
            continue
        r.raise_for_status()
        path.write_text(r.text, encoding="utf-8")
        return r.json()
    return None


def alphafold_model(accession: str) -> dict | None:
    meta = _get_json(AF_API + accession, f"{accession}.api.json")
    if not meta:
        return None
    m = meta[0]
    pdb = AF_DIR / f"{accession}.pdb"
    if not pdb.exists():
        r = requests.get(m["pdbUrl"], timeout=300, headers=UA)
        r.raise_for_status()
        pdb.write_bytes(r.content)
    return {"accession": accession, "path": pdb, "model_entity": m["modelEntityId"], "version": m["latestVersion"],
            "mean_plddt": m["globalMetricValue"], "created": m["modelCreatedDate"][:10],
            "fraction_plddt_confident": round(m["fractionPlddtConfident"] + m["fractionPlddtVeryHigh"], 3),
            "length": m["uniprotEnd"], "url": f"https://alphafold.ebi.ac.uk/entry/{accession}", "model_url": m["pdbUrl"]}


def read_structure(path: Path) -> dict:
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("af", str(path))
    chain = next(iter(next(iter(s))))
    res = [r for r in chain if r.id[0] == " " and "CA" in r]
    return {"numbers": np.array([r.id[1] for r in res]), "names": [r.get_resname() for r in res],
            "ca": np.array([r["CA"].coord for r in res], dtype=float),
            "plddt": np.array([float(np.mean([a.bfactor for a in r])) for r in res])}


def trimmed_pdb(path: Path, st: dict, out: Path) -> int:
    keep = {int(nr) for nr, p in zip(st["numbers"], st["plddt"]) if p >= MIN_PLDDT}
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines()
             if (ln.startswith(("ATOM", "HETATM")) and int(ln[22:26]) in keep)
             or ln.startswith(("HEADER", "TITLE", "SEQRES", "TER", "END"))]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(keep)


def _fit(x: float, lo: float, hi: float) -> float:
    if x <= 0:
        return 0.0
    if lo <= x <= hi:
        return 1.0
    span = math.log(lo / max(x, 1e-6)) if x < lo else math.log(x / hi)
    return float(max(0.0, 1 - span / math.log(3)))


def cavity_table(pdb: Path, st: dict) -> list[dict]:
    import pyKVFinder
    r = pyKVFinder.run_workflow(str(pdb), include_depth=True, include_hydropathy=True)
    if not r or not r.ncav:
        return []
    by_number = {int(n): i for i, n in enumerate(st["numbers"])}
    out = []
    for cav, vol in r.volume.items():
        nums = sorted({int(x[0]) for x in r.residues[cav] if int(x[0]) in by_number})
        if not nums:
            continue
        idx = np.array([by_number[n] for n in nums])
        comp = {"confidence": float(np.clip((st["plddt"][idx].mean() - 50) / 40, 0, 1)),
                "geometry": _fit(vol, *POCKET_VOLUME) * float(np.clip(r.max_depth[cav] / 10, 0, 1)),
                "hydrophobicity": float(np.clip((r.avg_hydropathy[cav] + 1) / 2, 0, 1))}
        out.append({"cavity": cav, "volume_a3": round(float(vol), 1), "area_a2": round(float(r.area[cav]), 1),
                    "max_depth_a": round(float(r.max_depth[cav]), 2), "avg_hydropathy": round(float(r.avg_hydropathy[cav]), 3),
                    "n_residues": len(nums), "residues": nums, "mean_plddt": round(float(st["plddt"][idx].mean()), 1),
                    "components": {k: round(v, 3) for k, v in comp.items()},
                    "pocket_score": round(sum(WEIGHTS[k] * comp[k] for k in WEIGHTS), 4)})
    out.sort(key=lambda c: -c["pocket_score"])
    for i, c in enumerate(out, 1):
        c["rank"] = i
    return out


def ligand_sites(accession: str) -> dict:
    data = _get_json(PDBE_LIGAND_SITES + accession, f"{accession}.ligand_sites.json")
    if not data or accession not in data:
        return {"ligands": 0, "residues": []}
    res, ligs = set(), []
    for lig in data[accession]["data"]:
        ligs.append(lig["accession"])
        for r in lig["residues"]:
            res.update(range(int(r["startIndex"]), int(r["endIndex"]) + 1))
    return {"ligands": len(ligs), "ligand_ids": sorted(set(ligs))[:10], "residues": sorted(res)}


def validate(cavities: list[dict], sites: dict, top: int) -> dict:
    if not sites["residues"] or not cavities:
        return {"checked": False}
    known = set(sites["residues"])
    overlaps = {c["rank"]: len(known & set(c["residues"])) for c in cavities}
    best = min((r for r, n in overlaps.items() if n), default=None)
    covered = set().union(*[set(c["residues"]) for c in cavities]) & known
    covered_top = set().union(*[set(c["residues"]) for c in cavities[:top]]) & known
    return {"checked": True, "n_ligands": sites["ligands"], "n_known_residues": len(known), "n_cavities": len(cavities),
            "recall": round(len(covered) / len(known), 3), "recall_in_top_pockets": round(len(covered_top) / len(known), 3),
            "best_cavity_rank": best}


def _uniprot() -> dict[str, str]:
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "uniprot_ids"], dtype=str).dropna()
    return {s: u.split("|")[0].strip() for s, u in zip(hg["symbol"], hg["uniprot_ids"])}


def targets(run: Path | None = None) -> list[str]:
    """Report candidates, the top 25 of the convergence ranking and the first 10 first-in-class genes."""
    from .runs import current_run
    out = []
    run = run or current_run()
    if run and (run / "findings.json").exists():
        f = json.loads((run / "findings.json").read_text(encoding="utf-8"))
        out += [t["gene"] for t in f.get("targets", [])]
        out += [r["gene"] for r in f["convergence"]["ranking"][:25]]
    d = pd.read_parquet(PROCESSED / "discovery.parquet")
    out += d[d.status == "first-in-class opportunity"].sort_values("meta_log2fc", ascending=False)["gene"].head(10).tolist()
    return list(dict.fromkeys(out))


def run(genes: list[str] | None = None) -> dict:
    genes = genes or targets()
    controls = [g for g in CONTROLS if g not in genes]
    uni = _uniprot()
    results, rows = [], []
    for gene in genes + controls:
        acc = uni.get(gene)
        if not acc:
            results.append({"gene": gene, "modeled": False, "reason": "no UniProt accession"})
            continue
        model = alphafold_model(acc)
        if not model:
            results.append({"gene": gene, "modeled": False, "reason": "no AlphaFold DB entry", "uniprot": acc})
            continue
        st = read_structure(model["path"])
        trimmed = AF_DIR / f"{acc}.plddt{MIN_PLDDT}.pdb"
        kept = trimmed_pdb(model["path"], st, trimmed)
        all_cav = cavity_table(trimmed, st) if kept >= 30 else []
        cav = all_cav[:TOP_POCKETS]
        entry = {"gene": gene, "modeled": True, "control": gene in controls, "uniprot": acc,
                 **{k: v for k, v in model.items() if k != "path"}, "residues_kept_above_plddt": kept,
                 "fraction_confident_residues": round(float((st["plddt"] >= CONFIDENT_PLDDT).mean()), 3),
                 "pockets": cav, "n_cavities": len(all_cav),
                 "ligand_site_validation": validate(all_cav, ligand_sites(acc), TOP_POCKETS),
                 "best_pocket_score": cav[0]["pocket_score"] if cav else None}
        results.append(entry)
        for c in cav:
            rows.append({"gene": gene, "control": gene in controls, "uniprot": acc, "rank": c["rank"],
                         "cavity": c["cavity"], "volume_a3": c["volume_a3"], "area_a2": c["area_a2"],
                         "max_depth_a": c["max_depth_a"], "avg_hydropathy": c["avg_hydropathy"],
                         "mean_plddt": c["mean_plddt"], "pocket_score": c["pocket_score"],
                         **{f"comp_{k}": v for k, v in c["components"].items()},
                         "residues": ",".join(map(str, c["residues"]))})
        print(gene, acc, entry["best_pocket_score"], flush=True)
    pd.DataFrame(rows).to_parquet(PROCESSED / "structure_pockets.parquet", index=False)
    summary = {"modeled": sum(1 for r in results if r["modeled"]), "targets": len(genes) + len(controls),
               "config": {"min_plddt_keep": MIN_PLDDT, "pocket_volume": list(POCKET_VOLUME), "weights": WEIGHTS,
                          "top_pockets": TOP_POCKETS, "controls": CONTROLS},
               "tool": "pyKVFinder (SES surface, probe_in 1.4 A, probe_out 4.0 A, step 0.6 A)", "results": results}
    (PROCESSED / "structure_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return {k: v for k, v in summary.items() if k != "results"}


def load_graph(store=None) -> int:
    """One protein_structure record per modelled protein (AlphaFold model, best pockets, ligand-site check)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    s = json.loads((PROCESSED / "structure_summary.json").read_text(encoding="utf-8"))
    n = 0
    with store.batch():
        store.add_node("DATASET:ALPHAFOLD", "Dataset", "AlphaFold DB models with pyKVFinder cavities (PDBe-KB ligand-site check)")
        for r in s["results"]:
            node = store.resolve(r["gene"])
            if not r["modeled"] or not node or node["type"] != "Gene":
                continue
            v = r["ligand_site_validation"]
            rec = {"gene": r["gene"], "accession": r["uniprot"], "model": r["model_entity"], "version": r["version"],
                   "mean_plddt": r["mean_plddt"], "fraction_confident_residues": r["fraction_confident_residues"],
                   "length": r["length"], "n_cavities": r["n_cavities"], "best_pocket_score": r["best_pocket_score"],
                   "pockets": [{k: c[k] for k in ("rank", "volume_a3", "max_depth_a", "mean_plddt", "pocket_score")}
                               for c in r["pockets"]],
                   "ligand_site_recall": v.get("recall"), "ligand_site_recall_top3": v.get("recall_in_top_pockets"),
                   "known_ligands": v.get("n_ligands", 0), "url": r["url"],
                   "interpretation": "pocket shape of a predicted monomer; not a binding prediction"}
            ev = store.add_evidence(evidence_id("protein_structure", r["gene"]), "protein_structure", "DATASET:ALPHAFOLD",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "has_structure_in", "DATASET:ALPHAFOLD", ev)
            n += 1
    return n


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
