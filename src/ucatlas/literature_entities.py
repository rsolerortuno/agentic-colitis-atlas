"""Deterministic entity resolution for extracted literature relations (strings -> graph node ids)."""
from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from .paths import LITERATURE, PROCESSED, RAW


@lru_cache(maxsize=1)
def symbol_lookup() -> dict[str, str]:
    """Approved symbols, dataset-era symbols and unique previous/alias symbols -> approved HGNC symbol."""
    hg = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "prev_symbol", "alias_symbol"],
                     dtype=str)
    lookup: dict[str, str] = {}
    ambiguous: set[str] = set()
    for col in ("alias_symbol", "prev_symbol"):
        for sym, al in zip(hg["symbol"], hg[col]):
            if isinstance(al, str):
                for a in al.split("|"):
                    a = a.strip().upper()
                    if a in lookup and lookup[a] != sym:
                        ambiguous.add(a)
                    lookup[a] = sym
    for a in ambiguous:
        lookup.pop(a, None)
    smp = PROCESSED / "symbol_map.parquet"  # dataset-era symbols seen in the cohorts (written by the expression step)
    if smp.exists():
        sm = pd.read_parquet(smp).dropna(subset=["symbol"])
        lookup.update({d.upper(): s for d, s in zip(sm["dataset_symbol"], sm["symbol"])})
    lookup.update({s.upper(): s for s in hg["symbol"]})
    return lookup


def resolve_gene(text: str) -> str | None:
    t = text.strip()
    t = re.sub(r"\s+(gene|protein|expression|mRNA|levels?)$", "", t, flags=re.I)
    if not re.fullmatch(r"[A-Za-z0-9\-]{2,15}", t):
        return None
    return symbol_lookup().get(t.upper())


def is_disease(text: str) -> bool:
    """Does the text name the disease of this atlas (config.yaml literature.disease_pattern)?"""
    from .config import load
    return bool(re.search(load()["literature"]["disease_pattern"], text, re.I))


def gene_literature_counts() -> pd.Series:
    """Number of validated relations in which each gene is subject or object."""
    path = LITERATURE / "relations.parquet"
    if not path.exists():
        return pd.Series(dtype=float)
    rel = pd.read_parquet(path)
    genes = [g for col in ("subject", "object") for g in rel[col].map(resolve_gene) if g]
    return pd.Series(genes).value_counts()
