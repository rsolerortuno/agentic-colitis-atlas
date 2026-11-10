"""Ulcerative colitis agentic multi-omics GraphRAG.

Linear algebra runs on one thread: multithreaded BLAS sums in a different order from run to run, which changed the
last binary digits of least-squares results between two rebuilds (2026-11-10, evals/history.md). The variables only
take effect if this package is imported before numpy, which is the case for every entry point (ucatlas CLI, scripts,
MCP servers); an explicit setting in the environment wins.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
