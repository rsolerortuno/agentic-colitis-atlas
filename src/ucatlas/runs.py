"""Run folders of the agent step: runs/<date>/ or ATLAS_RUN_DIR (for independent replicate runs)."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

from .paths import RUNS

CURRENT = RUNS / ".current"


def run_dir(day: str | None = None) -> Path:
    """runs/<date>/ by default; ATLAS_RUN_DIR (relative to runs/ or absolute) redirects a whole /informe run."""
    override = os.environ.get("ATLAS_RUN_DIR")
    if override and not day:
        d = Path(override) if Path(override).is_absolute() else RUNS / override
    else:
        d = RUNS / (day or date.today().isoformat())
    d.mkdir(parents=True, exist_ok=True)
    return d


def current_run() -> Path:
    """ATLAS_RUN_DIR, else the folder the integration step started (runs/.current), else the latest dated run.
    A step repeated on a later day must not silently switch to a new, empty folder."""
    if os.environ.get("ATLAS_RUN_DIR"):
        return run_dir()
    if CURRENT.exists():
        d = RUNS / CURRENT.read_text(encoding="utf-8").strip()
        if d.is_dir():
            return d
    runs = sorted(p for p in RUNS.glob("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]*") if p.is_dir())
    return runs[-1] if runs else run_dir()


def set_current(d: Path) -> None:
    CURRENT.write_text(d.name, encoding="utf-8")
