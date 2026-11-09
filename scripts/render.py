"""Render reports/uc_report.html from the latest runs/<date>/findings.json (or a given date).

With ATLAS_RUN_DIR set (replicate runs, see M17) the report is written next to that run's
findings.json instead, so a replicate never overwrites the published report."""
import os
import sys

from ucatlas.runs import current_run, run_dir
from ucatlas.render import render

if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(render(run_dir(sys.argv[1]) / "findings.json"))
    elif os.environ.get("ATLAS_RUN_DIR"):
        d = run_dir()
        print(render(d / "findings.json", out=d / "report.html"))
    else:
        print(render(current_run() / "findings.json"))
