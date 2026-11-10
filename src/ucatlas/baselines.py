"""Review (2026-11-10): what does the agentic layer add? The kidney atlas compared its pipeline with a closed-book LLM.

Closed-book: Claude Opus 5.5 with no tools, no MCP servers and no project files, asked the same question three times
for its 10 strongest candidates (JSON list). Pipeline: the target lists of the three team runs. Measures: mean
pairwise Jaccard of the lists, share of candidates that are UC clinical targets (the M8 gold set), share that are in
the data-only top 500 (target or biomarker profile), and statements citing checkable evidence.
Output: evals/baselines.json.
"""
from __future__ import annotations

import itertools
import json
import os
import subprocess
import tempfile

from .paths import ROOT

PROMPT = ("Question: what are the candidate therapeutic targets and biomarkers in ulcerative colitis with the strongest "
          "convergent evidence? Answer from your own knowledge only. Return ONLY a JSON list of the 10 strongest "
          "candidate genes as HGNC symbols, strongest first, for example [\"GENE1\", \"GENE2\"].")
TEAM_RUNS = ["2026-11-09_full", "2026-11-09_rep1", "2026-11-09_rep2"]


def closed_book(n: int = 3, model: str = "claude-opus-5-5") -> list[list[str]]:
    from .qa import find_claude
    claude = find_claude()
    out = []
    for _ in range(n):
        with tempfile.TemporaryDirectory() as tmp:
            mcp = os.path.join(tmp, "mcp.json")
            open(mcp, "w").write('{"mcpServers":{}}')
            p = subprocess.run([claude, "-p", PROMPT, "--model", model, "--tools", "", "--strict-mcp-config",
                                "--mcp-config", mcp, "--setting-sources", "", "--no-session-persistence"],
                               capture_output=True, text=True, encoding="utf-8", cwd=tmp, timeout=600)
        txt = p.stdout
        start, end = txt.find("["), txt.rfind("]")
        out.append([g.strip().upper() for g in json.loads(txt[start:end + 1])][:10])
    return out


def _jaccard(lists):
    pairs = list(itertools.combinations(lists, 2))
    return round(sum(len(set(a) & set(b)) / len(set(a) | set(b)) for a, b in pairs) / len(pairs), 3)


def run(lists_cb: list[list[str]] | None = None) -> dict:
    from .evaluate import clinical_targets
    clin = set(clinical_targets()["gene"])
    conv = json.loads((ROOT / "runs" / "2026-11-09_full" / "convergence.json").read_text(encoding="utf-8"))
    import pandas as pd
    from .integrate import compute
    from .profiles import components
    from .config import load
    t = components(compute())
    data_top = set()
    for name, w in load()["profiles"].items():
        s = sum(wt * t[f"c_{k}"] for k, wt in w.items())
        data_top |= set(t.assign(s=s).nlargest(500, "s")["gene"])
    lists_cb = lists_cb or closed_book()
    team = [[x["gene"] for x in json.loads((ROOT / "runs" / r / "narrative.json").read_text(encoding="utf-8"))["targets"]]
            for r in TEAM_RUNS]
    statements = [json.loads((ROOT / "runs" / r / "verification_final.json").read_text(encoding="utf-8")) for r in TEAM_RUNS]

    def describe(lists):
        return {"lists": lists, "mean_pairwise_jaccard": _jaccard(lists),
                "share_clinical_targets": [round(len(set(l) & clin) / len(l), 2) for l in lists],
                "share_in_data_top500": [round(len(set(l) & data_top) / len(l), 2) for l in lists]}

    out = {"closed_book": {**describe(lists_cb), "statements_citing_evidence": 0},
           "pipeline": {**describe(team), "statements_checked": [s["statements"] for s in statements],
                        "unsupported": [s["unsupported"] for s in statements]},
           "n_clinical_targets": len(clin)}
    (ROOT / "evals" / "baselines.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    r = run()
    for k in ("closed_book", "pipeline"):
        v = r[k]
        print(k, v["mean_pairwise_jaccard"], v["share_clinical_targets"], v["share_in_data_top500"])
        print("  ", v["lists"])
