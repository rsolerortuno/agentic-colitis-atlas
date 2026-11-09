"""M8 (c): known-answer questions (evals/qa.yaml) answered by Claude Code headless (`claude -p`) with the project
MCP servers, graded by regex. Results in evals/qa_answers.json."""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import yaml

from .paths import ROOT

EVALS = ROOT / "evals"


def find_claude() -> str | None:
    exe = shutil.which("claude")
    if exe:
        return exe
    import os
    base = Path(os.environ.get("LOCALAPPDATA", "")) / "Claude-Data"
    found = sorted(base.glob("*/claude-code/*/*/claude.exe"), key=lambda p: p.stat().st_mtime, reverse=True)
    return str(found[0]) if found else None


QA_PROMPT = ("You are answering a factual question about this project's ulcerative colitis analysis. Use ONLY the MCP tools "
             "(expression-mcp, kg-mcp, external-mcp) to look it up; do not answer from memory. "
             "Reply with one short sentence.\nQuestion: {q}")


def run_qa(model: str = "claude-sonnet-5-5", timeout: int = 300) -> dict:
    qa = yaml.safe_load((EVALS / "qa.yaml").read_text(encoding="utf-8"))
    claude = find_claude()
    allowed = "mcp__expression-mcp,mcp__kg-mcp,mcp__external-mcp"
    results = []
    for item in qa:
        answer, error = "", None
        if not claude:
            error = "claude CLI not found"
        else:
            try:
                p = subprocess.run([claude, "-p", QA_PROMPT.format(q=item["question"]), "--model", model,
                                    "--mcp-config", str(ROOT / ".mcp.json"), "--allowedTools", allowed],
                                   capture_output=True, text=True, timeout=timeout, cwd=ROOT, encoding="utf-8")
                answer = (p.stdout or "").strip()
                if p.returncode != 0 or "Not logged in" in answer:
                    error = answer[:200] or p.stderr[:200]
            except subprocess.TimeoutExpired:
                error = "timeout"
        correct = bool(not error and re.search(item["pattern"], answer, re.I))
        results.append({"id": item["id"], "question": item["question"], "expected": item["expected"],
                        "answer": answer[:500], "correct": correct, "error": error})
    return {"model": model, "correct": sum(r["correct"] for r in results), "total": len(results), "results": results}


def grade_answers(path: Path, model: str) -> dict:
    """Grade a saved answers file ([{"id", "answer"}]) with the same regexes as run_qa."""
    qa = {q["id"]: q for q in yaml.safe_load((EVALS / "qa.yaml").read_text(encoding="utf-8"))}
    answers = {a["id"]: a["answer"] for a in json.loads(Path(path).read_text(encoding="utf-8"))}
    results = []
    for qid, item in qa.items():
        ans = answers.get(qid, "")
        results.append({"id": qid, "question": item["question"], "expected": item["expected"], "answer": ans[:500],
                        "correct": bool(re.search(item["pattern"], ans, re.I)), "error": None if ans else "no answer"})
    return {"model": model, "correct": sum(r["correct"] for r in results), "total": len(results), "results": results}


if __name__ == "__main__":
    r = run_qa()
    (EVALS / "qa_answers.json").write_text(json.dumps(r, indent=1), encoding="utf-8")
    print(r["correct"], "/", r["total"])
