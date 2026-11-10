#!/usr/bin/env bash
# Run the orchestrated multi-agent analysis headless with Claude Code (subscription login, no API key).
# Linux / macOS equivalent of run_informe.ps1. Optional: ATLAS_RUN_DIR=replicates/rep1 to write a separate run.
set -euo pipefail
cd "$(dirname "$0")/.."
command -v uv >/dev/null || { echo "uv not found on PATH" >&2; exit 1; }
command -v claude >/dev/null || { echo "Claude Code CLI not found on PATH" >&2; exit 1; }
allowed="mcp__expression-mcp,mcp__kg-mcp,mcp__external-mcp,mcp__cell-mcp,mcp__pathway-mcp,Bash(uv run python scripts/*),Read,Write,Edit,Agent,Task"
exec claude -p "/informe" --model "${ATLAS_ORCHESTRATOR_MODEL:-claude-opus-5-5}" --mcp-config .mcp.json \
  --allowedTools "$allowed" "$@"
