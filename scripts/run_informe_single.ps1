# Single-agent arm of the agent benchmark (M17): one Opus agent with all MCP tools, no subagents.
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)
$env:Path = [Environment]::GetEnvironmentVariable("Path", "User") + ";" + [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" + $env:Path
$claude = (Get-Command claude -ErrorAction SilentlyContinue).Source
if (-not $claude) {
  $claude = Get-ChildItem "$env:LOCALAPPDATA\Claude-Data\*\claude-code\*\*\claude.exe" -ErrorAction SilentlyContinue |
            Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
}
$allowed = "mcp__expression-mcp,mcp__kg-mcp,mcp__external-mcp,mcp__cell-mcp,mcp__pathway-mcp,Bash(uv run python scripts/*),PowerShell(uv run python scripts/*),Read,Write,Edit"
& $claude -p "/informe_single" --model claude-opus-5-5 --mcp-config .mcp.json --allowedTools $allowed --output-format json | Out-File -Encoding utf8 (Join-Path (Split-Path -Parent $PSScriptRoot) "runs/last_informe_session.json")  # UTF-8, not the UTF-16 of ">"
"exit=$LASTEXITCODE"
