"""Immutable input inventory: verify bytes, never refresh expectations during a run."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
from fnmatch import fnmatchcase
import json
from pathlib import Path

from .paths import DATA, ROOT

MANIFEST = ROOT / "reproducibility" / "inputs.json"


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify(manifest: Path = MANIFEST, data: Path = DATA, project: Path = ROOT) -> dict:
    expected = json.loads(manifest.read_text(encoding="utf-8"))
    errors = []
    for key, name in (("config_sha256", "config.yaml"), ("uv_lock_sha256", "uv.lock")):
        if key in expected:
            path = project / name
            if not path.is_file() or sha256(path) != expected[key]:
                errors.append({"path":name, "error":"project configuration/lock checksum mismatch"})
    root = data.resolve()
    for entry in expected["files"]:
        path = (root / entry["path"]).resolve()
        if not path.is_relative_to(root):
            errors.append({"path": entry["path"], "error": "outside data root"})
        elif not path.is_file():
            errors.append({"path": entry["path"], "error": "missing"})
        elif path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            errors.append({"path": entry["path"], "error": "checksum mismatch"})
    # A new cached response can change the analysis just as a changed file can.
    # Declared input trees must match the snapshot inventory, not merely contain it.
    listed = {entry["path"] for entry in expected["files"]}
    for directory in expected.get("input_trees", []):
        tree = (root / directory).resolve()
        if not tree.is_relative_to(root):
            errors.append({"path": directory, "error": "outside data root"})
            continue
        for path in sorted(tree.rglob("*")):
            relative = path.relative_to(root).as_posix()
            generated = any(fnmatchcase(relative, pattern)
                            for pattern in expected.get("generated_exclusions", []))
            if path.is_file() and not generated and relative not in listed:
                errors.append({"path": path.relative_to(root).as_posix(),
                               "error": "unlisted input in frozen tree"})
    project_root = project.resolve()
    for entry in expected.get("project_files", []):
        path = (project_root / entry["path"]).resolve()
        if not path.is_relative_to(project_root):
            errors.append({"path": entry["path"], "error": "outside project root"})
        elif not path.is_file() or sha256(path) != entry["sha256"]:
            errors.append({"path": entry["path"], "error": "saved analysis input checksum mismatch"})
    project_listed = {entry["path"] for entry in expected.get("project_files", [])}
    for pattern in expected.get("project_input_globs", []):
        if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
            errors.append({"path": pattern, "error": "outside project root"})
            continue
        for path in sorted(project_root.glob(pattern)):
            relative = path.relative_to(project_root).as_posix()
            if path.is_file() and relative not in project_listed:
                errors.append({"path": relative, "error": "unlisted saved analysis input"})
    analysis_run = expected.get("analysis_run")
    if analysis_run:
        selected = (project_root / analysis_run).resolve()
        if not selected.is_relative_to(project_root / "runs"):
            errors.append({"path": analysis_run, "error": "analysis run outside runs directory"})
        elif (selected / "narrative.json").relative_to(project_root).as_posix() not in project_listed:
            errors.append({"path": analysis_run, "error": "selected narrative is not frozen"})
    return {"checked": len(expected["files"]),
            "project_checked": len(expected.get("project_files", [])),
            "analysis_run": analysis_run,
            "errors": errors, "ok": not errors}


def require_verified(manifest: Path = MANIFEST, data: Path = DATA, project: Path = ROOT) -> dict:
    result = verify(manifest, data, project)
    if not result["ok"]:
        raise RuntimeError("Frozen input verification failed; expectations were not changed: " + json.dumps(result["errors"]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--data", type=Path, default=DATA)
    args = parser.parse_args()
    result = verify(args.manifest, args.data)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["ok"] else 1)


if __name__ == "__main__":
    main()


@contextmanager
def offline_inputs():
    """Fail on a cache miss instead of silently fetching a newer remote response.

    All HTTP clients in this pipeline use requests. Patch its common transport
    entry point so retries cannot bypass the frozen-input policy. Restore it on
    exit, including on failure; this does not change system network settings.
    """
    from unittest.mock import patch

    def reject_request(self, method, url, *args, **kwargs):
        raise RuntimeError(
            f"Frozen reproduction attempted network access ({method} {url}). "
            "Restore the missing input from the published snapshot; "
            "do not refresh the expected hashes during reproduction."
        )

    with patch("requests.sessions.Session.request", reject_request):
        yield


def freeze(run: str = "2026-11-09_full") -> dict:
    """Write reproducibility/inputs.json from the download manifests (no re-hashing of raw files), the config and lock
    hashes and the frozen agent outputs of the published run. `run-all` checks it before rebuilding."""
    raw = DATA / "raw"
    files = []
    top = json.loads((raw / "manifest.json").read_text(encoding="utf-8"))
    files += [{"path": f"raw/{k}", "sha256": v["sha256"], "bytes": v["bytes"]} for k, v in sorted(top.items())]
    for sub in ("scp259", "scp3755"):
        m = json.loads((raw / sub / "manifest.json").read_text(encoding="utf-8"))
        files += [{"path": f"raw/{sub}/{k}", "sha256": v["sha256"], "bytes": v["bytes"]} for k, v in sorted(m["files"].items())]
    for gse in ("GSE312420", "GSE339323", "GSE189184", "GSE312415"):
        sp = raw / "spatial" / gse / "manifest.json"
        if sp.exists():
            m = json.loads(sp.read_text(encoding="utf-8"))
            files += [{"path": f"raw/spatial/{gse}/{k}", "sha256": v["sha256"], "bytes": v["bytes"]} for k, v in sorted(m.items())]
    run_dir = ROOT / "runs" / run
    project = [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha256(p)}
               for p in sorted(run_dir.glob("agent_*.json")) + [run_dir / "narrative.json"]
               if not p.name.endswith(".verified.json")]
    out = {"config_sha256": sha256(ROOT / "config.yaml"), "uv_lock_sha256": sha256(ROOT / "uv.lock"),
           "files": files, "project_files": project, "analysis_run": f"runs/{run}"}
    MANIFEST.parent.mkdir(exist_ok=True)
    MANIFEST.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return {"files": len(files), "project_files": len(project)}
