"""Record exactly what produced a set of outputs, so a run can be repeated or compared.

Every phase calls write(ctx, phase, modes) at start-up. It writes
<output_root>/metadata/provenance/<phase>_<mode>.json (and the merged config next to it) with:

  code_fingerprint   a short hash of every file that changes what the pipeline computes
                     (src/, scripts/, reference/, bin/). Works without git, e.g. for a ZIP
                     download from GitHub: two copies with the same fingerprint run the same code.
  git_commit         the commit, if the repo is a git checkout (plus "dirty" if edited)
  config_sha256      hash of the fully merged config (after `inherit:`), and a copy of it
  sample_sheet_sha256
  python, platform, and the versions of the packages in requirements.txt

10_report reads these files and shows them in a "Reproducibility" section, and warns if the
phases ran with different code or config (for example after editing the config between
phases without rerunning the earlier ones).
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import time
from importlib import metadata

import yaml

from . import __version__
from .config import REPO_ROOT

CODE_DIRS = ("src", "scripts", "reference", "bin")
CODE_SUFFIXES = {".py", ".yaml", ".yml", ".csv", ".sh", ".sbatch"}


def code_fingerprint() -> str:
    """sha256 over the repo-relative path and bytes of every code/reference file, 12 hex chars."""
    h = hashlib.sha256()
    for d in CODE_DIRS:
        base = REPO_ROOT / d
        if not base.exists():
            continue
        for p in sorted(base.rglob("*")):
            if not p.is_file() or p.suffix not in CODE_SUFFIXES or "__pycache__" in p.parts:
                continue
            h.update(str(p.relative_to(REPO_ROOT)).encode())
            h.update(b"\0")
            # normalize line endings so a Windows checkout gives the same fingerprint
            h.update(p.read_bytes().replace(b"\r\n", b"\n"))
            h.update(b"\0")
    return h.hexdigest()[:12]


def git_commit() -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        if r.returncode != 0:
            return None
        commit = r.stdout.strip()
        d = subprocess.run(["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--", *CODE_DIRS],
                           capture_output=True, text=True, timeout=5)
        return commit + (" (dirty)" if d.stdout.strip() else "")
    except Exception:
        return None


def package_versions() -> dict:
    out = {}
    req = REPO_ROOT / "requirements.txt"
    names = []
    if req.exists():
        for line in req.read_text().splitlines():
            line = line.split("#")[0].strip()
            if line:
                names.append(line.split("==")[0].split(">=")[0].strip())
    for n in names:
        try:
            out[n] = metadata.version(n)
        except metadata.PackageNotFoundError:
            out[n] = None
    return out


def sha256_file(p) -> str | None:
    try:
        return hashlib.sha256(p.read_bytes()).hexdigest()
    except OSError:
        return None


def collect(ctx) -> dict:
    cfg_text = yaml.safe_dump(ctx.cfg, sort_keys=True)
    return {
        "pipeline_version": __version__,
        "code_fingerprint": code_fingerprint(),
        "git_commit": git_commit(),
        "config_file": str(ctx.config_path),
        "config_sha256": hashlib.sha256(cfg_text.encode()).hexdigest()[:12],
        "sample_sheet_sha256": (sha256_file(ctx.sample_sheet) or "")[:12] or None,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": package_versions(),
    }


def write(ctx, phase: str, modes: list[str]) -> dict:
    rec = collect(ctx)
    rec.update({"phase": phase, "modes": list(modes), "started": time.strftime("%Y-%m-%d %H:%M:%S")})
    d = ctx.output_root / "metadata" / "provenance"
    d.mkdir(parents=True, exist_ok=True)
    # one record per (phase, mode), so per-mode HPC jobs don't overwrite each other and a
    # rerun of one mode replaces only that mode's record
    for m in (modes or [None]):
        name = phase + (f"_{m}" if m else "")
        (d / f"{name}.json").write_text(json.dumps(dict(rec, mode=m), indent=2))
        (d / f"{name}_config.yaml").write_text(yaml.safe_dump(ctx.cfg, sort_keys=False))
    return rec


def read_all(ctx) -> list[dict]:
    d = ctx.output_root / "metadata" / "provenance"
    if not d.exists():
        return []
    return [json.loads(p.read_text()) for p in sorted(d.glob("*.json"))]
