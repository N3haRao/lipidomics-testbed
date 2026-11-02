"""Load config.yaml, check it, and hand out output paths.

Every script starts with `ctx = load(args.config)`. The Context object knows where the data
is, where outputs go, and refuses to give out any write path inside the raw data folder.
That refusal is the safety net for the "raw data is read-only" rule: all writes in the
pipeline go through ctx.out(...) or ctx.cache(...), so nothing can land in data_root by
accident.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]   # src/cln3msi/config.py -> repo root

# Keys that must be present. Anything else has a default in the code or is optional.
REQUIRED = {
    "paths": ["data_root", "imzml_subdir", "output_root", "sample_sheet"],
    "modes": ["run", "definitions", "file_regex"],
    "qc": ["slide", "mode"],
}

OUTPUT_KINDS = ("metadata", "processed", "results", "figures", "logs", "summaries", "report")


class ConfigError(RuntimeError):
    """Raised for anything wrong in the config, with a message meant for a human."""


def _resolve(p: str | Path, base: Path) -> Path:
    """Expand ~ and $VARIABLES, make relative paths relative to `base` (the repo or the config file).

    $VARIABLES let one config work on a cluster where the work folder is only known at run
    time, e.g. data_root: "${CLN3_WORK}/data" (config.hpc.yaml does that).
    """
    raw = str(p)
    s = os.path.expandvars(raw)
    if "$" in s:
        raise ConfigError(f"path '{raw}' uses an environment variable that isn't set "
                          f"(e.g. export CLN3_WORK=/path/to/work before running)")
    p = Path(s).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _is_inside(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


@dataclass
class Context:
    cfg: dict
    config_path: Path
    data_root: Path
    imzml_root: Path
    output_root: Path
    cache_root: Path
    sample_sheet: Path

    # ---- convenient getters ---------------------------------------------------------
    def get(self, dotted: str, default: Any = None) -> Any:
        """ctx.get("peaks.snr") instead of ctx.cfg["peaks"]["snr"], with a default."""
        node: Any = self.cfg
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    @property
    def modes_to_run(self) -> list[str]:
        return list(self.cfg["modes"]["run"])

    def mode_def(self, mode: str) -> dict:
        defs = self.cfg["modes"]["definitions"]
        if mode not in defs:
            raise ConfigError(f"mode '{mode}' is not in modes.definitions (have: {list(defs)})")
        return defs[mode]

    def mode_dir(self, mode: str) -> Path:
        return self.imzml_root / self.mode_def(mode)["folder"]

    # ---- output paths ---------------------------------------------------------------
    def out(self, kind: str, *parts: str, mkdir: bool = True) -> Path:
        """Path under <output_root>/<kind>/<parts...>. Parent folders get created.

        Example: ctx.out("figures", "02_qc", "lip_neg", "tic_E.png")
        """
        if kind not in OUTPUT_KINDS:
            raise ValueError(f"unknown output kind {kind}, use one of {OUTPUT_KINDS}")
        p = self.output_root.joinpath(kind, *parts)
        self._check_writable(p)
        if mkdir:
            p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def cache(self, *parts: str, mkdir: bool = True) -> Path:
        p = self.cache_root.joinpath(*parts)
        self._check_writable(p)
        if mkdir:
            p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def _check_writable(self, p: Path) -> None:
        if _is_inside(p, self.data_root):
            raise ConfigError(
                f"refusing to write {p}: it is inside data_root ({self.data_root}). "
                "Raw data is read-only; set paths.output_root to somewhere else."
            )

    def phase_json(self, phase: str, mode: str, name: str) -> dict | None:
        """Read a small JSON result an earlier phase left behind, or None if it isn't there."""
        p = self.output_root / "results" / phase / mode / name
        if p.exists():
            return json.loads(p.read_text())
        return None


def _deep_merge(base: dict, over: dict) -> dict:
    out = dict(base)
    for k, v in over.items():
        out[k] = _deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def read_yaml(config_path: Path) -> dict:
    """Read a config. A config can say `inherit: other.yaml` and only list what it changes."""
    cfg = yaml.safe_load(config_path.read_text())
    if not isinstance(cfg, dict):
        raise ConfigError(f"{config_path} did not parse into a mapping")
    parent = cfg.pop("inherit", None)
    if parent:
        cfg = _deep_merge(read_yaml(_resolve(parent, config_path.parent)), cfg)
    return cfg


def load(config_path: str | Path, check_data: bool = False) -> Context:
    """Read and validate the YAML config. check_data=True also checks the data folder exists."""
    config_path = Path(config_path).expanduser().resolve()
    if not config_path.exists():
        raise ConfigError(f"config file not found: {config_path}")
    cfg = read_yaml(config_path)

    for section, keys in REQUIRED.items():
        if section not in cfg:
            raise ConfigError(f"config is missing the '{section}:' block")
        for k in keys:
            if k not in cfg[section]:
                raise ConfigError(f"config is missing '{section}.{k}'")

    paths = cfg["paths"]
    # Relative data/output paths are taken relative to the config file (handy for the test
    # config), the sample sheet relative to the repo (that's where metadata/ lives).
    base = config_path.parent
    data_root = _resolve(paths["data_root"], base)
    output_root = _resolve(paths["output_root"], base)
    cache_root = _resolve(paths["cache_root"], base) if paths.get("cache_root") else output_root / "cache"
    sheet = _resolve(paths["sample_sheet"], REPO_ROOT)
    if not sheet.exists():
        sheet = _resolve(paths["sample_sheet"], base)

    if _is_inside(output_root, data_root):
        raise ConfigError(
            f"output_root ({output_root}) is inside data_root ({data_root}). "
            "Raw data folders must stay untouched, so pick an output folder outside it."
        )

    defs = cfg["modes"]["definitions"]
    for m in cfg["modes"]["run"]:
        if m not in defs:
            raise ConfigError(f"modes.run lists '{m}' but modes.definitions has no entry for it")
        for k in ("folder", "polarity", "analyte"):
            if k not in defs[m]:
                raise ConfigError(f"modes.definitions.{m} is missing '{k}'")
        if defs[m]["polarity"] not in ("negative", "positive"):
            raise ConfigError(f"modes.definitions.{m}.polarity must be negative or positive")

    ctx = Context(cfg=cfg, config_path=config_path, data_root=data_root,
                  imzml_root=data_root / paths["imzml_subdir"], output_root=output_root,
                  cache_root=cache_root, sample_sheet=sheet)
    if check_data:
        check_data_layout(ctx)
    return ctx


EXPECTED_LAYOUT = """
Expected layout under data_root:

  <data_root>/
    {imzml_subdir}/
{mode_lines}
  e.g. slide_e_lip_neg.imzML + slide_e_lip_neg.ibd (same base name, same folder)
"""


def check_data_layout(ctx: Context) -> None:
    """Fail early, with the expected layout printed, if the data folder isn't where the config says."""
    lines = "\n".join(
        f"      {ctx.mode_def(m)['folder']}/   slide_<letter>_{m}.imzML + .ibd" for m in ctx.modes_to_run
    )
    layout = EXPECTED_LAYOUT.format(imzml_subdir=ctx.cfg["paths"]["imzml_subdir"], mode_lines=lines)
    if not ctx.data_root.is_dir():
        raise ConfigError(f"data_root does not exist or is not a folder:\n  {ctx.data_root}\n{layout}\n"
                          "Fix paths.data_root in your config.")
    if not ctx.imzml_root.is_dir():
        raise ConfigError(f"imzML folder not found:\n  {ctx.imzml_root}\n{layout}")
    missing = [m for m in ctx.modes_to_run if not ctx.mode_dir(m).is_dir()]
    if missing:
        raise ConfigError(f"mode folder(s) missing for {missing}:\n" +
                          "\n".join(f"  {ctx.mode_dir(m)}" for m in missing) + f"\n{layout}")
