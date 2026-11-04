"""Shared command-line and logging setup, so every phase script starts the same way.

    ctx, log, modes = cli.start("04_peaks_matrix", __doc__)

gives you the loaded config, a logger that writes to the screen AND to
<output_root>/logs/<phase>_<mode>_<timestamp>.log, and the list of modes to loop over
(from --mode, or modes.run in the config when --mode is not given).
"""
from __future__ import annotations

import argparse
import logging
import platform
import sys
import time

import psutil

from . import config as config_mod
from . import provenance


def memory_gb() -> float:
    """Resident memory of this process in GB (what Activity Monitor calls 'Memory')."""
    return psutil.Process().memory_info().rss / 1e9


class MemLogger:
    """Tiny helper: log.mem("after reading slide E") writes the current RSS to the log."""

    def __init__(self, log: logging.Logger):
        self.log = log
        self.peak = 0.0

    def __call__(self, where: str) -> None:
        gb = memory_gb()
        self.peak = max(self.peak, gb)
        self.log.info(f"[memory] {where}: {gb:.2f} GB (peak so far {self.peak:.2f} GB)")


def parser(description: str, slide: bool = False, mode: bool = True) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="path to config.yaml")
    if mode:
        p.add_argument("--mode", action="append", default=None,
                       help="ion mode to run (lip_neg, lip_pos, met_neg). Repeat for several. "
                            "Default: modes.run in the config")
    if slide:
        p.add_argument("--slide", default=None, help="slide letter (E, F, ...). Default: qc.slide in the config")
    return p


def start(phase: str, description: str, slide: bool = False, mode: bool = True,
          extra_args=None, argv=None, check_data: bool = True):
    """Parse args, load config, set up logging. Returns (ctx, log, modes, args)."""
    p = parser(description, slide=slide, mode=mode)
    if extra_args:
        extra_args(p)
    args = p.parse_args(argv)
    try:
        ctx = config_mod.load(args.config, check_data=check_data)
    except config_mod.ConfigError as e:
        print(f"\nCONFIG / DATA ERROR\n{e}\n", file=sys.stderr)
        sys.exit(2)

    modes = (args.mode if getattr(args, "mode", None) else ctx.modes_to_run) if mode else []
    for m in modes:
        ctx.mode_def(m)   # raises a clear error for a typo like lipneg

    tag = "_".join(modes) if modes and len(modes) < 3 else ("all" if modes else "")
    stamp = time.strftime("%Y%m%d-%H%M%S")
    logfile = ctx.out("logs", f"{phase}{'_' + tag if tag else ''}_{stamp}.log")

    log = logging.getLogger(phase)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(logfile)):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.propagate = False
    log.mem = MemLogger(log)   # type: ignore[attr-defined]

    log.info(f"phase {phase} | config {ctx.config_path}")
    log.info(f"python {platform.python_version()} on {platform.platform()} | "
             f"{psutil.cpu_count()} CPUs | {psutil.virtual_memory().total / 1e9:.1f} GB RAM")
    log.info(f"data_root   {ctx.data_root}")
    log.info(f"output_root {ctx.output_root}")
    if modes:
        log.info(f"modes: {modes}")
    log.info(f"log file: {logfile}")
    try:
        rec = provenance.write(ctx, phase, modes)
        log.info(f"code fingerprint {rec['code_fingerprint']}"
                 + (f" | git {rec['git_commit']}" if rec["git_commit"] else "")
                 + f" | config {rec['config_sha256']} | pipeline v{rec['pipeline_version']}")
    except Exception as e:      # provenance must never stop a run
        log.warning(f"could not write provenance record: {e}")
    log.mem("start")
    return ctx, log, modes, args


def finish(log: logging.Logger, t0: float) -> None:
    log.info(f"done in {(time.time() - t0) / 60:.1f} min, peak memory {log.mem.peak:.2f} GB")   # type: ignore
