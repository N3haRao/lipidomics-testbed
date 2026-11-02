"""Phase 00: check the config and sample sheet, and cache the reference data (once).

    python scripts/00_setup.py --config config.yaml                 # just validate
    python scripts/00_setup.py --config config.yaml --download      # + LIPID MAPS, HMDB, Allen atlas

If a download is blocked (some sites want you to click through a page), download the file
in your browser and point to it instead:

    python scripts/00_setup.py --config config.yaml --lipidmaps-file ~/Downloads/LMSD.sdf.zip \
                                                    --hmdb-file ~/Downloads/structures.zip

Everything is cached under <cache_root> (default <output_root>/cache). Later phases only
read the cache, they never go online. Nothing is ever uploaded anywhere.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, refdb, samples  # noqa: E402
from cln3msi.report import Summary  # noqa: E402


def extra(p):
    p.add_argument("--download", action="store_true", help="download LIPID MAPS, HMDB and the Allen atlas")
    p.add_argument("--lipidmaps-file", default=None, help="already downloaded LMSD .sdf or .zip")
    p.add_argument("--hmdb-file", default=None, help="already downloaded HMDB structures .sdf or .zip")
    p.add_argument("--skip-atlas", action="store_true", help="don't touch the atlas")


def cache_db(ctx, log, which: str, local: str | None, url: str, download: bool, S: Summary):
    target = ctx.cache(which, f"{which}.parquet")
    if local:
        src = Path(local).expanduser()
    elif download:
        src = ctx.cache(which, "raw_download" + (".zip" if url.lower().split("?")[0].endswith(".zip")
                                                   or "zip" in url.lower() else ".sdf"))
        try:
            refdb.download(url, src, ctx.get("downloads.timeout_s", 600), log)
        except Exception as e:
            S.warn(f"{which}: download failed ({e}). Download it by hand and rerun with "
                   f"--{which}-file <path>.", log)
            return
    elif target.exists():
        S.p(f"- {which}: cached at `{target}`")
        return
    else:
        S.p(f"- {which}: not cached yet (run with --download or --{which}-file)")
        return
    df = refdb.sdf_to_table(src, log)
    df.to_parquet(target)
    S.p(f"- {which}: {len(df)} entries with formulas cached at `{target}` (from `{src.name}`)")


def main():
    t0 = time.time()
    ctx, log, _, args = cli.start("00_setup", __doc__, mode=False, extra_args=extra, check_data=False)
    S = Summary(ctx, "00_setup", None, "Setup and validation")

    S.h("Config")
    S.kv({"config": ctx.config_path, "data_root": ctx.data_root, "output_root": ctx.output_root,
          "cache_root": ctx.cache_root, "modes to run": ", ".join(ctx.modes_to_run)})

    S.h("Data folder")
    from cln3msi.config import ConfigError, check_data_layout
    try:
        check_data_layout(ctx)
        S.p("Data folder and mode subfolders found.")
    except ConfigError as e:
        S.warn(f"Data folder check failed (fine if you're setting up before copying the data):\n{e}", log)

    S.h("Sample sheet")
    try:
        sheet = samples.load_sheet(ctx)
        S.p(f"{len(sheet)} animals on {sheet.slide.nunique()} slides; every slide has one Control (Top) "
            "and one Mutant (Bottom) according to the sheet.")
        S.table(samples.section_table(sheet))
        log.info("sample sheet OK")
    except samples.SampleSheetError as e:
        S.warn(str(e), log)
        S.write()
        sys.exit(1)

    S.h("Reference databases")
    cache_db(ctx, log, "lipidmaps", args.lipidmaps_file, ctx.get("downloads.lipidmaps_url"), args.download, S)
    cache_db(ctx, log, "hmdb", args.hmdb_file, ctx.get("downloads.hmdb_url"), args.download, S)

    S.h("Atlas")
    src = ctx.get("registration.atlas_source", "brainglobe")
    if args.skip_atlas:
        S.p("skipped (--skip-atlas)")
    elif src != "brainglobe":
        p = Path(src) if Path(src).is_absolute() else (ctx.config_path.parent / src)
        S.p(f"atlas_source is a file: `{p}` ({'found' if p.exists() else 'NOT FOUND'})")
    elif args.download:
        try:
            from brainglobe_atlasapi import BrainGlobeAtlas
            d = ctx.cache("atlas", "x").parent
            a = BrainGlobeAtlas(ctx.get("registration.atlas_name"), brainglobe_dir=d, check_latest=False)
            S.p(f"cached {getattr(a, 'atlas_name', ctx.get('registration.atlas_name'))}: annotation {a.annotation.shape}, resolution {a.resolution} um, in `{d}`")
        except Exception as e:
            S.warn(f"atlas download failed: {e}. 07_register will fall back to lipid WM/GM masks.", log)
    elif (ctx.cache_root / "atlas").is_dir() and any((ctx.cache_root / "atlas").iterdir()):
        S.p(f"cached in `{ctx.cache_root / 'atlas'}`")
    else:
        S.p("not downloaded yet (run with --download). 07_register falls back to lipid WM/GM masks without it.")

    S.h("Optional packages")
    for mod, why in [("uMAIA", "optional cross-slide normalization in 06"),
                     ("harmonypy", "optional batch-corrected embedding for plots in 06/08"),
                     ("leidenalg", "Leiden territories in 08 (fallback: KMeans)"),
                     ("SimpleITK", "registration in 07 (needed)"),
                     ("brainglobe_atlasapi", "Allen atlas in 07")]:
        try:
            __import__(mod)
            S.p(f"- {mod}: available ({why})")
        except Exception:
            S.p(f"- {mod}: not installed ({why})")

    log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
