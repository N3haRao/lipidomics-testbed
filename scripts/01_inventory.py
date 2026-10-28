"""Phase 01: what's in each imzML file? Read the headers, match the files to the sample sheet.

    python scripts/01_inventory.py --config config.yaml              # every mode in modes.run
    python scripts/01_inventory.py --config config.yaml --mode lip_neg

For each file: spectrum type, storage mode, polarity, instrument and software cvParams,
pixel grid, number of spectra, m/z range, data types. The XML is parsed for coordinates and
byte offsets (no intensities). To check what the header *claims* against the data, we read
the m/z arrays of ~200 spectra and the intensities of 30 (a few MB, not the slide).

Writes
  <output_root>/metadata/files.csv           one row per file (all modes run so far)
  <output_root>/metadata/header_<file>.txt   every header cvParam, as written by the exporter
  <output_root>/results/01_inventory/<mode>/inventory.json
  <output_root>/summaries/01_inventory_<mode>.md
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, samples  # noqa: E402
from cln3msi.imzml_io import SlideReader, detect_profile, read_header, summarize_header  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def inspect_file(ctx, path: Path, mode: str, log, rng) -> dict:
    row = {"mode": mode, "file": path.name, "path": str(path)}
    info = samples.parse_filename(path.name, ctx.cfg["modes"]["file_regex"])
    row["slide"] = info["slide"] if info else None
    row["name_matches_convention"] = info is not None and info["mode"] == mode
    ibd = path.with_suffix(".ibd")
    row["ibd_found"] = ibd.exists()
    row["imzml_MB"] = round(path.stat().st_size / 1e6, 1)
    row["ibd_GB"] = round(ibd.stat().st_size / 1e9, 3) if ibd.exists() else None

    hdr = read_header(path)
    s = summarize_header(hdr)
    row.update(s)
    dump = ctx.out("metadata", f"header_{path.stem}.txt")
    with open(dump, "w") as fh:
        for sec, params in hdr["sections"].items():
            fh.write(f"[{sec}]\n")
            for p in params:
                fh.write(f"  {p['accession']:<14} {p['name']}" + (f" = {p['value']}" if p["value"] else "")
                         + (f" ({p['unit']})" if p["unit"] else "") + "\n")
        fh.write("[first spectra]\n")
        for k, sp in enumerate(hdr["first_spectra"]):
            for p in sp:
                fh.write(f"  #{k} {p.get('where', ''):<30} {p['name']}" + (f" = {p['value']}" if p["value"] else "") + "\n")
    if not ibd.exists():
        return row

    t = time.time()
    r = SlideReader(path, log)
    row["parse_s"] = round(time.time() - t, 1)
    row["n_spectra"] = r.n
    row["grid_rows"], row["grid_cols"] = r.image_shape()
    row["x_min"], row["x_max"] = int(r.coords[:, 0].min()), int(r.coords[:, 0].max())
    row["y_min"], row["y_max"] = int(r.coords[:, 1].min()), int(r.coords[:, 1].max())
    row["fill_fraction"] = round(r.n / (row["grid_rows"] * row["grid_cols"]), 3)
    row["duplicate_coords"] = int(r.n - len(np.unique(r.coords, axis=0)))
    row["shares_mz_axis"] = r.shares_mz_axis
    row["memmap_ok"] = r._memmap is not None
    row["compressed"] = r.compressed
    row["mz_dtype"], row["int_dtype"] = str(r.mz_dtype), str(r.int_dtype)
    row["points_per_spectrum_median"] = int(np.median(r.int_lengths))

    sample = np.sort(rng.choice(r.n, min(200, r.n), replace=False))
    if r.shares_mz_axis:
        ax = r.mz_axis()
        row["mz_min"], row["mz_max"] = float(ax.min()), float(ax.max())
    else:
        lo, hi = [], []
        for i in sample:
            m = r.read_mz(int(i))
            if len(m):
                lo.append(m.min())
                hi.append(m.max())
        row["mz_min"], row["mz_max"] = (float(min(lo)), float(max(hi))) if lo else (np.nan, np.nan)

    # Is it really profile or centroid? Look at the 10 richest of 30 random spectra.
    probe = []
    for i in sample[:30]:
        mz, it = r.spectrum(int(i))
        probe.append((float(it.sum()), mz, it))
    probe.sort(key=lambda z: -z[0])
    guesses = [detect_profile(mz, it) for _, mz, it in probe[:10]]
    fr = np.nanmedian([g["shoulder_fraction"] for g in guesses])
    row["spectrum_type_data"] = "profile" if fr > 0.5 else "centroid"
    row["shoulder_fraction"] = round(float(fr), 3)
    row["median_spacing_ppm"] = round(float(np.nanmedian([g["median_spacing_ppm"] for g in guesses])), 2)
    row["_acquired"] = r.to_image(np.ones(r.n), fill=0)
    r.close()
    return row


def flags_for(row: dict, ctx, mode: str) -> list[str]:
    f = []
    d = ctx.mode_def(mode)
    if not row.get("name_matches_convention"):
        f.append("file name does not follow slide_<letter>_<mode>.imzML")
    if not row.get("ibd_found"):
        f.append("no .ibd file next to the .imzML")
        return f
    if row.get("polarity_header") and row["polarity_header"] != d["polarity"]:
        f.append(f"header polarity {row['polarity_header']} but folder says {d['polarity']}")
    if not row.get("polarity_header"):
        f.append("no polarity in the header")
    if row.get("spectrum_count_header") and row["spectrum_count_header"] != row.get("n_spectra"):
        f.append(f"spectrumList count {row['spectrum_count_header']} != spectra parsed {row.get('n_spectra')}")
    st_h, st_d = row.get("spectrum_type_header"), row.get("spectrum_type_data")
    if st_h and st_d and st_h != st_d:
        f.append(f"header says {st_h} spectra but the data look {st_d}; the data-based guess is used")
    if row.get("storage_mode") == "continuous" and not row.get("shares_mz_axis"):
        f.append("declared continuous but spectra do not share one m/z array")
    if row.get("storage_mode") == "processed" and row.get("shares_mz_axis"):
        f.append("declared processed but all spectra share one m/z array (fine, treated as continuous)")
    px = row.get("pixel_size_x_um")
    exp_px = ctx.get("acquisition.pixel_size_um")
    try:
        if px and exp_px and abs(float(px) - float(exp_px)) > 0.5:
            f.append(f"header pixel size {px} um, expected {exp_px} um")
    except ValueError:
        pass
    if not px:
        f.append("no pixel size in the header; using acquisition.pixel_size_um from the config")
    exp = d.get("expected_mz")
    if exp and np.isfinite(row.get("mz_min", np.nan)):
        if row["mz_min"] > exp[0] * 1.05 or row["mz_max"] < exp[1] * 0.95:
            f.append(f"m/z range {row['mz_min']:.1f} to {row['mz_max']:.1f} narrower than expected {exp}")
    if row.get("duplicate_coords"):
        f.append(f"{row['duplicate_coords']} duplicated pixel coordinates")
    if row.get("compressed"):
        f.append("zlib-compressed arrays: slower reading, no memmap fast path")
    return f


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("01_inventory", __doc__)
    rng = np.random.default_rng(ctx.get("runtime.seed", 0))
    sheet = samples.load_sheet(ctx)
    all_rows = []
    for mode in modes:
        S = Summary(ctx, "01_inventory", mode, f"Inventory: {mode}")
        folder = ctx.mode_dir(mode)
        found = sorted(p for p in folder.iterdir() if p.suffix.lower() == ".imzml")
        log.info(f"{mode}: {len(found)} imzML files in {folder}")
        rows = []
        for p in found:
            log.info(f"--- {p.name}")
            row = inspect_file(ctx, p, mode, log, rng)
            row["flags"] = " | ".join(flags_for(row, ctx, mode))
            rows.append(row)
            log.mem(f"after {p.name}")

        expected = samples.files_for_mode(ctx, mode, sheet)
        names = {r["file"] for r in rows}
        missing = [f"{s}: {p.name}" for s, p in expected.items() if p.name not in names]
        extra = sorted(names - {p.name for p in expected.values()})
        if missing:
            S.warn(f"files named in the sample sheet but not found: {missing}", log)
        if extra:
            S.warn(f"files in the folder that aren't in the sample sheet (ignored later): {extra}", log)

        df = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in rows])
        all_rows.append(df)

        S.h("Files")
        cols = ["file", "slide", "n_spectra", "grid_rows", "grid_cols", "fill_fraction", "ibd_GB",
                "storage_mode", "spectrum_type_header", "spectrum_type_data", "polarity_header",
                "mz_min", "mz_max", "median_spacing_ppm", "pixel_size_x_um"]
        S.table(df[[c for c in cols if c in df]])
        S.h("Instrument / software as written in the headers")
        S.p(f"Expected from the config: {ctx.get('acquisition.instrument_expected', 'not given')}")
        for _, r in df.iterrows():
            S.p(f"- **{r['file']}**: instrument: {r.get('instrument_params') or 'none listed'}; "
                f"software: {r.get('software') or 'none'}; processing: {r.get('data_processing') or 'none'}")
        S.h("Flags")
        any_flag = False
        for _, r in df.iterrows():
            if r["flags"]:
                any_flag = True
                S.p(f"- **{r['file']}**: {r['flags']}")
        if not any_flag:
            S.p("No flags.")
        types = df.get("spectrum_type_data", pd.Series(dtype=str)).dropna().unique().tolist()
        S.h("What the next phases will do")
        S.p(f"- spectrum type from the data: {types} -> "
            f"{'peak picking on mean spectra + window integration' if 'profile' in types else 'centroids used as they are (window sums)'}")
        S.p(f"- storage: {df['storage_mode'].dropna().unique().tolist()} "
            f"(memmap fast path: {df['memmap_ok'].tolist() if 'memmap_ok' in df else 'n/a'})")

        # quick look at which pixels were acquired on each slide
        acq = [r for r in rows if "_acquired" in r]
        if acq:
            fig, axes = plt.subplots(1, len(acq), figsize=(2.2 * len(acq), 3), squeeze=False)
            for ax, r in zip(axes[0], acq):
                ax.imshow(r["_acquired"], cmap="Greys", interpolation="nearest")
                ax.set_title(f"{r['slide']}  {r['n_spectra']} px", fontsize=8)
                ax.set_xticks([])
                ax.set_yticks([])
            fig.suptitle(f"{mode}: acquired pixels (black)", fontsize=9)
            S.fig(fig, "acquired_pixels.png")

        info = {"mode": mode, "spectrum_type": types[0] if len(types) == 1 else types,
                "storage": df["storage_mode"].dropna().unique().tolist(),
                "mz_min": float(df["mz_min"].min()) if "mz_min" in df else None,
                "mz_max": float(df["mz_max"].max()) if "mz_max" in df else None,
                "files": df["file"].tolist()}
        ctx.out("results", "01_inventory", mode, "inventory.json").write_text(json.dumps(info, indent=2, default=str))
        log.info(f"summary: {S.write()}")

    # files.csv collects every mode run so far (rerunning a mode replaces its rows)
    out = ctx.out("metadata", "files.csv")
    new = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    if out.exists() and len(new):
        old = pd.read_csv(out)
        new = pd.concat([old[~old["mode"].isin(modes)], new], ignore_index=True)
    new.to_csv(out, index=False)
    log.info(f"wrote {out}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
