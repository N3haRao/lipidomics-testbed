"""Phase 02: careful QC of ONE slide and mode, before committing to the full run.

    python scripts/02_qc_one_slide.py --config config.yaml                       # qc.slide / qc.mode
    python scripts/02_qc_one_slide.py --config config.yaml --slide F --mode lip_pos

What it does (see src/cln3msi/qc.py for the details):
  * TIC image, image without the DAN matrix ions, image of the matrix ions
  * tissue vs background, split into Top / Bottom sections, labelled figure to confirm
    which section is which animal
  * mean spectra of Top, Bottom and off-tissue
  * mass accuracy from calibrant lipids whose m/z are computed from their formulas
  * which DAN matrix ions are present (and how much stronger off tissue)
  * resolving power from peak widths in single spectra, if the spectra are profile
  * the suggested ppm tolerance for annotation (written to mass_accuracy.json)

Writes results/02_qc/<mode>/, figures/02_qc/<mode>/, summaries/02_qc_<mode>.md
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, qc, samples  # noqa: E402
from cln3msi.report import Summary  # noqa: E402


def main():
    t0 = time.time()
    ctx, log, modes, args = cli.start("02_qc", __doc__, slide=True)
    mode = args.mode[0] if args.mode else ctx.get("qc.mode")
    slide = (args.slide or ctx.get("qc.slide")).upper()
    sheet = samples.load_sheet(ctx)
    if slide not in set(sheet.slide):
        log.error(f"slide {slide} is not in the sample sheet ({sorted(sheet.slide.unique())})")
        sys.exit(2)

    S = Summary(ctx, "02_qc", mode, f"QC of slide {slide}, {mode}")
    res = qc.run_slide_qc(ctx, log, mode, slide, "02_qc", sheet, S)

    S.h("Data")
    S.kv({"file": res["file"], "spectra": res["n_spectra"], "spectrum type (from data)": res["spectrum_type"],
          "storage": res["storage"], "m/z range": f"{res['mz_min']:.2f} to {res['mz_max']:.2f}",
          "off-tissue pixels": res["n_off_tissue"]})
    S.h("Sections")
    S.p(f"Split: {res['section_split']} (split_axis = {res['split_axis']}). Check `sections_{slide}.png`: "
        "the Control animal must carry the 'Top' label. If Control and KO are swapped, set "
        "qc.sections.top_is_low_index: false (all slides) or list the slide in flip_slides.")
    import pandas as pd
    S.table(pd.DataFrame(res["sections"]))

    S.h("Mass accuracy (calibrants)")
    cal = res["_cal"]
    S.table(cal[["section", "name", "adduct", "theoretical_mz", "mz", "ppm_error", "snr"]].round(4))
    ma = res["mass_accuracy"]
    good = cal[cal.snr >= float(ctx.get("qc.calibrant_min_snr", 10))]
    tol = qc.suggest_tolerance(good["ppm_error"].to_numpy(), float(ctx.get("qc.ppm_floor", 2)))
    S.kv({"calibrant hits with SNR >= min": ma["n_calibrant_hits"], "median error (ppm)": ma["median_ppm"],
          "robust SD (ppm)": tol["robust_sd_ppm"],
          "suggested tolerance, raw (ppm)": tol["tol_raw_ppm"],
          "suggested tolerance after per-slide recalibration (ppm)": tol["tol_recalibrated_ppm"]})
    if ma["n_calibrant_hits"] < 3:
        S.warn("fewer than 3 calibrants found: tolerance can't be estimated from this slide; "
               "annotation will use annotation.ppm_fallback unless 03 finds more.", log)
    ctx.out("results", "02_qc", mode, "mass_accuracy.json").write_text(json.dumps(tol | {"slide": slide}, indent=2))

    S.h("DAN matrix ions")
    S.table(res["_mx"].round(3))
    S.h("Resolving power")
    S.kv(res["resolution"])
    S.h("Top-to-bottom gradient")
    S.kv({k: v for k, v in res["gradient"].items()})
    S.p("`*_pred_log2_bottom_vs_top` = how much brighter (log2) the Bottom section would look than the Top one "
        "if the off-tissue trend applied to tissue. Phase 09 uses it to flag hits.")
    log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
