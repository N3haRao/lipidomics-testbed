"""Phase 03: the 02 QC on every slide, plus cross-slide comparisons and the gradient check.

    python scripts/03_qc_all.py --config config.yaml                  # every mode in modes.run
    python scripts/03_qc_all.py --config config.yaml --mode lip_neg

Per mode it writes
  results/03_qc_all/<mode>/qc_<slide>.json, calibrants_<slide>.csv    (same as 02, every slide)
  results/03_qc_all/<mode>/mass_accuracy.json   per-slide m/z offset (used by 04 to recalibrate)
                                                and the suggested annotation tolerance (used by 05)
  results/03_qc_all/<mode>/gradient.csv         top-to-bottom trend per slide (used by 09)
  results/03_qc_all/<mode>/sections.csv         pixels / area / TIC per section
  processed/<mode>/qc/*.npz                     section labels + mean spectra (used by 04)
  figures/03_qc_all/<mode>/...                  incl. sections_all.png, the one to check first

About the gradient check: Control is always on top and Mutant always on the bottom, so
any technical top-to-bottom trend (matrix coating, laser drift during the run, ...) looks
exactly like a genotype effect. We can't correct for that. What we can do is measure the
trend where there's no biology (the off-tissue matrix signal) and flag hits that are not
clearly bigger than it.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, qc, samples  # noqa: E402
from cln3msi.report import Summary, decimate, plt  # noqa: E402
from cln3msi.tissue import BOTTOM, TOP  # noqa: E402


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("03_qc_all", __doc__)
    sheet = samples.load_sheet(ctx)
    floor = float(ctx.get("qc.ppm_floor", 2))
    min_snr = float(ctx.get("qc.calibrant_min_snr", 10))
    for mode in modes:
        S = Summary(ctx, "03_qc_all", mode, f"QC of all slides, {mode}")
        files = samples.files_for_mode(ctx, mode, sheet)
        results = {}
        for slide in sorted(files):
            if not files[slide].exists():
                S.warn(f"slide {slide}: {files[slide].name} not found, skipped", log)
                continue
            results[slide] = qc.run_slide_qc(ctx, log, mode, slide, "03_qc_all", sheet, S, make_figures=True)
            log.mem(f"after slide {slide}")
        if not results:
            S.warn("no slides processed", log)
            S.write()
            continue

        # ---- mass accuracy across slides ----------------------------------------------------
        cal = pd.concat([r["_cal"] for r in results.values()], ignore_index=True)
        cal.to_csv(ctx.out("results", "03_qc_all", mode, "calibrants_all.csv"), index=False)
        good = cal[cal.snr >= min_snr].copy()
        per_slide = good.groupby("slide")["ppm_error"].agg(["median", "count"]).rename(
            columns={"median": "offset_ppm", "count": "n"})
        good["recal_error"] = good["ppm_error"] - good["slide"].map(per_slide["offset_ppm"])
        raw = qc.suggest_tolerance(good["ppm_error"].to_numpy(), floor)
        rec = qc.suggest_tolerance(good["recal_error"].to_numpy(), floor)
        ma = {"per_slide_offset_ppm": per_slide["offset_ppm"].round(3).to_dict(),
              "n_per_slide": per_slide["n"].to_dict(),
              "tol_raw_ppm": raw["tol_raw_ppm"], "median_ppm": raw["median_ppm"],
              "robust_sd_raw_ppm": raw["robust_sd_ppm"],
              "tol_recalibrated_ppm": rec["tol_recalibrated_ppm"], "robust_sd_recal_ppm": rec["robust_sd_ppm"]}
        ctx.out("results", "03_qc_all", mode, "mass_accuracy.json").write_text(json.dumps(ma, indent=2))

        # ---- sections table ---------------------------------------------------------------
        secs = pd.DataFrame([s for r in results.values() for s in r["sections"]])
        secs.to_csv(ctx.out("results", "03_qc_all", mode, "sections.csv"), index=False)

        # ---- gradient table ----------------------------------------------------------------
        grad = pd.DataFrame([{"slide": s, **r["gradient"]} for s, r in results.items()])
        grad.to_csv(ctx.out("results", "03_qc_all", mode, "gradient.csv"), index=False)

        # ---- figures ---------------------------------------------------------------------
        slides = list(results)
        fig, axes = plt.subplots(1, len(slides), figsize=(2.6 * len(slides), 3.6), squeeze=False)
        for ax, s in zip(axes[0], slides):
            z = np.load(ctx.out("processed", mode, "qc", f"sections_{s}.npz", mkdir=False))
            img = np.zeros(z["positions"].shape)
            img[z["positions"] == TOP] = 1
            img[z["positions"] == BOTTOM] = 2
            ax.imshow(img, cmap="coolwarm", vmin=0, vmax=2, interpolation="nearest")
            sub = sheet[sheet.slide == s].set_index("position_on_slide")
            for code, posname in ((TOP, "Top"), (BOTTOM, "Bottom")):
                yy, xx = np.nonzero(z["positions"] == code)
                if len(yy):
                    g = sub.loc[posname, "genotype"]
                    lab = ctx.cfg["project"]["mutant_label"] if g == ctx.cfg["project"]["genotypes"]["mutant"] else g
                    ax.text(xx.mean(), yy.mean(), f"{sub.loc[posname, 'animal_id']}\n{lab}", ha="center",
                            va="center", fontsize=7, bbox=dict(facecolor="white", alpha=0.8, lw=0))
            ax.set_title(f"slide {s}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
        fig.suptitle(f"{mode}: section assignment on every slide (blue = Top, red = Bottom)", fontsize=10)
        S.fig(fig, "sections_all.png")

        fig, ax = plt.subplots(figsize=(7, 3.5))
        for s in slides:
            sub = good[good.slide == s]
            ax.scatter(sub.theoretical_mz, sub.ppm_error, s=14, label=s)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("m/z"); ax.set_ylabel("ppm error")  # noqa: E702
        ax.set_title(f"{mode}: calibrant mass error per slide", fontsize=9)
        ax.legend(fontsize=7, ncol=3, loc="upper right")
        S.fig(fig, "mass_error_all.png")

        fig, axes = plt.subplots(len(slides), 1, figsize=(12, 1.4 * len(slides)), sharex=True, squeeze=False)
        for ax, s in zip(axes[:, 0], slides):
            z = np.load(ctx.out("processed", mode, "qc", f"meanspec_{s}.npz", mkdir=False))
            for g, col in (("Top", "tab:blue"), ("Bottom", "tab:red")):
                if g in z:
                    sp = z[g] / max(z[g].sum(), 1e-12)     # TIC-scaled so slides are comparable
                    ax.plot(*decimate(z["axis"], sp), lw=0.5, color=col)
            ax.set_ylabel(s, rotation=0, labelpad=12)
            ax.set_yticks([])
        axes[-1, 0].set_xlabel("m/z")
        fig.suptitle(f"{mode}: TIC-scaled mean spectra per slide (blue Top, red Bottom)", fontsize=10)
        S.fig(fig, "mean_spectra_all.png")

        fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
        secs["label"] = secs["slide"] + " " + secs["position"].str[0]
        axes[0].bar(secs["label"], secs["median_tic"],
                    color=["tab:blue" if p == "Top" else "tab:red" for p in secs["position"]])
        axes[0].set_title("median TIC per section", fontsize=9)
        axes[0].tick_params(axis="x", rotation=90, labelsize=7)
        col = "matrix_off_pred_log2_bottom_vs_top"
        if col in grad:
            axes[1].bar(grad["slide"], grad[col], color="0.4", label="off-tissue DAN trend")
        if "tic_off_pred_log2_bottom_vs_top" in grad:
            axes[1].scatter(grad["slide"], grad["tic_off_pred_log2_bottom_vs_top"], color="k", zorder=3,
                            label="off-tissue TIC trend")
        axes[1].axhline(0, color="k", lw=0.5)
        axes[1].set_ylabel("predicted log2(Bottom / Top)")
        axes[1].set_title("top-to-bottom technical trend, per slide", fontsize=9)
        axes[1].legend(fontsize=7, loc="upper right")
        S.fig(fig, "gradient_all.png")

        # ---- summary --------------------------------------------------------------------
        S.h("Slides")
        rows = [{"slide": s, "spectra": r["n_spectra"], "type": r["spectrum_type"], "storage": r["storage"],
                 "split": r["section_split"], "tissue threshold": r.get("tissue_threshold_method", ""),
                 "off-tissue px": r["n_off_tissue"],
                 "median ppm": r["mass_accuracy"]["median_ppm"], "calibrants": r["mass_accuracy"]["n_calibrant_hits"],
                 "R median": r["resolution"].get("R_median", np.nan)} for s, r in results.items()]
        S.table(pd.DataFrame(rows))
        S.h("Sections (check sections_all.png first)")
        S.table(secs.drop(columns=["label"]))
        mut = ctx.cfg["project"]["genotypes"]["mutant"]
        ctl = ctx.cfg["project"]["genotypes"]["control"]
        piv = secs.pivot_table(index="slide", columns="genotype", values="area_mm2")
        if {mut, ctl} <= set(piv.columns):
            ratio = (piv[mut] / piv[ctl])
            S.p(f"Section area {ctx.cfg['project']['mutant_label']} / Control per slide: "
                + ", ".join(f"{s} {v:.2f}" for s, v in ratio.items())
                + f" (median {ratio.median():.2f}). Descriptive only: area depends on cutting and tissue "
                "detection too.")
        S.h("Mass accuracy")
        S.kv({"median error, all slides (ppm)": ma["median_ppm"], "robust SD raw (ppm)": ma["robust_sd_raw_ppm"],
              "tolerance raw (ppm)": ma["tol_raw_ppm"],
              "robust SD after per-slide offset (ppm)": ma["robust_sd_recal_ppm"],
              "tolerance after recalibration (ppm)": ma["tol_recalibrated_ppm"]})
        S.p("Per-slide offsets (ppm), applied in 04 when peaks.recalibrate_per_slide is true: "
            + ", ".join(f"{k} {v:+.2f}" for k, v in ma["per_slide_offset_ppm"].items()))
        S.h("Top-to-bottom gradient (position confound)")
        S.p("Control is always on top, so a technical trend down the slide looks like a genotype effect. "
            "Numbers are log2(Bottom/Top) predicted from the off-tissue trend; 0.1 means ~7%.")
        show = [c for c in grad.columns if c.endswith(("pred_log2_bottom_vs_top", "spearman"))] + ["n_off_pixels"]
        S.table(grad[["slide"] + show])
        big = grad[grad.get("matrix_off_pred_log2_bottom_vs_top", pd.Series(0, index=grad.index)).abs() > 0.2]
        if len(big):
            S.warn(f"slides with an off-tissue DAN trend above 0.2 log2 units Bottom vs Top: "
                   f"{big.slide.tolist()}. Hits of similar size on those slides will be flagged in 09.", log)
        if (grad["n_off_pixels"] < 20).any():
            S.warn("some slides have (almost) no off-tissue pixels, so the off-tissue gradient can't be "
                   "measured there; only the on-tissue TIC trend is reported.", log)
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
