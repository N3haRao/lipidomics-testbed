"""Phase 04: one peak list for all slides, then a pixels x peaks matrix per mode.

    python scripts/04_peaks_matrix.py --config config.yaml
    python scripts/04_peaks_matrix.py --config config.yaml --mode lip_neg

Needs 03 (it reuses 03's section labels, mean spectra and m/z offsets).

Steps
  1. Pick peaks on each SECTION's mean spectrum (12 per mode). If the data are profile,
     this is the centroiding step. If they're already centroided, the mean spectrum is a
     binned histogram of the centroids, and picking finds the consensus positions.
  2. Shift every slide's peaks by its own m/z offset from 03 (recalibration), then group
     peaks across sections when they're within align_ppm of each other.
  3. Keep a consensus peak if it was picked in >= min_section_fraction of the sections.
  4. Stream each slide once more and integrate every kept peak in every tissue pixel
     (sum of intensities inside the peak's window; the window moves with the slide offset).
  5. Drop peaks present (above local background) in < min_pixel_fraction of the pixels of
     most sections. DAN matrix peaks are moved out of X into obs columns (dan_...).

Writes processed/<mode>/matrix_raw.h5ad (X = raw integrated intensities, float32),
       processed/<mode>/peaks.csv, results/04_peaks_matrix/<mode>/..., figures, summary.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, qc, samples, spectra  # noqa: E402
from cln3msi.imzml_io import SlideReader  # noqa: E402
from cln3msi.report import Summary, decimate, plt  # noqa: E402
from cln3msi.tissue import BOTTOM, POSITION_NAME, TOP  # noqa: E402


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("04_peaks_matrix", __doc__)
    sheet = samples.load_sheet(ctx)
    P = ctx.cfg["peaks"]
    for mode in modes:
        S = Summary(ctx, "04_peaks_matrix", mode, f"Peak list and matrix, {mode}")
        files = samples.files_for_mode(ctx, mode, sheet)
        slides = [s for s in sorted(files) if ctx.out("processed", mode, "qc", f"meanspec_{s}.npz", mkdir=False).exists()]
        if len(slides) < len(files):
            S.warn(f"no 03 QC output for slides {sorted(set(files) - set(slides))}: run 03 first. Skipping them.", log)
        if not slides:
            S.write()
            continue

        ma = ctx.phase_json("03_qc_all", mode, "mass_accuracy.json") or {}
        offsets = {s: float(ma.get("per_slide_offset_ppm", {}).get(s, 0.0)) for s in slides} \
            if P.get("recalibrate_per_slide", True) else {s: 0.0 for s in slides}
        if P.get("recalibrate_per_slide", True) and not ma:
            S.warn("no 03 mass_accuracy.json: no recalibration applied", log)
        floor = float(ctx.get("qc.ppm_floor", 2))
        if P.get("align_ppm"):
            align = float(P["align_ppm"])
        else:
            sd = ma.get("robust_sd_recal_ppm") if P.get("recalibrate_per_slide", True) else ma.get("robust_sd_raw_ppm")
            align = max(floor, 3 * sd) if sd and np.isfinite(sd) else float(ctx.get("annotation.ppm_fallback", 5))

        # ---- 1. per-section peak picking ------------------------------------------------
        tables, meta = [], []
        stype = None
        for s in slides:
            z = np.load(ctx.out("processed", mode, "qc", f"meanspec_{s}.npz", mkdir=False))
            stype = str(z["spectrum_type"])
            binned = bool(z["binned"])
            if binned:
                # binned grid: width of one bin limits how precisely sections can agree
                bin_ppm = float(np.median(np.diff(np.log(z["axis"])))) * 1e6
                align = max(align, 2 * bin_ppm)
            for g in ("Top", "Bottom"):
                if g not in z:
                    continue
                pk = spectra.pick_peaks(z["axis"], z[g], snr=float(P["snr"]), noise_window=int(P["noise_window"]),
                                        smooth_bins=1.0 if (binned and stype == "centroid") else 0.0)
                mz_true = pk["mz"] / (1 + offsets[s] * 1e-6)
                tables.append(mz_true)
                meta.append({"slide": s, "section": g, "height": pk["height"], "snr": pk["snr"],
                             "fwhm_ppm": pk["fwhm_ppm"]})
                log.info(f"{s} {g}: {len(mz_true)} peaks (SNR >= {P['snr']})")
        n_sec = len(tables)
        log.info(f"align tolerance {align:.2f} ppm; spectrum type {stype}; {n_sec} sections")

        # ---- 2-3. alignment and section filter ------------------------------------------
        clusters = spectra.align_peaks(tables, align)
        rows = []
        for cl in clusters:
            secs = {k for k, _ in cl}
            mzs = np.array([tables[k][j] for k, j in cl])
            h = np.array([meta[k]["height"][j] for k, j in cl])
            rows.append({"mz": float(np.average(mzs, weights=np.clip(h, 1e-12, None))),
                         "n_sections": len(secs), "spread_ppm": float(np.ptp(mzs) / mzs.mean() * 1e6),
                         "median_height": float(np.median(h)),
                         "median_snr": float(np.median([meta[k]["snr"][j] for k, j in cl])),
                         "median_fwhm_ppm": float(np.median([meta[k]["fwhm_ppm"][j] for k, j in cl]))})
        cons = pd.DataFrame(rows)
        n_all = len(cons)
        need = int(np.ceil(float(P["min_section_fraction"]) * n_sec))
        cons = cons[cons.n_sections >= need].sort_values("mz").reset_index(drop=True)
        log.info(f"{n_all} aligned features, {len(cons)} found in >= {need}/{n_sec} sections")
        if len(cons) > int(P["max_peaks"]):
            cons = cons.nlargest(int(P["max_peaks"]), "median_height").sort_values("mz").reset_index(drop=True)
        if cons.empty:
            S.warn("no peaks survived; lower peaks.snr or peaks.min_section_fraction", log)
            S.write()
            continue

        if stype == "profile":
            half = float(P["integration_fwhm_factor"]) * cons["median_fwhm_ppm"].to_numpy()
        else:
            half = np.full(len(cons), align)
        lo, hi = spectra.ppm_windows(cons["mz"].to_numpy(), half)
        lo, hi = spectra.clip_windows(cons["mz"].to_numpy(), lo, hi)
        cons["mz_lo"], cons["mz_hi"] = lo, hi

        d = ctx.mode_def(mode)
        mx = qc.matrix_ions(ctx, d["polarity"], 0, 1e9)
        cons["is_matrix"] = False
        cons["matrix_ion"] = ""
        for _, m in mx.iterrows():
            close = (np.abs(cons["mz"] - m.mz) / m.mz * 1e6) <= max(align, 5)
            cons.loc[close, "is_matrix"] = True
            cons.loc[close, "matrix_ion"] = m["name"]

        # ---- 4. extraction --------------------------------------------------------------
        X_parts, obs_parts, present = [], [], {}
        chunk = int(ctx.get("runtime.spectra_chunk", 2000))
        nw = int(P["noise_window"])
        for s in slides:
            z = np.load(ctx.out("processed", mode, "qc", f"sections_{s}.npz", mkdir=False))
            ms = np.load(ctx.out("processed", mode, "qc", f"meanspec_{s}.npz", mkdir=False))
            pos_spec, core_spec = z["pos_spec"], z["core_spec"]
            tis = np.flatnonzero(pos_spec > 0)
            f = 1 + offsets[s] * 1e-6
            slo, shi = lo * f, hi * f
            r = SlideReader(files[s], log)
            vals = np.zeros((len(tis), len(cons)), dtype=np.float32)
            if r.shares_mz_axis:
                axis = r.mz_axis()
                k = 0
                for idx, block in r.iter_blocks(indices=tis, chunk=chunk):
                    vals[k:k + len(idx)] = spectra.window_sums_block(axis, block, slo, shi)
                    k += len(idx)
            else:
                row_of = {int(i): j for j, i in enumerate(tis)}
                for i, mzv, it in r.iter_spectra(indices=tis, log=log, every=int(ctx.get("runtime.log_memory_every", 10000))):
                    vals[row_of[i]] = spectra.window_sums(mzv, it, slo, shi)
            r.close()
            log.mem(f"extracted slide {s}: {vals.shape}")
            # background in each window from the section mean spectrum's baseline
            for code, g in ((TOP, "Top"), (BOTTOM, "Bottom")):
                sel = pos_spec[tis] == code
                if not sel.any() or g not in ms:
                    continue
                base, _ = spectra.noise_level(ms[g], nw)
                bg = spectra.window_sums(ms["axis"], base, slo, shi)
                present[(s, g)] = (vals[sel] > 2 * bg + 1e-9).mean(axis=0)
            sub = sheet[sheet.slide == s].set_index("position_on_slide")
            posname = np.array([POSITION_NAME.get(int(c), "") for c in pos_spec[tis]])
            o = pd.DataFrame({
                "slide": s, "section": posname,
                "animal_id": [sub.loc[p, "animal_id"] for p in posname],
                "genotype": [sub.loc[p, "genotype"] for p in posname],
                "sex": sub["sex"].iloc[0], "run_order": int(sub["run_order"].iloc[0]),
                "x": z["coords"][tis, 0], "y": z["coords"][tis, 1], "core": core_spec[tis],
                "tic_spectrum": z["tic"][tis], "dan_spectrum": z["matrix"][tis],
            })
            X_parts.append(vals)
            obs_parts.append(o)

        X = np.concatenate(X_parts)
        obs = pd.concat(obs_parts, ignore_index=True)
        obs.index = [f"{a}_{x}_{y}" for a, x, y in zip(obs.animal_id, obs.x, obs.y)]
        del X_parts
        log.mem(f"matrix assembled {X.shape}")

        # ---- 5. pixel-presence filter + matrix peaks out of X ---------------------------
        pres = pd.DataFrame(present).T              # (slide, section) x peaks
        frac_ok = (pres >= float(P["min_pixel_fraction"])).mean(axis=0).to_numpy()
        cons["pixel_presence_median"] = pres.median(axis=0).to_numpy()
        cons["sections_present_fraction"] = frac_ok
        keep = frac_ok >= float(P["min_section_fraction"])
        n_before = len(cons)
        if P.get("drop_matrix_peaks", True):
            for j in np.flatnonzero(cons["is_matrix"].to_numpy()):
                obs[f"dan_{cons.loc[j, 'mz']:.4f}"] = X[:, j]
            keep &= ~cons["is_matrix"].to_numpy()
        cons_kept = cons[keep].reset_index(drop=True)
        X = X[:, keep]
        log.info(f"pixel filter + matrix removal: {n_before} -> {len(cons_kept)} peaks")

        var = cons_kept.copy()
        var.index = [f"mz_{m:.4f}" for m in var["mz"]]
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["mode"] = mode
        adata.uns["spectrum_type"] = stype
        adata.uns["align_ppm"] = align
        adata.uns["recalibration_offsets_ppm"] = json.dumps(offsets)
        adata.uns["recalibrated"] = bool(P.get("recalibrate_per_slide", True) and ma)
        out = ctx.out("processed", mode, "matrix_raw.h5ad")
        adata.write_h5ad(out, compression="gzip")
        cons.to_csv(ctx.out("processed", mode, "peaks_all_candidates.csv"), index=False)
        var.to_csv(ctx.out("processed", mode, "peaks.csv"))
        log.info(f"wrote {out}")

        # ---- figures ----------------------------------------------------------------------
        z = np.load(ctx.out("processed", mode, "qc", f"meanspec_{slides[0]}.npz", mkdir=False))
        fig, ax = plt.subplots(figsize=(12, 3.8))
        tis_mean = np.mean([z[g] for g in ("Top", "Bottom") if g in z], axis=0)
        ax.plot(*decimate(z["axis"] / (1 + offsets[slides[0]] * 1e-6), tis_mean), lw=0.5, color="0.3")
        ymax = tis_mean.max()
        ax.vlines(var["mz"], 0, ymax * 0.05, color="tab:red", lw=0.6, label="kept peaks")
        dropped = cons[~keep]
        ax.vlines(dropped["mz"], 0, ymax * 0.03, color="tab:orange", lw=0.6, label="dropped (pixel filter / matrix)")
        ax.set_xlabel("m/z (recalibrated)")
        ax.set_title(f"{mode}: {len(var)} kept peaks on slide {slides[0]}'s tissue mean spectrum", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        S.fig(fig, "peaks_on_mean_spectrum.png")

        fig, axes = plt.subplots(1, 2, figsize=(9, 3))
        axes[0].hist(pd.DataFrame(rows)["n_sections"], bins=np.arange(0.5, n_sec + 1.5), color="0.5")
        axes[0].axvline(need - 0.5, color="tab:red")
        axes[0].set_xlabel("sections a feature was picked in")
        axes[1].hist(var["spread_ppm"], bins=30, color="0.5")
        axes[1].set_xlabel("m/z spread across sections (ppm), kept peaks")
        S.fig(fig, "alignment_qc.png")

        top = np.argsort(-np.asarray(X.mean(axis=0)).ravel())[:4]
        fig, axes = plt.subplots(len(top), len(slides), figsize=(1.9 * len(slides), 2.2 * len(top)), squeeze=False)
        for a, j in enumerate(top):
            vmax = np.percentile(X[:, j], 99)
            for b, s in enumerate(slides):
                sel = (obs.slide == s).to_numpy()
                yy, xx = obs.y.to_numpy()[sel], obs.x.to_numpy()[sel]
                img = np.full((yy.max() + 1, xx.max() + 1), np.nan)
                img[yy, xx] = X[sel, j]
                axes[a, b].imshow(img, vmin=0, vmax=vmax, cmap=ctx.get("report.cmap", "viridis"), interpolation="nearest")
                axes[a, b].set_xticks([]); axes[a, b].set_yticks([])  # noqa: E702
                if a == 0:
                    axes[a, b].set_title(s, fontsize=8)
            axes[a, 0].set_ylabel(f"{var['mz'].iloc[j]:.4f}", fontsize=8)
        fig.suptitle(f"{mode}: 4 strongest peaks, raw intensity, shared colour scale per row", fontsize=9)
        S.fig(fig, "top_peaks_images.png")

        S.h("Peak list")
        S.kv({"spectrum type": stype, "alignment tolerance (ppm)": align,
              "recalibration offsets (ppm)": ", ".join(f"{k} {v:+.2f}" for k, v in offsets.items()),
              "sections": n_sec, "aligned features": n_all,
              f"in >= {need} sections": n_before, "after pixel filter and matrix removal": len(var),
              "matrix (DAN) peaks moved to obs": int(cons["is_matrix"].sum())})
        S.h("Matrix")
        S.kv({"pixels": adata.n_obs, "peaks": adata.n_vars,
              "pixels per section (median)": int(obs.groupby(["slide", "section"]).size().median()),
              "file": out.relative_to(ctx.output_root)})
        S.table(obs.groupby(["slide", "section", "animal_id", "genotype"]).size().rename("pixels").reset_index())
        S.h("Strongest peaks")
        S.table(var.nlargest(15, "median_height")[["mz", "n_sections", "spread_ppm", "median_snr",
                                                    "median_fwhm_ppm", "pixel_presence_median"]].round(4))
        log.info(f"summary: {S.write()}")
        del adata, X
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
