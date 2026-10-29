"""Phase 05: name the peaks at sum-composition level and tag the hypothesis classes.

    python scripts/05_annotate.py --config config.yaml

Adducts come from annotation.adducts in the config ([M-H]- for negative mode; [M+H]+,
[M+Na]+, [M+K]+ for positive), plus a few class-specific ones (e.g. [M-2H]2- for
gangliosides), which are ranked below the default adducts.

ppm tolerance, in order of preference:
  1. annotation.ppm_tolerance if you set a number
  2. 03's suggestion (after recalibration if 04 recalibrated, raw otherwise)
  3. 02's suggestion
  4. annotation.ppm_fallback

Writes processed/<mode>/annotations.csv (one row per peak in matrix_raw.h5ad),
results/05_annotate/<mode>/hypothesis_peaks.csv, figures, summary.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import annotate, cli, refdb  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def pick_tolerance(ctx, mode: str, recalibrated: bool, log) -> tuple[float, str]:
    A = ctx.cfg["annotation"]
    if A.get("ppm_tolerance"):
        return float(A["ppm_tolerance"]), "annotation.ppm_tolerance (config)"
    m3 = ctx.phase_json("03_qc_all", mode, "mass_accuracy.json")
    key = "tol_recalibrated_ppm" if recalibrated else "tol_raw_ppm"
    if m3 and m3.get(key) and np.isfinite(m3[key]):
        return float(m3[key]), f"03 mass accuracy ({key})"
    m2 = ctx.phase_json("02_qc", mode, "mass_accuracy.json")
    if m2 and m2.get("tol_raw_ppm") and np.isfinite(m2["tol_raw_ppm"]):
        return float(m2["tol_raw_ppm"]), "02 mass accuracy (tol_raw_ppm, one slide)"
    log.warning("no mass-accuracy estimate found, using annotation.ppm_fallback")
    return float(A.get("ppm_fallback", 5)), "annotation.ppm_fallback"


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("05_annotate", __doc__)
    A = ctx.cfg["annotation"]
    lmsd = refdb.cached_table(ctx, "lipidmaps")
    if lmsd is None and A.get("lipidmaps_filter", True):
        log.warning("no LIPID MAPS cache (run 00_setup.py --download): compositions are NOT filtered by LMSD")
    lib_lip = annotate.build_lipid_library(lmsd, log)

    for mode in modes:
        S = Summary(ctx, "05_annotate", mode, f"Annotation, {mode}")
        p = ctx.out("processed", mode, "matrix_raw.h5ad", mkdir=False)
        if not p.exists():
            S.warn("no matrix_raw.h5ad, run 04 first", log)
            S.write()
            continue
        a = ad.read_h5ad(p, backed="r")     # only var/uns needed, so don't load X
        var = a.var.copy()
        recal = bool(a.uns.get("recalibrated", False))
        a.file.close()

        d = ctx.mode_def(mode)
        lib = lib_lip
        if d["analyte"] == "metabolite":
            lib = pd.concat([lib_lip, annotate.build_metabolite_library(ctx, log)], ignore_index=True)
        ions = annotate.ion_table(lib, d["polarity"], A["adducts"])
        tol, src = pick_tolerance(ctx, mode, recal, log)
        log.info(f"{mode}: {len(var)} peaks, {len(ions)} candidate ions, tolerance {tol:.2f} ppm ({src})")

        peaks = var[["mz"]].copy()
        ann = annotate.match_peaks(peaks, ions, tol, bool(A.get("lipidmaps_filter", True)) and lmsd is not None,
                                   int(A.get("max_candidates_reported", 5)))
        iso = annotate.isotope_flags(var["mz"].to_numpy(), var["median_height"].to_numpy(), max(tol, 3.0))
        iso.index = var.index
        out = pd.concat([var[["mz", "median_height", "n_sections"]], ann, iso], axis=1)
        out["hypothesis"] = [annotate.hypothesis_direction(t, ctx.cfg.get("hypotheses", {})) for t in out["tag"]]
        # Names for the plots. A peak that looks like an isotope keeps its candidate list, but
        # it's labelled after its monoisotopic peak (we follow M+2 -> M+1 -> M back down).
        base = np.where(out["annotation"] != "", out["annotation"] + " " + out["adduct"],
                        [f"mz {m:.4f}" for m in out["mz"]])
        base = pd.Series(base, index=out.index)
        base[out["ambiguous"]] += " *"
        by_mz = dict(zip(out["mz"].round(6), out.index))
        labels = []
        for pid, row in out.iterrows():
            k, cur = 0, row
            while cur["is_isotope"] and k < 5:
                k += 1
                cur = out.loc[by_mz[round(cur["isotope_of_mz"], 6)]]
            labels.append(f"M+{k} of {base[cur.name]}" if k else base[pid])
        out["label"] = labels
        out["tolerance_ppm"] = tol
        out.index.name = "peak_id"
        out.to_csv(ctx.out("processed", mode, "annotations.csv"))

        hyp = out[(out["hypothesis"] != "") & ~out["is_isotope"]].sort_values(["tag", "mz"])
        hyp.to_csv(ctx.out("results", "05_annotate", mode, "hypothesis_peaks.csv"))

        # figures
        named = out[(out["annotation"] != "") & ~out["is_isotope"]]
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
        axes[0].hist(named["ppm_error"], bins=30, color="0.5")
        axes[0].axvline(0, color="k", lw=0.5)
        axes[0].set_xlabel("ppm error of the best match")
        axes[0].set_title(f"{len(named)} annotated peaks (tolerance {tol:.1f} ppm)", fontsize=9)
        md = out["mz"] - np.floor(out["mz"])
        for cls, sub in out.assign(md=md).groupby(out["lipid_class"].replace("", "unannotated")):
            axes[1].scatter(sub["mz"], sub["md"], s=10, label=cls)
        axes[1].set_xlabel("m/z")
        axes[1].set_ylabel("mass defect")
        axes[1].set_title("Kendrick-style overview (colour = class)", fontsize=9)
        axes[1].legend(fontsize=6, ncol=2, loc="upper left")
        S.fig(fig, "annotation_overview.png")

        S.h("Settings")
        S.kv({"tolerance (ppm)": tol, "tolerance source": src, "adducts": ", ".join(A["adducts"][d["polarity"]]),
              "LIPID MAPS filter": "on" if (lmsd is not None and A.get("lipidmaps_filter", True)) else "off (no cache)",
              "candidate ions": len(ions)})
        S.h("Counts")
        S.kv({"peaks": len(out), "annotated (best match)": int((out["annotation"] != "").sum()),
              "ambiguous (>1 composition within tolerance)": int(out["ambiguous"].sum()),
              "look like M+1 isotopes": int(out["is_isotope"].sum()),
              "type II overlap risk (>10% of signal)": int(out["typeII_flag"].sum()),
              "hypothesis-class peaks": len(hyp)})
        S.h("Hypothesis classes found")
        S.table(hyp[["mz", "label", "tag", "hypothesis", "ppm_error", "n_candidates", "typeII_flag", "candidates"]]
                .round(4), max_rows=60)
        S.h("All peaks")
        S.table(out[["mz", "label", "lipid_class", "ppm_error", "is_isotope", "typeII_fraction", "candidates"]]
                .round(4), max_rows=80)
        S.p("\n`*` = ambiguous: more than one composition fits; see `candidates`. BMP and PG have identical "
            "formulas, so they're always reported together as BMP/PG. GalCer and GlcCer are HexCer.")
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
