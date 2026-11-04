"""Phase 09: which lipids change in CLN3 KO, and where. Per-animal summaries, paired by slide.

    python scripts/09_differential.py --config config.yaml

Units we test in:
  regions      from 07: whole_section, wm_lipid / gm_lipid (always), and the atlas regions
               (white_matter, corpus_callosum, cortex, hippocampus, thalamus) for sections
               whose registration passed QC
  territories  from 08 (differential.run_territories), each territory as a region

For each unit and lipid: mean normalized intensity per animal (edge pixels left out,
at least differential.min_pixels_per_unit pixels), log2, then Mutant minus Control within
each slide. See src/cln3msi/stats.py for the tests. Isotope peaks aren't tested (they just
mirror their monoisotopic peak).

FDR (Benjamini-Hochberg) families, within each mode and unit:
  q_targeted   only the peaks tagged with a hypothesis class (Cer, HexCer, SHexCer, ...)
  q_all        every tested peak (the untargeted view)

Position flag: Control is always on top. From 03 we know, per slide, how much brighter the
bottom of the slide was than the top in the off-tissue DAN signal. If a hit's log2FC has
the same sign and is less than qc.gradient.flag_ratio times that, it's flagged
"could_be_position". It's a warning, not a correction.

Then, from the same per-animal means (see src/cln3msi/classes.py):
  class tests        totals per lipid class and within-animal ratios (config class_tests),
                     one BH family per mode
  slide consistency  do the slides' KO/Control difference profiles agree? A slide that is
                     the mirror image of the others is flagged as a possible label swap
                     (config pair_check). A warning, never an automatic fix.

Writes results/09_differential/<mode>/differential.csv, per_animal_means.csv,
class_tests.csv, class_members.csv, pair_consistency.csv, pair_correlation_<unit>.csv,
volcanoes and class/consistency figures, summary (hypothesis classes first, then untargeted).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import classes, cli, stats  # noqa: E402
from cln3msi.util import read_h5ad  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def gradient_prediction(ctx, mode: str, obs: pd.DataFrame, mutant: str) -> tuple[float, pd.Series]:
    """Median over slides of the predicted log2(Mutant/Control) from the off-tissue trend."""
    p = ctx.output_root / "results" / "03_qc_all" / mode / "gradient.csv"
    if not p.exists():
        return np.nan, pd.Series(dtype=float)
    g = pd.read_csv(p).set_index("slide")
    col = next((c for c in ("matrix_off_pred_log2_bottom_vs_top", "tic_off_pred_log2_bottom_vs_top") if c in g), None)
    if col is None:
        return np.nan, pd.Series(dtype=float)
    mut_pos = obs[obs.genotype == mutant].drop_duplicates("slide").set_index("slide")["section"]
    sign = pd.Series({s: (1.0 if mut_pos.get(s) == "Bottom" else -1.0) for s in g.index})
    pred = g[col] * sign
    return float(pred.median()), pred


SLIDE_COLOURS = ["#c0392b", "#e67e22", "#d4ac0d", "#2874a6", "#7d3c98", "#17a589", "#7f8c8d", "#1abc9c"]


def class_level(ctx, mode, ML, animals, ctl, mut, S, log):
    """Class totals, ratios and the slide consistency check, from the per-animal means."""
    C = ctx.cfg.get("class_tests", {}) or {}
    PC = ctx.cfg.get("pair_check", {}) or {}
    out = lambda name: ctx.out("results", "09_differential", mode, name)  # noqa: E731
    min_pairs = int(ctx.cfg["differential"].get("min_pairs", 5))
    alpha = float(ctx.cfg["differential"]["fdr_alpha"])

    # ---- slide consistency first: the class summary refers to it --------------------------
    flagged = []
    if PC.get("enabled", True):
        units = [u for u in PC.get("units", ["whole_section"]) if u in set(ML["unit"])]
        mats, T, flagged = classes.consistency(ML, animals, ctl, mut, units, float(PC.get("flag_below", -0.3)),
                                                  int(PC.get("min_peaks", 20)))
        S.h("Slide consistency (do the KO/Control differences agree between slides?)")
        if T.empty:
            S.p(f"Not enough slides or peaks to check (needs 3 slides and pair_check.min_peaks = {int(PC.get('min_peaks', 20))} peaks).")
        else:
            T.to_csv(out("pair_consistency.csv"))
            for u, M in mats.items():
                M.to_csv(out(f"pair_correlation_{u}.csv"))
            S.p("Spearman correlation, over all tested peaks, between the slides' log2(Mutant/Control) "
                "profiles. With correct labels the slides share the genotype effect, so the correlations "
                "should be positive or near 0. The table gives each slide's median correlation with the "
                f"other slides; a slide is flagged when that is below {float(PC.get('flag_below', -0.3)):g} "
                "in at least half of the regions checked.")
            S.table(T.reset_index().round(3))
            if flagged:
                S.warn(f"slide(s) {', '.join(flagged)}: the KO/Control difference profile is the mirror image "
                       "of the other slides. This looks like the two sections are swapped (sample sheet, or "
                       "slide loaded the other way round). Check the sectioning/mounting records; only if they "
                       "confirm it, add the slide to qc.sections.flip_slides and rerun from 03. Do not relabel "
                       "because of this result alone.", log)
            else:
                S.p("No slide flagged.")
            fig, axes = plt.subplots(1, len(mats), figsize=(3.3 * len(mats), 3.0), squeeze=False)
            for ax, (u, M) in zip(axes[0], mats.items()):
                im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
                ax.set_xticks(range(len(M)), M.columns, fontsize=7)
                ax.set_yticks(range(len(M)), M.index, fontsize=7)
                for i in range(len(M)):
                    for j in range(len(M)):
                        ax.text(j, i, f"{M.iat[i, j]:.2f}", ha="center", va="center", fontsize=6)
                ax.set_title(u, fontsize=8)
            fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, label="Spearman r")
            fig.suptitle(f"{mode}: do the slides' KO/Control differences agree?", fontsize=9)
            S.fig(fig, "pair_consistency.png")

    # ---- class and ratio tests --------------------------------------------------------------
    if not C.get("enabled", True):
        return
    class_cfg = (C.get("classes") or {}).get(mode) or {}
    ratio_cfg = (C.get("ratios") or {}).get(mode) or []
    if not class_cfg:
        log.info(f"{mode}: no class_tests.classes for this mode, class tests skipped")
        return
    ap = ctx.out("processed", mode, "annotations.csv", mkdir=False)
    if not ap.exists():
        S.warn("no annotations.csv (run 05): class tests skipped", log)
        return
    ann = classes.usable_peaks(pd.read_csv(ap), int(C.get("min_sections", 0)), C.get("ambiguity", "tags"))
    regions = [u for u in C.get("regions", ["whole_section"]) if u in set(ML["unit"])]
    vals, mem = classes.totals(ML, ann, class_cfg, ratio_cfg, regions)
    if vals.empty:
        S.warn("none of the configured classes were found among the tested peaks", log)
        return
    R = classes.test_features(vals, animals, ctl, mut, min_pairs)
    if PC.get("sensitivity", True) and flagged:
        R = R.merge(classes.sensitivity(vals, animals, ctl, mut, flagged, min_pairs), on=["unit", "feature"], how="left")
    R.to_csv(out("class_tests.csv"), index=False)
    mem.to_csv(out("class_members.csv"), index=False)
    vals.to_csv(out("class_values_per_animal.csv"), index=False)

    S.h("Class-level and ratio tests (pre-specified, one BH family for this mode)")
    S.p(f"{len(R)} tests ({R['feature'].nunique()} classes/ratios x {R['unit'].nunique()} regions). "
        f"Species per class are in class_members.csv. log2fc is the mean over slides of "
        "log2(Mutant/Control); pair_<slide> columns show each slide on its own.")
    pair_cols = [c for c in R.columns if c.startswith("pair_")]
    show = ["unit", "feature", "kind", "n_species", "log2fc", "pct_change", "ci_low", "ci_high", "p_ttest",
            "p_wilcoxon", "p_sign", "q_bh", "n_pairs_up", "n_pairs", "p_sex_difference_exploratory"] + pair_cols
    S.table(R.sort_values("p_ttest")[[c for c in show if c in R]].round(4), max_rows=30)
    S.kv({f"q < {alpha}": int((R["q_bh"] < alpha).sum()), "lowest q": float(R["q_bh"].min())})
    if flagged:
        tag = "+".join(flagged)
        S.p(f"Sensitivity (exploratory, slide(s) {tag} flagged above): the CSV has extra columns with "
            f"{tag} left out (`without_{tag}_*`) and flipped (`flipped_{tag}_*`). Use them to see how much "
            "the flagged slide matters, not as the result: choosing the analysis after seeing the data "
            "is circular. Only the records can decide whether the slide is swapped.")
        fc = f"flipped_{tag}_p_ttest"
        if fc in R:
            S.table(R.sort_values(fc)[["unit", "feature", "log2fc", "p_ttest", f"flipped_{tag}_log2fc", fc,
                                       f"flipped_{tag}_n_pairs_up", f"flipped_{tag}_q_bh", f"without_{tag}_p_ttest"]]
                    .round(4), max_rows=15)

    # dot plot: each slide's log2FC for a few chosen features, per region
    feats = [f for f in ((C.get("plot_features") or {}).get(mode) or []) if f in set(R["feature"])]
    if feats and pair_cols:
        regs = [u for u in regions if u in set(R["unit"])]
        fig, axes = plt.subplots(len(regs), 1, figsize=(1.3 * len(feats) + 2.5, 2.1 * len(regs)), squeeze=False, sharex=True)
        slides = [c.replace("pair_", "") for c in pair_cols]
        for ax, u in zip(axes[:, 0], regs):
            sub = R[R["unit"] == u].set_index("feature")
            for i, f in enumerate(feats):
                if f not in sub.index:
                    continue
                for j, s in enumerate(slides):
                    v = sub.at[f, f"pair_{s}"]
                    mk = "D" if s in flagged else "o"
                    ax.scatter(i + (j - (len(slides) - 1) / 2) * 0.09, v, s=22, marker=mk,
                               color=SLIDE_COLOURS[j % len(SLIDE_COLOURS)], edgecolor="k" if s in flagged else "none",
                               label=s if (i == 0 and u == regs[0]) else None, zorder=3)
            ax.axhline(0, color="0.5", lw=0.7)
            ax.set_ylabel(u.replace("_", " "), fontsize=7)
        axes[-1, 0].set_xticks(range(len(feats)), feats, fontsize=7, rotation=20, ha="right")
        fig.suptitle(f"{mode}: log2(Mutant/Control) per slide pair" +
                     (f" (diamonds = flagged slide {', '.join(flagged)})" if flagged else ""), fontsize=8, y=0.995)
        fig.legend(*axes[0, 0].get_legend_handles_labels(), fontsize=7,
                   ncol=len(slides), loc="upper center", bbox_to_anchor=(0.5, 0.975), frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        S.fig(fig, "class_tests_per_slide.png")

    # forest plot for whole section + white matter
    F = R[R["unit"].isin(["whole_section", "white_matter"]) & R["p_ttest"].notna()].sort_values(["unit", "kind", "feature"])
    if len(F):
        fig, ax = plt.subplots(figsize=(6.5, 0.24 * len(F) + 1.2))
        yy = np.arange(len(F))
        ax.errorbar(F["log2fc"], yy, xerr=[F["log2fc"] - F["ci_low"], F["ci_high"] - F["log2fc"]],
                    fmt="o", ms=3, color="k", ecolor="0.6")
        ax.set_yticks(yy, F["feature"] + "  [" + F["unit"] + "]", fontsize=6)
        ax.axvline(0, color="k", lw=0.5)
        ax.invert_yaxis()
        ax.set_xlabel("log2(Mutant/Control), 95% CI over slide pairs")
        ax.set_title(f"{mode}: class totals and ratios", fontsize=9)
        S.fig(fig, "class_tests_forest.png")


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("09_differential", __doc__)
    D = ctx.cfg["differential"]
    G = ctx.cfg["project"]["genotypes"]
    ctl, mut = G["control"], G["mutant"]
    alpha = float(D["fdr_alpha"])
    flag_ratio = float(ctx.get("qc.gradient.flag_ratio", 2.0))
    for mode in modes:
        S = Summary(ctx, "09_differential", mode, f"Differential analysis, {mode}: {ctx.cfg['project']['mutant_label']} vs Control")
        p = ctx.out("processed", mode, "matrix_norm.h5ad", mkdir=False)
        if not p.exists():
            S.warn("no matrix_norm.h5ad, run 06 first", log)
            S.write()
            continue
        a = read_h5ad(p)
        obs, var = a.obs, a.var
        iso = var.get("is_isotope", pd.Series(False, index=var.index)).fillna(False).astype(bool).to_numpy()
        test_cols = np.flatnonzero(~iso)
        family = D.get("untargeted_peaks", "annotated_or_structured")
        if family != "all":
            annotated = var.get("annotation", pd.Series("", index=var.index)).fillna("").astype(str) != ""
            mp = ctx.output_root / "results" / "08_explore" / mode / "morans_i.csv"
            structured = pd.Series(False, index=var.index)
            if mp.exists():
                mor = pd.read_csv(mp)
                structured = var.index.isin(mor.loc[mor["selected"], "peak"])
            else:
                S.warn("no 08 Moran's I results: the untargeted family is the annotated peaks only", log)
            keep = (~iso) & (annotated.to_numpy() | np.asarray(structured))
            log.info(f"{mode}: untargeted family = {int(keep.sum())} peaks "
                     f"({int((annotated & ~pd.Series(iso, index=var.index)).sum())} annotated, "
                     f"{int(np.asarray(structured).sum())} spatially structured) of {int((~iso).sum())} non-isotope peaks")
            test_cols = np.flatnonzero(keep)
        X = np.asarray(a.X[:, test_cols], dtype=np.float64)
        peaks = var.index[test_cols]
        core = obs["core"].to_numpy().astype(bool)

        units = {}
        rp = ctx.out("processed", mode, "regions.parquet", mkdir=False)
        if rp.exists():
            reg = pd.read_parquet(rp).loc[obs.index]
            atlas_cols = [c for c in reg.columns if c.startswith("region_")
                          and c not in ("region_whole_section", "region_wm_lipid", "region_gm_lipid")]
            have_atlas = any(reg[c].any() for c in atlas_cols)
            for c in [c for c in reg.columns if c.startswith("region_")]:
                if c in ("region_wm_lipid", "region_gm_lipid") and have_atlas \
                        and D.get("lipid_wm_gm", "fallback_only") == "fallback_only":
                    continue
                units[("region", c.replace("region_", ""))] = reg[c].to_numpy().astype(bool)
        else:
            S.warn("no regions.parquet (run 07): only whole sections are tested", log)
            units[("region", "whole_section")] = np.ones(len(obs), bool)
        tp = ctx.out("processed", mode, "territories.parquet", mkdir=False)
        if D.get("run_territories", True) and tp.exists():
            ter = pd.read_parquet(tp).loc[obs.index, "territory"].to_numpy()
            for t in np.unique(ter):
                units[("territory", f"T{t}")] = ter == t

        animals = obs.drop_duplicates("animal_id").set_index("animal_id")[["slide", "genotype", "sex", "section"]]
        g_med, g_slide = gradient_prediction(ctx, mode, obs, mut)
        log.info(f"{mode}: {len(peaks)} peaks x {len(units)} units; gradient-predicted log2(Mut/Ctl) = {g_med:.3f}")

        results, means_long = [], []
        for (kind, name), m in units.items():
            mask = m & core
            if mask.sum() == 0:
                continue
            means, npx = stats.per_animal_means(X, obs, mask, int(D["min_pixels_per_unit"]))
            means.columns = peaks
            if means.notna().any(axis=1).sum() < 4:
                log.info(f"  {kind} {name}: fewer than 4 animals with enough pixels, skipped")
                continue
            c = stats.pseudocount(means)
            log2 = np.log2(means + c)
            res = stats.paired_tests(log2, animals, ctl, mut, run_mixed=bool(D.get("run_mixedlm", False)),
                                     min_pairs=int(D.get("min_pairs", 5)))
            res.insert(0, "unit", name)
            res.insert(0, "unit_kind", kind)
            res["median_pixels_per_animal"] = float(npx.median())
            results.append(res)
            ml = means.stack().rename("mean_norm_intensity").reset_index()
            ml.columns = ["animal_id", "peak", "mean_norm_intensity"]
            ml.insert(0, "unit", name)
            ml["n_pixels"] = ml["animal_id"].map(npx)
            means_long.append(ml)
        if not results:
            S.warn("nothing could be tested", log)
            S.write()
            continue
        R = pd.concat(results, ignore_index=True)
        ann_cols = [c for c in ("mz", "label", "annotation", "adduct", "lipid_class", "tag", "hypothesis",
                                "ambiguous", "typeII_flag", "candidates") if c in var]
        R = R.join(var[ann_cols], on="peak")
        R["targeted"] = R.get("hypothesis", pd.Series("", index=R.index)).fillna("").astype(str) != ""
        R["q_all"] = np.nan
        R["q_targeted"] = np.nan
        for _, idx in R.groupby("unit").groups.items():
            R.loc[idx, "q_all"] = stats.bh(R.loc[idx, "p_ttest"])
            t_idx = [i for i in idx if R.loc[i, "targeted"]]
            if t_idx:
                R.loc[t_idx, "q_targeted"] = stats.bh(R.loc[t_idx, "p_ttest"])
        exp_dir = R.get("hypothesis", pd.Series("", index=R.index)).fillna("")
        R["matches_hypothesis"] = np.where(exp_dir == "accumulated", R["log2fc"] > 0,
                                           np.where(exp_dir == "depleted", R["log2fc"] < 0,
                                                    np.where(exp_dir == "both", True, np.nan)))
        R["gradient_pred_log2fc"] = g_med
        if np.isfinite(g_med) and g_med != 0:
            same = np.sign(R["log2fc"]) == np.sign(g_med)
            R["could_be_position"] = same & (R["log2fc"].abs() < flag_ratio * abs(g_med))
        else:
            R["could_be_position"] = np.nan
        R["significant"] = (np.where(R["targeted"], R["q_targeted"], R["q_all"]) < alpha)
        R = R.sort_values(["unit_kind", "unit", "p_ttest"])
        R.to_csv(ctx.out("results", "09_differential", mode, "differential.csv"), index=False)
        ML = pd.concat(means_long, ignore_index=True)
        ML.to_csv(ctx.out("results", "09_differential", mode, "per_animal_means.csv"), index=False)

        # ---- volcanoes ------------------------------------------------------------------
        reg_units = [u for k, u in units if k == "region" and u in set(R.unit)]
        ncol = min(4, len(reg_units))
        nrow = int(np.ceil(len(reg_units) / max(ncol, 1)))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.2 * nrow), squeeze=False)
        for ax, u in zip(axes.ravel(), reg_units):
            sub = R[R.unit == u]
            col = np.where(sub["significant"], np.where(sub["log2fc"] > 0, "tab:red", "tab:blue"), "0.7")
            ax.scatter(sub["log2fc"], -np.log10(sub["p_ttest"].clip(lower=1e-300)), c=col, s=12,
                       edgecolor=np.where(sub["targeted"], "k", "none"), linewidth=0.5)
            if np.isfinite(g_med):
                ax.axvspan(min(0, flag_ratio * g_med), max(0, flag_ratio * g_med), color="orange", alpha=0.15)
            for _, r in sub[sub["significant"]].nsmallest(6, "p_ttest").iterrows():
                ax.annotate(str(r.get("label", r["peak"]))[:22], (r["log2fc"], -np.log10(r["p_ttest"])), fontsize=5)
            ax.axvline(0, color="k", lw=0.5)
            ax.set_title(u, fontsize=9)
            ax.set_xlabel(f"log2FC ({ctx.cfg['project']['mutant_label']} / Control)", fontsize=7)
            ax.set_ylabel("-log10 p (paired t)", fontsize=7)
        for ax in axes.ravel()[len(reg_units):]:
            ax.axis("off")
        fig.suptitle(f"{mode}: volcano per region (black edge = hypothesis class, orange band = within "
                     f"{flag_ratio:g}x the position gradient)", fontsize=8)
        S.fig(fig, "volcano_regions.png")

        # hypothesis classes: forest plot in whole section + white matter
        hyp = R[R["targeted"] & R.unit.isin(["whole_section", "white_matter", "corpus_callosum"])
                & R["p_ttest"].notna()].sort_values(["unit", "lipid_class", "mz"] if "lipid_class" in R else ["unit"])
        if len(hyp):
            fig, ax = plt.subplots(figsize=(7, 0.3 * len(hyp) + 1.5))
            yl = hyp["label"].astype(str) + "  [" + hyp["unit"] + "]"
            yy = np.arange(len(hyp))
            ax.errorbar(hyp["log2fc"], yy, xerr=[hyp["log2fc"] - hyp["ci_low"], hyp["ci_high"] - hyp["log2fc"]],
                        fmt="o", ms=4, color="k", ecolor="0.5")
            ax.set_yticks(yy, yl, fontsize=6)
            ax.axvline(0, color="k", lw=0.5)
            if np.isfinite(g_med):
                ax.axvspan(min(0, flag_ratio * g_med), max(0, flag_ratio * g_med), color="orange", alpha=0.15)
            ax.set_xlabel("log2FC with 95% CI (paired over slides)")
            ax.set_title(f"{mode}: hypothesis classes", fontsize=9)
            S.fig(fig, "forest_hypothesis.png")

        # ---- summary --------------------------------------------------------------------
        S.h("Design")
        pairs = R["pairs"].iloc[0]
        S.p(f"Pairs (slides): {pairs}. Unit = animal; mean of TIC-normalized intensity over the region's "
            "non-edge pixels, log2, Mutant minus Control within slide.")
        S.p(f"Position gradient: median predicted log2(Mutant/Control) from the off-tissue DAN trend = "
            f"{g_med:.3f}" + ("" if np.isfinite(g_med) else " (not available)") +
            f". Hits with the same sign and |log2FC| < {flag_ratio:g}x that are flagged `could_be_position`.")
        show = ["unit", "label", "hypothesis", "log2fc", "ci_low", "ci_high", "p_ttest", "q_targeted",
                "p_wilcoxon", "p_mixedlm", "n_pairs", "n_pairs_up", "matches_hypothesis", "could_be_position",
                "ambiguous"]
        show = [c for c in show if c in R]
        S.h("1. Hypothesis classes (targeted, BH within region over tagged peaks)")
        T = R[R["targeted"]].sort_values(["p_ttest"])
        S.table(T[show], max_rows=40)
        S.h("2. Untargeted: top hits per region (BH over all peaks in the region)")
        U = R[R.unit_kind == "region"].sort_values("p_ttest").groupby("unit").head(5)
        S.table(U[[c if c != "q_targeted" else "q_all" for c in show]], max_rows=60)
        if (R.unit_kind == "territory").any():
            S.h("3. Territories: top hits")
            Tt = R[R.unit_kind == "territory"].sort_values("p_ttest").head(15)
            S.table(Tt[[c if c != "q_targeted" else "q_all" for c in show]])
        S.h("Counts")
        sig = R[R["significant"]]
        S.kv({"tests": len(R), f"significant (q < {alpha})": len(sig),
              "of which flagged could_be_position": int(sig["could_be_position"].fillna(False).astype(bool).sum()),
              "units": R.unit.nunique()})
        if "p_sex_difference_exploratory" in R:
            S.h("Sex (exploratory only)")
            S.p("Females were run in runs 1-3, all males in run 4, so sex and run are completely confounded. "
                "`p_sex_difference_exploratory` compares the per-slide log2FC of 3 female vs 3 male pairs; "
                "treat it as a hint for follow-up, not a result.")
            Sx = R[R.unit == "whole_section"].nsmallest(10, "p_sex_difference_exploratory")
            S.table(Sx[[c for c in ("label", "log2fc", "log2fc_Female", "log2fc_Male",
                                     "p_sex_difference_exploratory") if c in Sx]].round(4))
        class_level(ctx, mode, ML, animals, ctl, mut, S, log)
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
