"""Phase 11 (optional): compare this run with another run of the pipeline, mode by mode.

    python scripts/11_compare_runs.py --config config.flipJ.yaml --baseline config.yaml
    python scripts/11_compare_runs.py --config config.flipJ.yaml --baseline /Users/neharao/cln3_msi_outputs

Made for the slide J sensitivity run (config.flipJ.yaml: same data, slide J's sections
swapped), but it works for any two runs, e.g. before/after a config change. --baseline is
either a config file (its output_root is used) or an output folder. The baseline outputs
are only read; everything is written into THIS run's output folder:

  results/11_compare/<mode>/class_tests_compared.csv   every class/ratio test, both runs side by side
  results/11_compare/<mode>/differential_compared.csv  every single-peak test, both runs
  results/11_compare/<mode>/slide_consistency.csv      each slide's agreement with the others, both runs
  results/11_compare/<mode>/counts.csv                 significant results per family, both runs
  figures/11_compare/<mode>/...                        scatter plots of the fold changes
  summaries/11_compare_<mode>.md

Nothing here is a new test: it lines up the tests the two runs already did. Keep in mind
that choosing between the runs because one gives nicer results would be circular; the
labels have to come from the records.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, provenance  # noqa: E402
from cln3msi.config import ConfigError, load  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def add_args(p):
    p.add_argument("--baseline", required=True,
                   help="the run to compare with: a config file (its output_root) or an output folder")
    p.add_argument("--labels", default="baseline,this run",
                   help="names for the two runs in tables and figures, comma separated")


def baseline_root(arg: str) -> Path:
    p = Path(arg).expanduser()
    if p.suffix in (".yaml", ".yml") and p.is_file():
        return load(p).output_root
    if p.is_dir():
        return p.resolve()
    raise ConfigError(f"--baseline {arg}: not a config file or an existing output folder")


def read(root: Path, *parts) -> pd.DataFrame | None:
    f = root.joinpath(*parts)
    return pd.read_csv(f) if f.exists() else None


def side_by_side(A, B, keys, cols, la, lb):
    cols = [c for c in cols if c in A.columns or c in B.columns]
    a = A[keys + [c for c in cols if c in A]].rename(columns={c: f"{c}__{la}" for c in cols})
    b = B[keys + [c for c in cols if c in B]].rename(columns={c: f"{c}__{lb}" for c in cols})
    return a.merge(b, on=keys, how="outer")


def scatter(S, x, y, lab, la, lb, title, name, color_by=None):
    ok = np.isfinite(x) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    if color_by is None:
        ax.scatter(x[ok], y[ok], s=10, color="0.3")
    else:
        for g in pd.unique(color_by[ok]):
            sel = ok & (color_by == g)
            ax.scatter(x[sel], y[sel], s=12, label=str(g))
        ax.legend(fontsize=6, frameon=False, loc="upper left")
    lim = np.nanmax(np.abs(np.r_[x[ok], y[ok], 0.1])) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.7)
    ax.axhline(0, color="0.8", lw=0.5); ax.axvline(0, color="0.8", lw=0.5)  # noqa: E702
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)  # noqa: E702
    ax.set_xlabel(f"log2 fold change, {la}"); ax.set_ylabel(f"log2 fold change, {lb}")  # noqa: E702
    r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else np.nan
    ax.set_title(f"{title}\nPearson r = {r:.2f}, n = {int(ok.sum())}", fontsize=8)
    S.fig(fig, name)
    return r


def main():
    t0 = time.time()
    ctx, log, modes, args = cli.start("11_compare_runs", __doc__, extra_args=add_args)
    base = baseline_root(args.baseline)
    la, lb = [s.strip() for s in args.labels.split(",")][:2]
    if base.resolve() == ctx.output_root.resolve():
        raise SystemExit("--baseline is this run's own output folder; nothing to compare")
    alpha = float(ctx.cfg["differential"]["fdr_alpha"])
    log.info(f"comparing {la} = {base}  with  {lb} = {ctx.output_root}")
    have_base_prov = (base / "metadata" / "provenance").exists()
    this_fp = provenance.code_fingerprint()

    for mode in modes:
        S = Summary(ctx, "11_compare", mode, f"Run comparison, {mode}: {la} vs {lb}")
        out = lambda n: ctx.out("results", "11_compare", mode, n)  # noqa: E731
        S.h("Runs")
        flip_now = ctx.get("qc.sections.flip_slides") or []
        S.kv({la: str(base), lb: str(ctx.output_root), f"flip_slides in {lb}": ", ".join(flip_now) or "none"})
        if have_base_prov:
            import json
            def fp09(root):
                return {json.loads(p.read_text()).get("code_fingerprint")
                        for p in (root / "metadata" / "provenance").glob(f"09_differential_{mode}.json")}
            fa, fb = fp09(base), fp09(ctx.output_root) or {this_fp}
            same = bool(fa) and fa == fb
            S.p(f"Code used for 09: {la} {', '.join(sorted(map(str, fa))) or 'unknown'}, "
                f"{lb} {', '.join(sorted(map(str, fb)))}"
                + (" (same code)." if same else ". The code differs, so differences between the runs can come "
                   "from the code as well as from the settings; rerun the baseline from 09 with this code to isolate "
                   "the effect of the settings."))
        else:
            S.p(f"The {la} run has no provenance records (it ran before they were added), so its code version is "
                "not known. Rerun it from 09 with the current code if the comparison should isolate the settings.")

        # ---- slide consistency ----------------------------------------------------------
        A = read(base, "results", "09_differential", mode, "pair_consistency.csv")
        B = read(ctx.output_root, "results", "09_differential", mode, "pair_consistency.csv")
        S.h("Slide consistency (median correlation of each slide's KO/Control profile with the others)")
        if A is not None and B is not None:
            T = side_by_side(A, B, ["slide"], [c for c in B.columns if c != "slide"], la, lb)
            T.to_csv(out("slide_consistency.csv"), index=False)
            S.table(T.round(3))
        else:
            S.warn("pair_consistency.csv missing in one of the runs (09 older than the slide check?)", log)

        # ---- class tests ----------------------------------------------------------------
        A = read(base, "results", "09_differential", mode, "class_tests.csv")
        B = read(ctx.output_root, "results", "09_differential", mode, "class_tests.csv")
        counts = []
        S.h("Class-level and ratio tests")
        if A is not None and B is not None:
            cols = ["kind", "n_species", "primary", "log2fc", "ci_low", "ci_high", "p_ttest", "q_bh", "q_primary",
                    "p_gradient_corrected", "n_pairs_up", "n_pairs"]
            T = side_by_side(A, B, ["unit", "feature"], cols, la, lb)
            T["delta_log2fc"] = T[f"log2fc__{lb}"] - T[f"log2fc__{la}"]
            T = T.sort_values(f"p_ttest__{lb}")
            T.to_csv(out("class_tests_compared.csv"), index=False)
            for lab, D in ((la, A), (lb, B)):
                if "q_primary" in D and D["q_primary"].notna().any():
                    counts.append({"run": lab, "family": "primary endpoints", "tests": int(D["q_primary"].notna().sum()),
                                   f"q < {alpha}": int((D["q_primary"] < alpha).sum()), "lowest q": float(D["q_primary"].min())})
                counts.append({"run": lab, "family": "class and ratio tests", "tests": int(D["p_ttest"].notna().sum()),
                               f"q < {alpha}": int((D["q_bh"] < alpha).sum()), "lowest q": float(D["q_bh"].min())})
            show = ["unit", "feature", f"log2fc__{la}", f"p_ttest__{la}", f"q_bh__{la}", f"n_pairs_up__{la}",
                    f"log2fc__{lb}", f"p_ttest__{lb}", f"q_bh__{lb}", f"n_pairs_up__{lb}"]
            S.p(f"Sorted by the p-value in {lb}. n_pairs_up = pairs with KO above Control.")
            S.table(T[[c for c in show if c in T]].round(4), max_rows=25)
            scatter(S, T[f"log2fc__{la}"].to_numpy(float), T[f"log2fc__{lb}"].to_numpy(float), None, la, lb,
                    f"{mode}: class and ratio tests", "class_tests_scatter.png", T["unit"].to_numpy())
        else:
            S.warn("class_tests.csv missing in one of the runs", log)

        # ---- single peaks ---------------------------------------------------------------
        A = read(base, "results", "09_differential", mode, "differential.csv")
        B = read(ctx.output_root, "results", "09_differential", mode, "differential.csv")
        S.h("Single-peak tests")
        if A is not None and B is not None:
            keys = ["unit_kind", "unit", "peak"]
            cols = ["label", "targeted", "log2fc", "p_ttest", "p_primary", "q_targeted", "q_all", "significant",
                    "n_pairs_up", "could_be_position"]
            T = side_by_side(A, B, keys, cols, la, lb)
            pk = f"p_primary__{lb}" if f"p_primary__{lb}" in T else f"p_ttest__{lb}"
            T = T.sort_values(pk)
            T.to_csv(out("differential_compared.csv"), index=False)
            for lab, D in ((la, A), (lb, B)):
                sig = D["significant"].fillna(False).astype(bool)
                tg = D["targeted"].fillna(False).astype(bool)
                counts.append({"run": lab, "family": "single peaks, targeted", "tests": int((tg & D["p_ttest"].notna()).sum()),
                               f"q < {alpha}": int((sig & tg).sum()), "lowest q": float(D.loc[tg, "q_targeted"].min())})
                counts.append({"run": lab, "family": "single peaks, all", "tests": int(D["p_ttest"].notna().sum()),
                               f"q < {alpha}": int(sig.sum()), "lowest q": float(D["q_all"].min())})
            reg = T[T["unit_kind"] == "region"]
            r = scatter(S, reg[f"log2fc__{la}"].to_numpy(float), reg[f"log2fc__{lb}"].to_numpy(float), None, la, lb,
                        f"{mode}: single peaks, atlas regions", "differential_scatter.png", reg["unit"].to_numpy())
            S.p(f"Correlation of single-peak fold changes between the runs (atlas regions): r = {r:.2f}.")
            either = T[T[f"significant__{la}"].fillna(False).astype(bool) | T[f"significant__{lb}"].fillna(False).astype(bool)]
            S.p(f"Significant in either run: {len(either)}.")
            S.table(either[[c for c in ["unit", f"label__{lb}", f"log2fc__{la}", f"q_all__{la}", f"log2fc__{lb}",
                                        f"q_all__{lb}", f"n_pairs_up__{lb}"] if c in either]].round(4), max_rows=30)
        else:
            S.warn("differential.csv missing in one of the runs", log)

        if counts:
            C = pd.DataFrame(counts)
            C.to_csv(out("counts.csv"), index=False)
            S.h("Counts")
            S.table(C)

        # ---- registration (should be unchanged unless the sections were reassigned) ------
        A = read(base, "results", "07_register", mode, "registration.csv")
        B = read(ctx.output_root, "results", "07_register", mode, "registration.csv")
        if A is not None and B is not None and {"slide", "section", "dice"} <= set(A.columns) & set(B.columns):
            R = A[["slide", "section", "animal_id", "dice"]].merge(
                B[["slide", "section", "animal_id", "dice"]], on=["slide", "section"], suffixes=(f"__{la}", f"__{lb}"))
            S.h("Registration per section")
            S.p("Swapping a slide's labels changes which animal each section belongs to, not the images, so the "
                "Dice values should be the same and only animal_id should differ on the swapped slide.")
            S.table(R.round(3))
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
