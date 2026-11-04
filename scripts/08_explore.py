"""Phase 08: spatially structured lipids, lipid programs and lipid territories.

    python scripts/08_explore.py --config config.yaml

Same recipe as the tutorial notebooks 05-06, adapted to 12 sections:
  1. Moran's I of every lipid (isotope peaks left out), separately in each section, on a
     k-nearest-neighbour graph of pixel positions, with a permutation p-value.
     Keep lipids with I >= min_I in >= min_section_fraction of the sections.
  2. Seeded NMF on Control pixels only (values = normalized intensity, rescaled 0..1 per
     lipid with percentiles taken on Control), k picked at the elbow of the error curve
     unless explore.nmf.k is a number. Then applied to every pixel with the programs fixed.
  3. Leiden clustering of the NMF activities of a Control subsample = lipid territories,
     transferred to all pixels (rest of Control and all Mutant) by a kNN vote.

Nothing here is a statistical test of genotype. 09 uses the territories as extra regions.

Writes processed/<mode>/territories.parquet, results/08_explore/<mode>/{morans_i.csv,
nmf_programs.csv, territory_composition.csv}, figures, summary.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, spatial  # noqa: E402
from cln3msi.util import read_h5ad  # noqa: E402
from cln3msi.images import section_image  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("08_explore", __doc__)
    E = ctx.cfg["explore"]
    if E.get("mode"):
        modes = [E["mode"]] if isinstance(E["mode"], str) else list(E["mode"])
    seed = int(ctx.get("runtime.seed", 0))
    rng = np.random.default_rng(seed)
    ctl = ctx.cfg["project"]["genotypes"]["control"]
    for mode in modes:
        S = Summary(ctx, "08_explore", mode, f"Spatial exploration, {mode}")
        p = ctx.out("processed", mode, "matrix_norm.h5ad", mkdir=False)
        if not p.exists():
            S.warn("no matrix_norm.h5ad, run 06 first", log)
            S.write()
            continue
        a = read_h5ad(p)
        obs = a.obs
        iso = a.var.get("is_isotope", pd.Series(False, index=a.var_names)).fillna(False).astype(bool).to_numpy()
        feats = np.flatnonzero(~iso)
        X = np.asarray(a.X[:, feats], dtype=np.float64)
        labels = a.var["label"].iloc[feats].astype(str).to_numpy() if "label" in a.var else a.var_names[feats]
        log.info(f"{mode}: {len(feats)} non-isotope peaks, {a.n_obs} pixels")

        # ---- 1. Moran's I per section --------------------------------------------------
        M = E["morans_i"]
        recs = []
        keys = obs.groupby(["slide", "section"]).groups
        for (slide, section), idx in keys.items():
            rows = obs.index.get_indexer(idx)
            core = obs["core"].to_numpy()[rows].astype(bool)
            rows = rows[core] if core.sum() > 30 else rows
            xy = obs[["x", "y"]].to_numpy()[rows].astype(float)
            W = spatial.knn_weights(xy, int(M["k_neighbors"]))
            I, pv = spatial.morans_i(X[rows], W, int(M["n_perm"]), rng)
            for j in range(len(feats)):
                recs.append({"slide": slide, "section": section, "peak": a.var_names[feats[j]],
                             "label": labels[j], "I": I[j], "p_perm": pv[j]})
        mor = pd.DataFrame(recs)
        mor["q"] = multipletests(mor["p_perm"].fillna(1), method="fdr_bh")[1]
        mor.to_csv(ctx.out("results", "08_explore", mode, "morans_i_per_section.csv"), index=False)
        agg = mor.groupby(["peak", "label"]).agg(I_median=("I", "median"), I_min=("I", "min"),
                                                 frac_above=("I", lambda v: (v >= float(M["min_I"])).mean()),
                                                 frac_q05=("q", lambda v: (v < 0.05).mean())).reset_index()
        agg["selected"] = agg["frac_above"] >= float(M["min_section_fraction"])
        agg = agg.sort_values("I_median", ascending=False)
        agg.to_csv(ctx.out("results", "08_explore", mode, "morans_i.csv"), index=False)
        sel_peaks = agg.loc[agg.selected, "peak"].tolist()
        if len(sel_peaks) < 3:
            S.warn(f"only {len(sel_peaks)} lipids pass the Moran filter; using the top 10 by median I instead", log)
            sel_peaks = agg["peak"].head(min(10, len(agg))).tolist()
        cols = [list(a.var_names[feats]).index(pk) for pk in sel_peaks]
        V_all = X[:, cols].copy()
        if E.get("per_slide_scaling", True):
            # each lipid divided by its median on that slide (both sections pooled, so the
            # Control vs KO difference within a slide is kept), to stop territories from
            # simply sorting pixels by slide
            for s_ in obs["slide"].unique():
                rows_s = (obs["slide"] == s_).to_numpy()
                med = np.median(V_all[rows_s], axis=0)
                V_all[rows_s] /= np.where(med > 0, med, 1)
        is_ctl = (obs["genotype"] == ctl).to_numpy()
        lo = np.percentile(V_all[is_ctl], 0.5, axis=0)
        hi = np.percentile(V_all[is_ctl], 99.5, axis=0)
        V_all = np.clip((V_all - lo) / np.where(hi > lo, hi - lo, 1), 0, 1)

        # ---- 2. seeded NMF on Control --------------------------------------------------
        Nc = E["nmf"]
        ctl_idx = np.flatnonzero(is_ctl)
        tr = rng.choice(ctl_idx, min(len(ctl_idx), int(Nc["pixel_subsample"])), replace=False)
        grid = [k for k in Nc["k_grid"] if k < len(cols)] or [max(2, len(cols) - 1)]
        errs = {}
        for k in grid:
            m, _, _ = spatial.seeded_nmf(V_all[tr], k, int(Nc["max_iter"]), seed)
            errs[k] = float(m.reconstruction_err_)
        k_use = spatial.choose_k(errs) if Nc.get("k", "auto") == "auto" else int(Nc["k"])
        k_use = min(k_use, len(cols) - 1) if len(cols) > 2 else 1
        model, _, seeds = spatial.seeded_nmf(V_all[tr], k_use, int(Nc["max_iter"]), seed)
        W_all = model.transform(V_all)
        H = pd.DataFrame(model.components_, columns=[labels[list(a.var_names[feats]).index(pk)] for pk in sel_peaks],
                         index=[f"program_{i + 1}" for i in range(k_use)])
        H.to_csv(ctx.out("results", "08_explore", mode, "nmf_programs.csv"))
        order = np.argsort(-W_all.var(axis=0))

        # ---- 3. territories ------------------------------------------------------------
        T = E["territories"]
        trt = rng.choice(ctl_idx, min(len(ctl_idx), int(T["pixel_subsample"])), replace=False)
        lab_tr = spatial.territories(W_all[trt], int(T["k_neighbors"]), float(T["resolution"]), seed, log)
        # merge clusters with < 0.5% of the training pixels: drop them before transfer
        vals, cnt = np.unique(lab_tr, return_counts=True)
        big = vals[cnt >= max(5, 0.005 * len(lab_tr))]
        keep = np.isin(lab_tr, big)
        terr = spatial.knn_transfer(W_all[trt][keep], lab_tr[keep], W_all, k=15)
        # renumber 1..n by size
        vals, cnt = np.unique(terr, return_counts=True)
        remap = {v: i + 1 for i, v in enumerate(vals[np.argsort(-cnt)])}
        terr = np.array([remap[v] for v in terr])
        out = pd.DataFrame(W_all, index=obs.index, columns=[f"nmf_{i + 1}" for i in range(k_use)])
        out["territory"] = terr
        out["trained_on"] = False
        out.iloc[trt, out.columns.get_loc("trained_on")] = True
        out.to_parquet(ctx.out("processed", mode, "territories.parquet"))

        comp = pd.crosstab(obs["animal_id"], terr, normalize="index")
        comp.columns = [f"T{c}" for c in comp.columns]
        meta = obs.drop_duplicates("animal_id").set_index("animal_id")[["slide", "genotype", "sex"]]
        comp = meta.join(comp)
        comp.to_csv(ctx.out("results", "08_explore", mode, "territory_composition.csv"))

        # ---- figures ---------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(7, max(3, 0.18 * min(40, len(agg)))))
        top = agg.head(40)
        ax.barh(top["label"][::-1], top["I_median"][::-1],
                color=["tab:green" if s else "0.6" for s in top["selected"][::-1]])
        ax.axvline(float(M["min_I"]), color="k", lw=0.5, ls="--")
        ax.set_xlabel("median Moran's I across sections (green = selected)")
        ax.tick_params(axis="y", labelsize=6)
        S.fig(fig, "morans_i.png")

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot(list(errs), list(errs.values()), marker="o")
        ax.axvline(k_use, color="tab:red", ls="--")
        ax.set_xlabel("k (programs)"); ax.set_ylabel("NMF reconstruction error")  # noqa: E702
        S.fig(fig, "nmf_k.png")

        fig, ax = plt.subplots(figsize=(min(16, 0.35 * H.shape[1] + 3), 0.4 * k_use + 2))
        im = ax.imshow(H.to_numpy() / (H.to_numpy().max(axis=1, keepdims=True) + 1e-12), aspect="auto", cmap="viridis")
        ax.set_yticks(range(k_use), H.index, fontsize=7)
        ax.set_xticks(range(H.shape[1]), H.columns, rotation=90, fontsize=6)
        plt.colorbar(im, ax=ax, fraction=0.03)
        ax.set_title("NMF programs (rows scaled to max 1)", fontsize=9)
        S.fig(fig, "nmf_programs.png")

        secs = list(keys)
        show_p = order[:min(4, k_use)]
        fig, axes = plt.subplots(len(show_p) + 1, len(secs), figsize=(1.5 * len(secs), 1.4 * (len(show_p) + 1)),
                                 squeeze=False)
        n_terr = int(terr.max())
        cmap_t = plt.get_cmap("tab20", max(n_terr, 2))
        for b, (slide, section) in enumerate(secs):
            for r_, pi in enumerate(show_p):
                img, *_ = section_image(obs, W_all[:, pi], slide, section)
                axes[r_, b].imshow(img, cmap="magma", vmin=0, vmax=np.percentile(W_all[:, pi], 99), interpolation="nearest")
                if b == 0:
                    axes[r_, b].set_ylabel(f"prog {pi + 1}", fontsize=7)
            img, *_ = section_image(obs, terr.astype(float), slide, section)
            axes[-1, b].imshow(img, cmap=cmap_t, vmin=0.5, vmax=n_terr + 0.5, interpolation="nearest")
            g = obs.loc[(obs.slide == slide) & (obs.section == section), "genotype"].iloc[0]
            axes[0, b].set_title(f"{slide} {section[0]}\n{'KO' if g != ctl else 'Ctl'}", fontsize=7)
            if b == 0:
                axes[-1, b].set_ylabel("territory", fontsize=7)
        for ax in axes.ravel():
            ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
        fig.suptitle(f"{mode}: top NMF programs and territories (learned on Control, transferred)", fontsize=9)
        S.fig(fig, "programs_territories.png")

        S.h("Moran's I")
        S.kv({"lipids tested (no isotopes)": len(feats), "selected (spatially structured)": int(agg.selected.sum()),
              "threshold": f"I >= {M['min_I']} in >= {M['min_section_fraction']:.0%} of sections"})
        S.table(agg.head(20).round(3))
        S.h("NMF")
        S.kv({"k grid errors": ", ".join(f"{k}: {v:.3g}" for k, v in errs.items()), "k used": k_use,
              "seed lipids": ", ".join(H.columns[i] for i in seeds)})
        for i in order[:min(6, k_use)]:
            top_l = H.iloc[i].sort_values(ascending=False).head(5)
            S.p(f"- program {i + 1}: " + ", ".join(top_l.index))
        S.h("Territories")
        S.p(f"{n_terr} territories. Fraction of each animal's pixels per territory (descriptive):")
        S.table(comp.reset_index().round(3))
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
