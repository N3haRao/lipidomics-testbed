"""Phase 06: per-pixel normalization (TIC, RMS, median), compared by how much slide effect is left.

    python scripts/06_normalize.py --config config.yaml

Each method divides every pixel by a size factor:
  tic     sum of all kept (non-matrix) peaks in that pixel, or the full-spectrum TIC
          (normalization.tic_source: peaks | full_spectrum)
  rms     root mean square of the pixel's peak intensities
  median  median of the pixel's non-zero peak intensities
The factors are rescaled to a median of 1, so normalized values stay in "raw-ish" units.

How we compare them: for every peak we average each SECTION (12 values), take log2, and ask
how much of the variation between sections is explained by slide, after allowing for
genotype (partial R2 of slide). Lower = less slide effect. Plus PCA of the section profiles.

The method used for statistics is fixed in advance (normalization.for_stats, default tic),
so we don't pick whichever gives the nicest results later.
uMAIA (cross-slide) runs only if normalization.umaia.enabled and uMAIA imports; it's added
as a separate layer for comparison. A batch-corrected embedding (Harmony) is only ever
stored in obsm for plotting, and nothing in 09 reads it.

Writes processed/<mode>/matrix_norm.h5ad (X = the for_stats normalization; obs has all size
factors so any other method can be recomputed from matrix_raw.h5ad), figures, summary.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def size_factors(X: np.ndarray, obs: pd.DataFrame, tic_source: str) -> dict:
    f = {}
    f["tic"] = obs["tic_spectrum"].to_numpy(float) if tic_source == "full_spectrum" else X.sum(axis=1, dtype=np.float64)
    f["rms"] = np.sqrt(np.mean(X.astype(np.float64) ** 2, axis=1))
    Xn = np.where(X > 0, X, np.nan)
    with np.errstate(all="ignore"):
        f["median"] = np.nanmedian(Xn, axis=1)
    for k, v in f.items():
        v = np.where(np.isfinite(v) & (v > 0), v, np.nan)
        med = np.nanmedian(v)
        f[k] = np.where(np.isfinite(v), v / med, 1.0)   # a pixel with no signal keeps factor 1
    return f


def section_means(X: np.ndarray, obs: pd.DataFrame) -> pd.DataFrame:
    key = obs["slide"].astype(str) + "|" + obs["section"].astype(str)
    df = pd.DataFrame(X, index=key.to_numpy())
    return df.groupby(level=0).mean()


def slide_partial_r2(means: pd.DataFrame, sec_info: pd.DataFrame) -> np.ndarray:
    """Per peak: partial R2 of slide in log2(mean) ~ slide + genotype across the 12 sections."""
    y = np.log2(means.to_numpy() + 1e-9)
    info = sec_info.loc[means.index]
    G = pd.get_dummies(info["genotype"], drop_first=True).to_numpy(float)
    Sd = pd.get_dummies(info["slide"], drop_first=True).to_numpy(float)
    one = np.ones((len(info), 1))

    def sse(D):
        beta, *_ = np.linalg.lstsq(D, y, rcond=None)
        r = y - D @ beta
        return (r ** 2).sum(axis=0)

    sse_g = sse(np.hstack([one, G]))
    sse_full = sse(np.hstack([one, G, Sd]))
    with np.errstate(all="ignore"):
        return np.where(sse_g > 0, 1 - sse_full / sse_g, np.nan)


def pca2(M: np.ndarray) -> np.ndarray:
    Z = np.log2(M + 1e-9)
    Z = Z - Z.mean(axis=0)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    return U[:, :2] * s[:2]


def run_umaia(X, obs, cfg, log):
    """Optional. Shapes follow the tutorial: (N_pixels, S_sections, V) + boolean mask, log space."""
    try:
        import uMAIA  # noqa: F401
        from uMAIA import norm as unorm
    except Exception as e:
        log.info(f"uMAIA not available ({e.__class__.__name__}); skipped")
        return None
    try:
        key = (obs["slide"].astype(str) + "|" + obs["section"].astype(str)).to_numpy()
        secs = sorted(set(key))
        n_max = max((key == s).sum() for s in secs)
        eps = 2e-4
        x = np.zeros((n_max, len(secs), X.shape[1]), np.float32)
        mask = np.zeros((n_max, len(secs)), bool)
        for j, s in enumerate(secs):
            idx = np.flatnonzero(key == s)
            x[:len(idx), j] = np.log(X[idx] + eps)
            mask[:len(idx), j] = True
        init = unorm.initialize(x, mask, subsample=True)
        params = unorm.normalize(x, mask, *init if isinstance(init, tuple) else init,
                                 num_steps=int(cfg.get("num_steps", 2000)), seed=42, covariate_vector=None)
        xt = unorm.transform(x, mask, params)
        out = np.zeros_like(X, dtype=np.float32)
        for j, s in enumerate(secs):
            idx = np.flatnonzero(key == s)
            out[idx] = np.exp(np.asarray(xt)[:len(idx), j]) - eps
        return np.clip(out, 0, None)
    except Exception as e:
        log.warning(f"uMAIA failed ({e}); continuing without it. Its API may differ from the tutorial "
                    "version this was written for; send me the error.")
        return None


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("06_normalize", __doc__)
    N = ctx.cfg["normalization"]
    for mode in modes:
        S = Summary(ctx, "06_normalize", mode, f"Normalization, {mode}")
        p = ctx.out("processed", mode, "matrix_raw.h5ad", mkdir=False)
        if not p.exists():
            S.warn("no matrix_raw.h5ad, run 04 first", log)
            S.write()
            continue
        a = ad.read_h5ad(p)
        annp = ctx.out("processed", mode, "annotations.csv", mkdir=False)
        if annp.exists():
            ann = pd.read_csv(annp, index_col=0)
            a.var = a.var.join(ann.drop(columns=[c for c in ann.columns if c in a.var.columns]), how="left")
        else:
            S.warn("no annotations.csv (run 05); peaks will only have m/z labels", log)
            a.var["label"] = [f"mz {m:.4f}" for m in a.var["mz"]]
        X = np.asarray(a.X, dtype=np.float32)
        log.mem(f"loaded {X.shape}")

        f = size_factors(X, a.obs, N.get("tic_source", "peaks"))
        for m in N["methods"]:
            a.obs[f"size_factor_{m}"] = f[m]
        um = run_umaia(X, a.obs, N.get("umaia", {}), log) if N.get("umaia", {}).get("enabled") else None

        # One normalized copy at a time (the matrix can be a few GB), only section means are kept
        def layer(name, X=X, f=f, um=um):
            if name == "raw":
                return X
            if name == "umaia":
                return um
            return (X / f[name][:, None]).astype(np.float32)

        names = ["raw"] + list(N["methods"]) + (["umaia"] if um is not None else [])
        sec_info = a.obs.assign(key=a.obs["slide"].astype(str) + "|" + a.obs["section"].astype(str)) \
            .drop_duplicates("key").set_index("key")[["slide", "section", "genotype", "animal_id"]]
        rows, pcs = [], {}
        for name in names:
            means = section_means(layer(name), a.obs)
            r2 = slide_partial_r2(means, sec_info)
            cv = (means.std(axis=0) / means.mean(axis=0)).to_numpy()
            rows.append({"method": name, "median slide partial R2": float(np.nanmedian(r2)),
                         "median CV of section means": float(np.nanmedian(cv))})
            pcs[name] = (means.index, pca2(means.to_numpy()))
        cmp = pd.DataFrame(rows)
        cmp.to_csv(ctx.out("results", "06_normalize", mode, "method_comparison.csv"), index=False)

        # figures: section PCA before/after, size factors per section
        fig, axes = plt.subplots(1, len(pcs), figsize=(3.2 * len(pcs), 3.2), squeeze=False)
        slides = sorted(a.obs["slide"].unique())
        cmap = plt.get_cmap("tab10")
        for ax, (name, (idx, pc)) in zip(axes[0], pcs.items()):
            info = sec_info.loc[idx]
            for k, s in enumerate(slides):
                for g, mk in ((ctx.cfg["project"]["genotypes"]["control"], "o"),
                              (ctx.cfg["project"]["genotypes"]["mutant"], "^")):
                    sel = ((info["slide"] == s) & (info["genotype"] == g)).to_numpy()
                    ax.scatter(pc[sel, 0], pc[sel, 1], color=cmap(k), marker=mk, s=40,
                               label=s if g == ctx.cfg["project"]["genotypes"]["control"] else None)
            r2v = cmp.set_index("method").loc[name, "median slide partial R2"]
            ax.set_title(f"{name}\nslide partial R2 {r2v:.2f}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
        axes[0, 0].legend(fontsize=6, title="slide (o Control, ^ KO)", title_fontsize=6, loc="upper right")
        fig.suptitle(f"{mode}: PCA of section mean profiles (colour = slide)", fontsize=10)
        S.fig(fig, "slide_effect_pca.png")

        fig, axes = plt.subplots(1, len(N["methods"]), figsize=(4 * len(N["methods"]), 3), squeeze=False)
        key = a.obs["slide"].astype(str) + " " + a.obs["section"].str[0]
        for ax, m in zip(axes[0], N["methods"]):
            groups = sorted(key.unique())
            ax.boxplot([np.log2(f[m][(key == g).to_numpy()]) for g in groups], showfliers=False)
            ax.set_xticks(range(1, len(groups) + 1), groups, rotation=90, fontsize=7)
            ax.set_title(f"log2 size factor: {m}", fontsize=9)
        S.fig(fig, "size_factors.png")

        # write the normalized object: X = method for stats
        use = N.get("for_stats", "tic")
        if use not in names:
            S.warn(f"normalization.for_stats = {use} isn't available, falling back to tic", log)
            use = "tic"
        Xn = layer(use)
        out = ad.AnnData(X=Xn, obs=a.obs, var=a.var, uns=dict(a.uns))
        out.uns["normalization"] = use
        if um is not None and use != "umaia":
            out.layers["umaia"] = um
        emb = N.get("embedding_batch_correction", "none")
        Z = np.log2(Xn + 1e-3)
        Z = Z - Z.mean(axis=0)
        k = min(20, Z.shape[1] - 1)
        if k >= 2:
            # PCA axes learned on a random subset of pixels (fast), then applied to all of them
            sub = np.random.default_rng(0).choice(len(Z), min(len(Z), 20000), replace=False)
            comps = np.linalg.svd(Z[sub], full_matrices=False)[2][:k]
            out.obsm["X_pca"] = (Z @ comps.T).astype(np.float32)
            if emb == "harmony":
                try:
                    import harmonypy
                    ho = harmonypy.run_harmony(out.obsm["X_pca"], out.obs, "slide")
                    out.obsm["X_pca_harmony_PLOTS_ONLY"] = np.asarray(ho.Z_corr).T.astype(np.float32)
                except Exception as e:
                    S.warn(f"harmony requested but failed/unavailable: {e}", log)
        path = ctx.out("processed", mode, "matrix_norm.h5ad")
        out.write_h5ad(path, compression="gzip")
        log.info(f"wrote {path}")

        S.h("Method comparison (lower = less slide effect)")
        S.table(cmp.round(3))
        best = cmp[cmp.method != "raw"].sort_values("median slide partial R2").iloc[0]["method"]
        S.p(f"Used for statistics: **{use}** (fixed in the config). Lowest slide effect here: {best}. "
            "If they differ a lot, tell me before running 09; we shouldn't switch after seeing results.")
        S.p("Note: with Control always on top and Mutant always on the bottom, no per-pixel "
            "normalization can remove a top-to-bottom technical trend; see the 03 gradient check.")
        S.kv({"pixels": out.n_obs, "peaks": out.n_vars, "embedding batch correction (plots only)": emb})
        log.info(f"summary: {S.write()}")
        del a, X, Xn, Z, out
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
