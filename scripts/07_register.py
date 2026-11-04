"""Phase 07: register every section to the Allen CCF and make region masks.

    python scripts/07_register.py --config config.yaml

Per mode and per section (12 per mode):
  * guide image: the optical TIFF if registration.optical_tiffs gives one for the slide,
    otherwise a lipid image with clear anatomy (auto: strongest sulfatide, then HexCer,
    then the first principal component of all lipid images)
  * search AP levels around Bregma -1.58 to -2.18 (+ a margin), affine, then B-spline
  * warp atlas labels onto the section -> per-pixel Allen structure
  * region masks from registration.regions (descendants included)
  * fallback (always computed): white vs grey matter from the guide image itself

Bregma to CCF is approximate (registration.bregma_ccf_ap_mm). The AP level that fits best is
reported per section; if all 12 sections land at very different levels, tell me.

Writes processed/<mode>/regions.parquet (one row per pixel of matrix_norm.h5ad, boolean
region columns), results/07_register/<mode>/registration.csv, overlays, summary.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, register  # noqa: E402
from cln3msi.util import read_h5ad  # noqa: E402
from cln3msi.images import scaled01, section_image  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402


def choose_guide(a: ad.AnnData, how, log) -> tuple[np.ndarray, str]:
    """Per-pixel guide values for the whole mode, plus a description."""
    var = a.var
    if how not in (None, "auto"):
        j = int(np.argmin(np.abs(var["mz"].to_numpy() - float(how))))
        return np.asarray(a.X[:, j]).ravel(), f"m/z {var['mz'].iloc[j]:.4f} (from config)"
    if "tag" in var:
        ok = ~var.get("is_isotope", pd.Series(False, index=var.index)).fillna(False).astype(bool)
        for tag in ("SHexCer", "HexCer"):
            cand = var[ok & (var["tag"] == tag)]
            if len(cand):
                pid = cand["median_height"].idxmax()
                j = var.index.get_loc(pid)
                return np.asarray(a.X[:, j]).ravel(), f"{var.loc[pid, 'label']} (strongest {tag})"
    Z = np.log2(np.asarray(a.X) + 1e-3)
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    sub = np.random.default_rng(0).choice(len(Z), min(len(Z), 20000), replace=False)
    v = np.linalg.svd(Z[sub], full_matrices=False)[2][0]
    log.info("no sulfatide/HexCer annotated: guide = PC1 of all peaks")
    pc = Z @ v
    return pc - pc.min(), "PC1 of all peaks (no sulfatide/HexCer annotated)"


def tiff_guide(ctx, slide: str, shape, log) -> np.ndarray | None:
    tiffs = ctx.get("registration.optical_tiffs") or {}
    t = tiffs.get(slide)
    if not t:
        return None
    p = Path(t)
    if not p.is_absolute():
        p = ctx.data_root / p
    if not p.exists():
        log.warning(f"TIFF for slide {slide} not found: {p}")
        return None
    import tifffile
    from skimage import transform
    img = tifffile.imread(p).astype(np.float32)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=-1) if img.shape[-1] in (3, 4) else img[0]
    um = ctx.get("registration.tiff_um_per_px")
    pix = float(ctx.get("acquisition.pixel_size_um", 50))
    if um:
        img = transform.rescale(img, float(um) / pix, anti_aliasing=True, preserve_range=True)
        img = img[: shape[0], : shape[1]]
        pad = [(0, max(0, shape[0] - img.shape[0])), (0, max(0, shape[1] - img.shape[1]))]
        img = np.pad(img, pad)
    else:
        img = transform.resize(img, shape, anti_aliasing=True, preserve_range=True)
    log.info(f"slide {slide}: using TIFF {p.name} as guide (assumed to cover the MSI field of view)")
    return img.astype(np.float32)


def reuse_registration(ctx, log, S, mode: str, src: str, obs: pd.DataFrame, gvals: np.ndarray,
                       gdesc: str, R: dict) -> bool:
    """Copy another mode's per-pixel regions when both modes share the same pixel grid.

    On these slides the three modes were exported on identical rasters (same pixel count and
    coordinates per slide), so one good registration, from the mode with the clearest
    white-matter marker, serves all three. Pixels are matched by (animal, x, y). If fewer
    than 90% of this mode's pixels have a partner, we register this mode on its own instead.
    Returns True if the regions were reused.
    """
    rp = ctx.output_root / "processed" / src / "regions.parquet"
    if not rp.exists():
        S.warn(f"registration.reuse_from = {src}, but {src} has no regions.parquet yet; "
               "registering this mode on its own", log)
        return False
    src_reg = pd.read_parquet(rp)
    matched = obs.index.isin(src_reg.index)
    frac = float(matched.mean())
    if frac < 0.9:
        S.warn(f"only {frac:.0%} of the {mode} pixels exist in {src}: the pixel grids differ, "
               "so this mode is registered on its own", log)
        return False

    regions = src_reg.reindex(obs.index)
    for c in regions.columns:
        if c.startswith("region_") or c == "registration_ok":
            regions[c] = regions[c].astype("boolean").fillna(False).astype(bool)
    regions["atlas_id"] = regions["atlas_id"].fillna(0).astype(np.int64)
    regions["atlas_acronym"] = regions["atlas_acronym"].fillna("").astype(str)
    regions["region_whole_section"] = True
    regions.to_parquet(ctx.out("processed", mode, "regions.parquet"))

    sp = ctx.output_root / "results" / "07_register" / src / "registration.csv"
    reg = pd.read_csv(sp) if sp.exists() else pd.DataFrame()
    if len(reg):
        reg.insert(0, "reused_from", src)
    reg.to_csv(ctx.out("results", "07_register", mode, "registration.csv"), index=False)

    # overlays: this mode's guide image + the reused region outlines
    names = [c for c in regions.columns if c.startswith("region_")
             and c not in ("region_whole_section", "region_wm_lipid", "region_gm_lipid")]
    keys = list(obs.groupby(["slide", "section"]).groups)
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow), squeeze=False)
    colors = plt.get_cmap("tab10")
    for ax, (slide, section) in zip(axes.ravel(), keys):
        img, mask, box, _, _ = section_image(obs, gvals, slide, section, pad=4, fill=0)
        g = np.zeros_like(img)
        g[mask] = scaled01(img[mask])
        ax.imshow(np.where(mask, g, np.nan), cmap="Greys_r", interpolation="nearest")
        for k, name in enumerate(names):
            m, *_ = section_image(obs, regions[name].to_numpy().astype(float), slide, section, pad=4, fill=0, box=box)
            if np.nanmax(m) > 0:
                ax.contour(m, levels=[0.5], colors=[colors(k)], linewidths=0.8)
        ax.set_title(f"{slide} {section}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
    for ax in axes.ravel()[len(keys):]:
        ax.axis("off")
    handles = [plt.Line2D([], [], color=colors(k), label=n.replace("region_", "")) for k, n in enumerate(names)]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 7), fontsize=7)
    fig.suptitle(f"{mode}: regions reused from {src}, drawn on {mode}'s own guide ({gdesc})", fontsize=9)
    S.fig(fig, "overlays.png")

    S.h("Settings")
    S.kv({"regions": f"reused from {src} (registration.reuse_from)",
          f"{mode} pixels with a partner in {src}": f"{frac:.1%}"})
    S.h("Per section (copied from the source mode)")
    S.table(reg.round(3))
    cnt = regions.filter(like="region_").groupby(obs["animal_id"]).sum()
    S.h("Pixels per region and animal")
    S.table(cnt.reset_index())
    log.info(f"{mode}: regions reused from {src} ({frac:.1%} of pixels matched)")
    return True


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("07_register", __doc__)
    R = ctx.cfg["registration"]
    atlas = register.load_atlas(ctx, log)
    pix = float(ctx.get("acquisition.pixel_size_um", 50))
    if atlas is not None:
        b0, b1 = R["ap_search_bregma_mm"]
        steps = np.arange(min(b0, b1), max(b0, b1) + 1e-9, float(R["ap_step_mm"]))
        ap_idx = sorted({int(round((float(R["bregma_ccf_ap_mm"]) - b) * 1000 / atlas.res)) for b in steps})
        ap_idx = [i for i in ap_idx if 0 <= i < atlas.annotation.shape[0]]
        log.info(f"AP indices searched: {ap_idx}")
        region_ids = {name: atlas.ids_under(acrs) for name, acrs in R["regions"].items()}
        for name, ids in region_ids.items():
            if not ids:
                log.warning(f"region {name}: none of {R['regions'][name]} found in the atlas")

    for mode in modes:
        S = Summary(ctx, "07_register", mode, f"Registration and regions, {mode}")
        p = ctx.out("processed", mode, "matrix_norm.h5ad", mkdir=False)
        if not p.exists():
            S.warn("no matrix_norm.h5ad, run 06 first", log)
            S.write()
            continue
        a = read_h5ad(p)
        obs = a.obs
        guide_cfg = R.get("guide_ion", "auto")
        how = guide_cfg.get(mode, "auto") if isinstance(guide_cfg, dict) else guide_cfg
        gvals, gdesc = choose_guide(a, how, log)
        log.info(f"{mode}: guide = {gdesc}")
        src = R.get("reuse_from")
        if src and src != mode and reuse_registration(ctx, log, S, mode, src, obs, gvals, gdesc, R):
            log.info(f"summary: {S.write()}")
            continue

        regions = pd.DataFrame(index=obs.index)
        regions["atlas_id"] = 0
        regions["atlas_acronym"] = ""
        regions["registration_ok"] = False
        for name in R["regions"]:
            regions[f"region_{name}"] = False
        regions["region_wm_lipid"] = False
        regions["region_gm_lipid"] = False
        regions["region_whole_section"] = True

        rows, panels = [], []
        for slide in sorted(obs["slide"].unique()):
            z = np.load(ctx.out("processed", mode, "qc", f"sections_{slide}.npz", mkdir=False))
            tif = tiff_guide(ctx, slide, z["positions"].shape, log)
            for section in ("Top", "Bottom"):
                if not ((obs.slide == slide) & (obs.section == section)).any():
                    continue
                img, mask, box, rows_idx, (yy, xx) = section_image(obs, gvals, slide, section, pad=4, fill=0)
                if tif is not None:
                    r0, r1, c0, c1 = box
                    g = np.zeros((r1 - r0, c1 - c0), np.float32)      # box can stick out of the image
                    part = tif[r0:min(r1, tif.shape[0]), c0:min(c1, tif.shape[1])]
                    g[:part.shape[0], :part.shape[1]] = part
                    g = scaled01(g)
                    g[~mask] = 0
                    src = "TIFF"
                else:
                    g = np.zeros_like(img)
                    g[mask] = scaled01(img[mask])
                    src = gdesc
                wm, gm = register.wm_gm_from_guide(np.where(mask, img, 0), mask)
                pix_ids = regions.index[rows_idx]
                regions.loc[pix_ids, "region_wm_lipid"] = wm[yy, xx]
                regions.loc[pix_ids, "region_gm_lipid"] = gm[yy, xx]
                rec = {"slide": slide, "section": section, "animal_id": obs.loc[pix_ids[0], "animal_id"],
                       "genotype": obs.loc[pix_ids[0], "genotype"], "guide": src,
                       "wm_lipid_fraction": float(wm[mask].mean())}
                res = None
                if atlas is not None:
                    log.info(f"registering {mode} {slide} {section} ...")
                    res = register.register_section(g, mask, atlas, ap_idx, [float(r) for r in R["try_rotations_deg"]],
                                                    pix, float(R["bspline_grid_mm"]), int(R["bspline_iterations"]), log)
                    if res["ok"]:
                        ok = res["dice"] >= float(R["min_dice"])
                        lab = res["labels"]
                        ids = lab[yy, xx]
                        regions.loc[pix_ids, "atlas_id"] = ids
                        regions.loc[pix_ids, "atlas_acronym"] = [atlas.id_to_acronym.get(int(i), "") for i in ids]
                        regions.loc[pix_ids, "registration_ok"] = ok
                        if ok:
                            for name, rid in region_ids.items():
                                regions.loc[pix_ids, f"region_{name}"] = np.isin(ids, list(rid))
                        ap_mm = res["ap_index"] * atlas.res / 1000
                        rec.update({"ap_index": res["ap_index"],
                                    "bregma_estimate_mm": round(float(R["bregma_ccf_ap_mm"]) - ap_mm, 3),
                                    "rotation_start": res["rotation"], "ncc_affine": res["ncc_affine"],
                                    "ncc_final": res["ncc_final"], "dice": res["dice"], "bspline_used": res["bspline_used"],
                                    "scale_section_vs_atlas": res["scale_section_vs_atlas"], "used": ok})
                    else:
                        rec.update({"used": False, "note": res.get("reason")})
                rows.append(rec)
                panels.append((slide, section, g, mask, wm, res, rec))
            log.mem(f"after slide {slide}")

        reg = pd.DataFrame(rows)
        reg.to_csv(ctx.out("results", "07_register", mode, "registration.csv"), index=False)
        regions.to_parquet(ctx.out("processed", mode, "regions.parquet"))

        # overlays: guide image + contours of the region groups + lipid WM outline
        ncol = 4
        nrow = int(np.ceil(len(panels) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow), squeeze=False)
        colors = plt.get_cmap("tab10")
        for ax, (slide, section, g, mask, wm, res, rec) in zip(axes.ravel(), panels):
            ax.imshow(np.where(mask, g, np.nan), cmap="Greys_r", interpolation="nearest")
            if res is not None and res.get("ok"):
                lab = res["labels"]
                for k, (name, rid) in enumerate(region_ids.items()):
                    m = np.isin(lab, list(rid)).astype(float)
                    if m.any():
                        ax.contour(m, levels=[0.5], colors=[colors(k)], linewidths=0.9)
                ax.contour((lab > 0).astype(float), levels=[0.5], colors="yellow", linewidths=0.6, linestyles=":")
                t = f"Dice {res['dice']:.2f}, Bregma~{rec.get('bregma_estimate_mm')}"
            else:
                ax.contour(wm.astype(float), levels=[0.5], colors="cyan", linewidths=0.8)
                t = "no atlas: lipid WM outline (cyan)"
            ax.set_title(f"{slide} {section} {rec['animal_id']}\n{t}", fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
        for ax in axes.ravel()[len(panels):]:
            ax.axis("off")
        if atlas is not None:
            handles = [plt.Line2D([], [], color=colors(k), label=n) for k, n in enumerate(region_ids)]
            handles.append(plt.Line2D([], [], color="yellow", ls=":", label="atlas brain outline"))
            fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=7)
        fig.suptitle(f"{mode}: atlas regions warped onto each section (guide: {gdesc})", fontsize=9)
        S.fig(fig, "overlays.png")

        fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.5 * nrow), squeeze=False)
        for ax, (slide, section, g, mask, wm, res, rec) in zip(axes.ravel(), panels):
            show = np.full(mask.shape, np.nan)
            show[mask] = 0.4
            show[wm] = 1.0
            ax.imshow(show, cmap="Greys", vmin=0, vmax=1, interpolation="nearest")
            ax.set_title(f"{slide} {section}: WM {rec['wm_lipid_fraction']:.0%}", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
        for ax in axes.ravel()[len(panels):]:
            ax.axis("off")
        fig.suptitle(f"{mode}: fallback white (black) vs grey matter from the guide image", fontsize=9)
        S.fig(fig, "wm_gm_fallback.png")

        S.h("Settings")
        S.kv({"atlas": atlas.name if atlas is not None else "none (fallback WM/GM only)", "guide": gdesc,
              "min Dice to use atlas regions": R["min_dice"]})
        S.h("Per section")
        S.table(reg.round(3))
        if atlas is not None and "used" in reg:
            n_ok = int(reg["used"].fillna(False).sum())
            S.p(f"{n_ok}/{len(reg)} sections passed the Dice threshold and get atlas regions.")
            if n_ok < len(reg):
                S.warn(f"{len(reg) - n_ok} sections failed registration QC; they only get the lipid WM/GM "
                       "regions (and whole_section).", log)
            if "scale_section_vs_atlas" in reg:
                sc = reg.groupby("genotype")["scale_section_vs_atlas"].median()
                S.p("Median scale section/atlas by genotype (descriptive, <1 = smaller than atlas): "
                    + ", ".join(f"{k} {v:.2f}" for k, v in sc.items()))
        cnt = regions.filter(like="region_").groupby(obs["animal_id"]).sum()
        S.h("Pixels per region and animal")
        S.table(cnt.reset_index())
        S.p("Look at overlays.png: contours should follow the corpus callosum and the hippocampus. If the "
            "AP levels are off, adjust registration.ap_search_bregma_mm or bregma_ccf_ap_mm.")
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
