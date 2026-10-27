"""Registering MSI sections to an Allen CCF coronal slice, and turning that into region masks.

The idea, step by step, for one section:
  1. FIXED image = the section, at MSI pixel size. Either a lipid image with clear anatomy
     (a sulfatide/HexCer peak, which lights up white matter) or the optical TIFF if you
     gave one for this slide.
  2. MOVING image = an atlas coronal slice, rendered as a "pseudo lipid image": white matter
     bright (1.0), grey matter dim (0.35), outside 0. Resampled to the MSI pixel size.
  3. For each candidate AP level (around Bregma -1.58 to -2.18) and starting rotation:
     moments-based start, similarity transform (rotate, scale, shift), then a full affine.
     We keep the AP level / rotation with the best normalized cross-correlation.
     Mutual information is the metric during the fit, so it doesn't matter whether the
     guide image has white matter bright or dark (TIFFs often have it dark).
  4. On the best one, a B-spline (non-linear) step on top of the affine, to soak up local
     differences (the KO brain may be atrophied, ventricles may be bigger, ...).
  5. Warp the atlas LABELS (nearest neighbour) onto the section grid. Each MSI pixel now
     has an Allen structure id; region groups (cortex, hippocampus, ...) include all
     descendants of the listed acronyms.
  6. QC: Dice overlap between the tissue mask and the warped atlas brain. Below
     registration.min_dice, the atlas regions for that section are not used.

There is always a fallback, too: white vs grey matter from the guide image itself (Otsu
split inside the tissue, the smaller class is white matter).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import filters, transform

WM_VALUE, GM_VALUE = 1.0, 0.35


class Atlas:
    def __init__(self, annotation, resolution_um: float, structures: pd.DataFrame, name: str):
        self.annotation = annotation          # (AP, DV, ML), structure ids
        self.res = float(resolution_um)
        self.structures = structures          # id, acronym, name, structure_id_path (list)
        self.name = name
        self._anc = {int(r.id): set(int(i) for i in r.structure_id_path) for r in structures.itertuples()}
        acr = dict(zip(structures.acronym, structures.id))
        self.acronym_to_id = acr
        self.id_to_acronym = dict(zip(structures.id, structures.acronym))
        wm_root = acr.get("fiber tracts")
        self.wm_ids = {i for i, anc in self._anc.items() if wm_root in anc} if wm_root else set()
        vs = acr.get("VS")   # ventricular systems
        self.vs_ids = {i for i, anc in self._anc.items() if vs in anc} if vs else set()

    def ids_under(self, acronyms: list[str]) -> set[int]:
        roots = {self.acronym_to_id[a] for a in acronyms if a in self.acronym_to_id}
        return {i for i, anc in self._anc.items() if anc & roots}

    def coronal(self, ap_index: int) -> np.ndarray:
        ap_index = int(np.clip(ap_index, 0, self.annotation.shape[0] - 1))
        return np.asarray(self.annotation[ap_index])


def load_atlas(ctx, log) -> Atlas | None:
    src = ctx.get("registration.atlas_source", "brainglobe")
    if src == "brainglobe":
        try:
            from brainglobe_atlasapi import BrainGlobeAtlas
            d = ctx.cache_root / "atlas"
            if not d.exists():
                log.warning("no cached atlas (run 00_setup.py --download)")
                return None
            a = BrainGlobeAtlas(ctx.get("registration.atlas_name"), brainglobe_dir=d, check_latest=False)
            st = pd.DataFrame([{"id": v["id"], "acronym": v["acronym"], "name": v["name"],
                                "structure_id_path": v["structure_id_path"]} for v in a.structures.values()])
            name = getattr(a, "atlas_name", ctx.get("registration.atlas_name"))
            log.info(f"atlas {name}: {a.annotation.shape}, {a.resolution} um, orientation {a.orientation}")
            if a.orientation != "asr":
                log.warning(f"atlas orientation is {a.orientation}, code assumes asr (AP, DV, ML)")
            return Atlas(a.annotation, a.resolution[0], st, name)
        except Exception as e:
            log.warning(f"could not load the brainglobe atlas ({e})")
            return None
    p = Path(src)
    if not p.is_absolute():
        p = ctx.config_path.parent / p
    if not p.exists():
        log.warning(f"atlas file {p} not found")
        return None
    z = np.load(p, allow_pickle=False)
    st = pd.DataFrame(json.loads(str(z["structures"])))
    log.info(f"atlas from file {p.name}: {z['annotation'].shape}, {float(z['resolution_um'])} um")
    return Atlas(z["annotation"], float(z["resolution_um"]), st, p.stem)


def pseudo_image(atlas: Atlas, labels: np.ndarray) -> np.ndarray:
    img = np.zeros(labels.shape, np.float32)
    img[labels > 0] = GM_VALUE
    if atlas.wm_ids:
        img[np.isin(labels, list(atlas.wm_ids))] = WM_VALUE
    if atlas.vs_ids:
        img[np.isin(labels, list(atlas.vs_ids))] = 0.0
    return img


def to_msi_scale(img: np.ndarray, atlas_res: float, pixel_um: float, order: int) -> np.ndarray:
    f = atlas_res / pixel_um
    out = transform.rescale(img.astype(np.float64), f, order=order, preserve_range=True,
                            anti_aliasing=(order > 0 and f < 1), mode="constant")
    return out.astype(img.dtype if order == 0 else np.float32)


def crop_to_content(img: np.ndarray, pad: int) -> tuple[np.ndarray, tuple[int, int]]:
    rr, cc = np.nonzero(img > 0)
    if not len(rr):
        return img, (0, 0)
    r0, r1 = max(0, rr.min() - pad), min(img.shape[0], rr.max() + pad + 1)
    c0, c1 = max(0, cc.min() - pad), min(img.shape[1], cc.max() + pad + 1)
    return img[r0:r1, c0:c1], (r0, c0)


def ncc(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    x, y = a[mask].astype(float), b[mask].astype(float)
    if x.std() == 0 or y.std() == 0:
        return -1.0
    return float(np.corrcoef(x, y)[0, 1])


def dice(a: np.ndarray, b: np.ndarray) -> float:
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else 0.0


def _register_affine(sitk, fixed, moving, angle_deg: float, iters: int = 200):
    init = sitk.CenteredTransformInitializer(fixed, moving, sitk.Similarity2DTransform(),
                                             sitk.CenteredTransformInitializerFilter.MOMENTS)
    init = sitk.Similarity2DTransform(init)
    init.SetAngle(np.deg2rad(angle_deg))

    def method():
        R = sitk.ImageRegistrationMethod()
        R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=24)
        R.SetMetricSamplingStrategy(R.NONE)
        R.SetInterpolator(sitk.sitkLinear)
        R.SetOptimizerAsRegularStepGradientDescent(learningRate=1.0, minStep=1e-4, numberOfIterations=iters,
                                                   relaxationFactor=0.6)
        R.SetOptimizerScalesFromPhysicalShift()
        R.SetShrinkFactorsPerLevel([2, 1])
        R.SetSmoothingSigmasPerLevel([1.0, 0.0])
        return R

    R = method()
    R.SetInitialTransform(init, inPlace=False)
    sim = R.Execute(fixed, moving)
    aff = sitk.AffineTransform(2)
    s = sitk.Similarity2DTransform(sim.GetNthTransform(0) if hasattr(sim, "GetNthTransform") else sim)
    aff.SetMatrix(s.GetMatrix())
    aff.SetTranslation(s.GetTranslation())
    aff.SetCenter(s.GetCenter())
    R = method()
    R.SetInitialTransform(aff, inPlace=False)
    out = R.Execute(fixed, moving)
    return out


def _register_bspline(sitk, fixed, moving, affine, grid_px: float, iters: int):
    size = fixed.GetSize()
    mesh = [max(1, int(round(s / max(grid_px, 2)))) for s in size]
    bs = sitk.BSplineTransformInitializer(fixed, mesh)
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=24)
    R.SetMetricSamplingStrategy(R.NONE)
    R.SetInterpolator(sitk.sitkLinear)
    R.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-5, numberOfIterations=iters,
                           maximumNumberOfCorrections=5, maximumNumberOfFunctionEvaluations=iters * 4)
    R.SetMovingInitialTransform(affine)
    R.SetInitialTransform(bs, inPlace=True)
    R.Execute(fixed, moving)
    comp = sitk.CompositeTransform(2)
    comp.AddTransform(affine)
    comp.AddTransform(bs)
    return comp


def resample(sitk, moving_np, fixed_shape, tfm, nearest: bool):
    mov = sitk.GetImageFromArray(moving_np.astype(np.float32 if not nearest else np.float64))
    ref = sitk.Image(int(fixed_shape[1]), int(fixed_shape[0]), sitk.sitkFloat32)
    out = sitk.Resample(mov, ref, tfm, sitk.sitkNearestNeighbor if nearest else sitk.sitkLinear, 0.0)
    arr = sitk.GetArrayFromImage(out)
    return arr.astype(np.int64) if nearest else arr


def register_section(guide: np.ndarray, tissue: np.ndarray, atlas: Atlas, ap_indices: list[int],
                     rotations: list[float], pixel_um: float, bspline_grid_mm: float, bspline_iters: int,
                     log=None) -> dict:
    """guide: 2D float image of ONE section (0 outside), tissue: its boolean mask."""
    import SimpleITK as sitk

    g = guide.astype(np.float32)
    fixed = sitk.GetImageFromArray(g)
    best = None
    tried = []
    for ap in ap_indices:
        lab = atlas.coronal(ap)
        lab_s = to_msi_scale(lab, atlas.res, pixel_um, order=0)
        lab_c, (r0, c0) = crop_to_content(lab_s, pad=3)      # same crop for labels and image
        mov_np = pseudo_image(atlas, lab_c)
        moving = sitk.GetImageFromArray(mov_np)
        for rot in rotations:
            try:
                tfm = _register_affine(sitk, fixed, moving, rot)
            except Exception as e:
                if log:
                    log.info(f"   AP {ap} rot {rot}: affine failed ({e})")
                continue
            warped = resample(sitk, mov_np, g.shape, tfm, nearest=False)
            score = ncc(g, warped, tissue | (warped > 0))
            tried.append({"ap_index": ap, "rotation": rot, "ncc": score})
            if best is None or score > best["ncc"]:
                best = {"ap_index": ap, "rotation": rot, "ncc": score, "affine": tfm, "moving": mov_np,
                        "labels": lab_c}
    if best is None:
        return {"ok": False, "reason": "all affine fits failed", "tried": tried}

    # affine matrix -> scale (sqrt det): how big the atlas had to be made to fit
    try:
        A = np.array(sitk.AffineTransform(best["affine"].GetNthTransform(0) if hasattr(best["affine"], "GetNthTransform")
                                          else best["affine"]).GetMatrix()).reshape(2, 2)
    except Exception:
        A = np.eye(2)
    # The transform maps fixed -> moving coordinates, so section size / atlas size = 1/sqrt(det)
    scale = float(1 / np.sqrt(abs(np.linalg.det(A)))) if np.linalg.det(A) != 0 else np.nan

    final = best["affine"]
    bs_ok = False
    try:
        grid_px = bspline_grid_mm * 1000 / pixel_um
        final = _register_bspline(sitk, fixed, sitk.GetImageFromArray(best["moving"]), best["affine"],
                                  grid_px, bspline_iters)
        bs_ok = True
    except Exception as e:
        if log:
            log.info(f"   B-spline step failed ({e}); keeping the affine result")
    lab_w = resample(sitk, best["labels"].astype(np.float64), g.shape, final, nearest=True)
    lab_aff = resample(sitk, best["labels"].astype(np.float64), g.shape, best["affine"], nearest=True)
    img_w = resample(sitk, best["moving"], g.shape, final, nearest=False)
    d = dice(tissue, lab_w > 0)
    d_aff = dice(tissue, lab_aff > 0)
    if d_aff > d + 0.02:      # B-spline made it worse: keep the affine one
        lab_w, d, bs_ok = lab_aff, d_aff, False
        img_w = resample(sitk, best["moving"], g.shape, best["affine"], nearest=False)
    return {"ok": True, "ap_index": best["ap_index"], "rotation": best["rotation"], "ncc_affine": best["ncc"],
            "ncc_final": ncc(g, img_w, tissue | (img_w > 0)), "dice": d, "dice_affine": d_aff,
            "bspline_used": bs_ok, "scale_section_vs_atlas": scale, "labels": lab_w, "warped_image": img_w,
            "tried": tried}


def wm_gm_from_guide(guide: np.ndarray, tissue: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fallback: Otsu split of the guide inside the tissue; the smaller class is white matter."""
    v = np.log1p(np.clip(guide[tissue], 0, None))
    if v.size < 10 or np.ptp(v) == 0:
        return np.zeros_like(tissue), tissue.copy()
    t = filters.threshold_otsu(v)
    hi = np.zeros_like(tissue)
    hi[tissue] = v > t
    hi = ndimage.binary_opening(hi, iterations=1) & tissue
    lo = tissue & ~hi
    wm, gm = (hi, lo) if hi.sum() <= lo.sum() else (lo, hi)
    return wm, gm
