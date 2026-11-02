"""Tissue vs background, and splitting a slide into its Top and Bottom sections.

The recipe:
  1. Take a per-pixel signal image (by default the TIC without the DAN matrix ions, because
     off-tissue pixels can be full of matrix signal and look "bright" in a plain TIC).
  2. log, smooth a little, Otsu threshold on the acquired pixels.
  3. Clean up: close small gaps, drop specks, fill holes.
  4. Work out the split direction (qc.sections.split_axis): "rows" if the two sections sit
     one above the other in the image, "cols" if they sit side by side. "auto" looks at the
     two biggest pieces (or the shape of the tissue if it's one piece) and picks the
     direction they're separated along.
  5. Split into sections along that direction:
       - if the config gives a manual split index for this slide, use it
       - if there are exactly 2 big pieces, those are the sections
       - if there are more (a section broke apart), group pieces by position
       - if there's only 1 (sections touch), cut at the emptiest line in the middle
  6. Name them: the section with the smaller index along the split direction (higher up for
     rows, further left for cols) is "Top" by default. Flip with
     qc.sections.top_is_low_index, or per slide with qc.sections.flip_slides.

"Top" and "Bottom" are the sample sheet's names for the two positions on the slide. When the
slide was imaged rotated, they show up as left/right in the image, so always check the
labelled figure (sections_<slide>.png) against what Pilar says.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from skimage import filters, measure, morphology

TOP, BOTTOM = 1, 2
POSITION_NAME = {TOP: "Top", BOTTOM: "Bottom"}


def _clean(m: np.ndarray, acquired: np.ndarray, min_object_px: int, close_radius: int) -> np.ndarray:
    if close_radius > 0:
        m = morphology.binary_closing(m, morphology.disk(close_radius))
    # drop specks smaller than min_object_px (done by hand: the skimage API changed in 0.26)
    lab, n = ndimage.label(m)
    if n:
        sizes = np.bincount(lab.ravel())
        m = (sizes >= max(1, min_object_px))[lab] & (lab > 0)
    return ndimage.binary_fill_holes(m) & acquired


def tissue_mask(img: np.ndarray, smooth_sigma: float, min_object_px: int, close_radius: int,
                method: str = "auto") -> tuple[np.ndarray, float, str]:
    """img: 2D signal, NaN where nothing was acquired. Returns (mask, log threshold, method used).

    Why "auto": Otsu splits the pixels into the two most different groups. When there's a
    good amount of empty glass, those groups are tissue vs glass, as intended. But when the
    raster hugs the tissue (only a thin margin of glass), the two most different groups are
    bright vs dim TISSUE, e.g. white matter vs grey matter in lipid negative mode, and whole
    regions of tissue get thrown out as "background" (holes, rings).
    So "auto" also tries a lower threshold (the lowest of a 3-group Otsu). If the extra pixels
    it adds lie mostly INSIDE the sections (inside their convex outline, i.e. they fill the
    holes), the lower threshold is the right one; if they lie outside, they're glass and plain
    Otsu stays. Force one with qc.tissue.threshold_method: otsu | multiotsu | triangle.
    """
    acquired = np.isfinite(img)
    x = np.log1p(np.where(acquired, np.clip(img, 0, None), 0))
    fill = np.nanmin(np.where(acquired, x, np.nan)) if acquired.any() else 0
    x = np.where(acquired, x, fill)
    if smooth_sigma > 0:
        x = ndimage.gaussian_filter(x, smooth_sigma)
    vals = x[acquired]
    if np.ptp(vals) == 0:
        thr = float(vals.max())
        return _clean((x > thr) & acquired, acquired, min_object_px, close_radius), thr, "constant"

    def low_threshold() -> tuple[float, str]:
        try:
            return float(filters.threshold_multiotsu(vals, classes=3)[0]), "multiotsu"
        except Exception:
            return float(filters.threshold_triangle(vals)), "triangle"

    if method == "otsu":
        thr = float(filters.threshold_otsu(vals))
        return _clean((x > thr) & acquired, acquired, min_object_px, close_radius), thr, "otsu"
    if method == "multiotsu":
        thr, used = low_threshold()
        return _clean((x > thr) & acquired, acquired, min_object_px, close_radius), thr, used
    if method == "triangle":
        thr = float(filters.threshold_triangle(vals))
        return _clean((x > thr) & acquired, acquired, min_object_px, close_radius), thr, "triangle"

    # auto
    t_otsu = float(filters.threshold_otsu(vals))
    m_otsu = _clean((x > t_otsu) & acquired, acquired, min_object_px, close_radius)
    t_low, low_name = low_threshold()
    if not m_otsu.any() or t_low >= t_otsu:
        return m_otsu, t_otsu, "otsu"
    m_low = _clean((x > t_low) & acquired, acquired, min_object_px, close_radius)
    added = m_low & ~m_otsu
    if added.sum() < 0.05 * m_otsu.sum():
        return m_otsu, t_otsu, "otsu"
    hull = morphology.convex_hull_object(m_otsu, connectivity=2)
    inside = float((added & hull).sum() / added.sum())
    if inside > 0.6:
        return m_low, t_low, f"auto: {low_name} (otsu left holes inside the tissue; {inside:.0%} of the extra pixels were inside the sections)"
    return m_otsu, t_otsu, "otsu"


def choose_axis(mask: np.ndarray, n_expected: int = 2) -> str:
    """'rows' if the sections are stacked vertically in the image, 'cols' if side by side."""
    lab = measure.label(mask, connectivity=2)
    props = sorted(measure.regionprops(lab), key=lambda r: -r.area)
    if not props:
        return "rows"
    big = [p for p in props if p.area >= 0.05 * props[0].area]
    if len(big) >= n_expected:
        # spread of the big pieces' centroids, area-weighted, in each direction
        cy = np.array([p.centroid[0] for p in big])
        cx = np.array([p.centroid[1] for p in big])
        w = np.array([p.area for p in big], float)
        sy = np.sqrt(np.average((cy - np.average(cy, weights=w)) ** 2, weights=w))
        sx = np.sqrt(np.average((cx - np.average(cx, weights=w)) ** 2, weights=w))
        return "rows" if sy >= sx else "cols"
    r0, c0, r1, c1 = big[0].bbox            # one piece: cut across its longer side
    return "rows" if (r1 - r0) >= (c1 - c0) else "cols"


def _split_rows(mask: np.ndarray, n_expected: int, manual: int | None) -> tuple[np.ndarray, str]:
    """Split along rows. Labels 1..n by increasing row (top first)."""
    labels = np.zeros(mask.shape, dtype=np.int8)
    rows = np.arange(mask.shape[0])[:, None]
    if manual is not None:
        labels[mask & (rows < manual)] = 1
        labels[mask & (rows >= manual)] = 2
        return labels, f"manual split at index {manual}"

    lab = measure.label(mask, connectivity=2)
    props = sorted(measure.regionprops(lab), key=lambda r: -r.area)
    if not props:
        return labels, "no tissue found"
    big = [p for p in props if p.area >= 0.05 * props[0].area]

    if len(big) >= n_expected:
        if len(big) == n_expected:
            groups = [[p] for p in big]
            how = f"{n_expected} connected pieces"
        else:
            # A section broke into pieces. Sort the pieces by position and cut at the
            # largest gap, which is the space between the two sections.
            big.sort(key=lambda p: p.centroid[0])
            cy = np.array([p.centroid[0] for p in big])
            cut = int(np.argmax(np.diff(cy))) + 1
            groups = [big[:cut], big[cut:]]
            how = f"{len(big)} pieces grouped by position (cut after piece {cut})"
        groups.sort(key=lambda g: np.average([p.centroid[0] for p in g], weights=[p.area for p in g]))
        for k, g in enumerate(groups[:n_expected], start=1):
            for p in g:
                labels[lab == p.label] = k
        # small leftovers (< 5% of the largest piece) are left as background on purpose
        return labels, how

    # Only one big piece: the two sections touch. Cut at the emptiest line in its middle half.
    p = big[0]
    r0, _, r1, _ = p.bbox
    profile = mask[r0:r1].sum(axis=1).astype(float)
    a, b = int(0.25 * len(profile)), int(0.75 * len(profile))
    cut = r0 + a + int(np.argmin(ndimage.uniform_filter1d(profile, 3)[a:b]))
    piece = lab == p.label
    labels[piece & (rows < cut)] = 1
    labels[piece & (rows >= cut)] = 2
    return labels, f"one connected piece, cut at the emptiest line ({cut})"


def split_sections(mask: np.ndarray, n_expected: int = 2, manual_split: int | None = None,
                   axis: str = "rows") -> tuple[np.ndarray, str]:
    """Return (labels, how): 0 background, 1..n by increasing index along `axis`."""
    if axis == "cols":
        lab, how = _split_rows(mask.T, n_expected, manual_split)
        return lab.T.copy(), how
    return _split_rows(mask, n_expected, manual_split)


def assign_positions(labels: np.ndarray, top_is_low_index: bool, flip: bool) -> np.ndarray:
    """Map labels 1/2 (by increasing index) to TOP/BOTTOM codes, honouring the config flips."""
    out = np.zeros_like(labels)
    first_is_top = top_is_low_index != flip         # XOR: either switch turns it around
    out[labels == 1] = TOP if first_is_top else BOTTOM
    out[labels == 2] = BOTTOM if first_is_top else TOP
    return out


def eroded(labels: np.ndarray, px: int) -> np.ndarray:
    """Same labels with `px` pixels peeled off each section's edge (edge pixels are often odd)."""
    if px <= 0:
        return labels.copy()
    out = np.zeros_like(labels)
    for v in np.unique(labels[labels > 0]):
        out[ndimage.binary_erosion(labels == v, iterations=px)] = v
    return out


def segment_slide(signal_img: np.ndarray, slide: str, qc_cfg: dict) -> dict:
    """Run the whole recipe for one slide with the qc: block of the config."""
    t = qc_cfg["tissue"]
    s = qc_cfg["sections"]
    mask, thr, method = tissue_mask(signal_img, t["smooth_sigma_px"], t["min_object_px"], t["close_radius_px"],
                                    t.get("threshold_method", "auto"))
    axis = s.get("split_axis", "auto")
    if axis not in ("rows", "cols"):
        axis = choose_axis(mask, s.get("n_expected", 2))
    # new key names, with the old ones still accepted
    manual = (s.get("manual_split") or s.get("manual_split_row") or {}).get(slide)
    top_low = s.get("top_is_low_index", s.get("top_is_low_row_index", True))
    lab, how = split_sections(mask, s.get("n_expected", 2), manual, axis)
    how = f"{how}, sections {'one above the other' if axis == 'rows' else 'side by side'}"
    pos = assign_positions(lab, top_low, slide in (s.get("flip_slides") or []))
    core = eroded(pos, t.get("erode_edge_px", 0))
    n_found = len(np.unique(lab[lab > 0]))
    return {"mask": mask, "threshold_log": thr, "positions": pos, "core": core, "how": how,
            "n_sections": n_found, "axis": axis, "threshold_method": method}
