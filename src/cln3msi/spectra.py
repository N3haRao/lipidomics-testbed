"""Spectrum-level number crunching: streaming sums, mean spectra, peak picking, alignment.

Everything here works one spectrum (or one block of spectra) at a time, so memory stays
at "a few mean spectra + a few per-pixel vectors", never "all spectra of a slide".

Words used below
  axis      the m/z values a mean spectrum lives on. For continuous files that's the file's
            own shared axis; for processed files we bin onto a log-spaced grid (constant
            width in ppm), because every spectrum has its own m/z values.
  window    an [lo, hi] m/z interval. Peak intensity in a pixel = sum of the intensities
            that fall inside the window. Works the same for profile and centroid data.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage, signal


# ---------------------------------------------------------------------------------------
# Binning grid for processed data
# ---------------------------------------------------------------------------------------
class PpmGrid:
    """Log-spaced bins: bin k covers [lo*(1+b)^k, lo*(1+b)^(k+1)), b = bin_ppm*1e-6."""

    def __init__(self, mz_lo: float, mz_hi: float, bin_ppm: float):
        self.lo, self.hi = float(mz_lo), float(mz_hi)
        self.step = np.log1p(bin_ppm * 1e-6)
        self.n = int(np.ceil(np.log(self.hi / self.lo) / self.step)) + 1
        self.centers = self.lo * np.exp((np.arange(self.n) + 0.5) * self.step)

    def index(self, mz: np.ndarray) -> np.ndarray:
        k = np.floor(np.log(np.maximum(mz, 1e-9) / self.lo) / self.step).astype(np.int64)
        return k

    def add(self, acc: np.ndarray, mz: np.ndarray, it: np.ndarray) -> None:
        k = self.index(mz)
        ok = (k >= 0) & (k < self.n)
        acc += np.bincount(k[ok], weights=it[ok], minlength=self.n)

    def add_profile(self, acc: np.ndarray, mz: np.ndarray, it: np.ndarray) -> None:
        """Profile data: interpolate the curve onto the grid instead of summing into bins.

        Summing would make a comb of spikes when the grid is finer than the data points.
        Processed profile spectra leave out empty stretches, so we never interpolate across
        a gap bigger than 3x the typical point spacing (those grid points stay 0).
        Only grid points that fall between two neighbouring data points are touched, so the
        cost scales with the data, not with the size of the grid.
        """
        if len(mz) < 2:
            return
        gap = np.diff(mz)
        ok = (gap > 0) & (gap <= 3 * np.median(gap))
        m0, m1 = mz[:-1][ok], mz[1:][ok]
        y0, y1 = it[:-1][ok].astype(np.float64), it[1:][ok].astype(np.float64)
        # first grid centre >= m0 and last grid centre < m1, in grid-index units
        u0 = np.log(m0 / self.lo) / self.step - 0.5
        u1 = np.log(m1 / self.lo) / self.step - 0.5
        k_lo = np.ceil(u0).astype(np.int64)
        k_hi = np.ceil(u1).astype(np.int64) - 1
        n = np.clip(k_hi - k_lo + 1, 0, None)
        if n.sum() == 0:
            return
        seg = np.repeat(np.arange(len(n)), n)
        start = np.repeat(np.cumsum(n) - n, n)
        ks = k_lo[seg] + (np.arange(n.sum()) - start)
        keep = (ks >= 0) & (ks < self.n)
        ks, seg = ks[keep], seg[keep]
        t = (self.centers[ks] - m0[seg]) / (m1[seg] - m0[seg])
        acc += np.bincount(ks, weights=y0[seg] + t * (y1[seg] - y0[seg]), minlength=self.n)


# ---------------------------------------------------------------------------------------
# Window sums: the core "how much of this ion is in this pixel" operation
# ---------------------------------------------------------------------------------------
def window_sums(mz: np.ndarray, it: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Sum of intensities with lo <= m/z <= hi, for many windows at once (cumsum trick)."""
    cs = np.concatenate([[0.0], np.cumsum(it, dtype=np.float64)])
    a = np.searchsorted(mz, lo, side="left")
    b = np.searchsorted(mz, hi, side="right")
    return cs[b] - cs[a]


def window_sums_block(axis: np.ndarray, block: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Same as window_sums for a (n_pixels, L) block that shares one m/z axis."""
    a = np.searchsorted(axis, lo, side="left")
    b = np.searchsorted(axis, hi, side="right")
    cs = np.zeros((block.shape[0], block.shape[1] + 1), dtype=np.float64)
    np.cumsum(block, axis=1, dtype=np.float64, out=cs[:, 1:])
    return cs[:, b] - cs[:, a]


def ppm_windows(centers: np.ndarray, ppm_half: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    half = np.asarray(centers) * np.asarray(ppm_half) * 1e-6
    return centers - half, centers + half


# ---------------------------------------------------------------------------------------
# Streaming pass 1: per-pixel TIC and per-pixel window images
# ---------------------------------------------------------------------------------------
def stream_pixel_sums(reader, lo: np.ndarray, hi: np.ndarray, chunk: int, log=None, every=10000):
    """One pass over a slide: returns (tic[n], win[n, n_windows])."""
    n = reader.n
    tic = np.zeros(n, dtype=np.float64)
    win = np.zeros((n, len(lo)), dtype=np.float64)
    if reader.shares_mz_axis:
        axis = reader.mz_axis()
        done = 0
        for idx, block in reader.iter_blocks(chunk=chunk):
            tic[idx] = block.sum(axis=1, dtype=np.float64)
            if len(lo):
                win[idx] = window_sums_block(axis, block, lo, hi)
            done += len(idx)
            if log is not None and (done // chunk) % max(1, every // chunk) == 0:
                log.mem(f"{reader.path.name}: {done}/{n} spectra")
    else:
        for i, mz, it in reader.iter_spectra(log=log, every=every):
            tic[i] = it.sum(dtype=np.float64)
            if len(lo):
                win[i] = window_sums(mz, it, lo, hi)
    return tic, win


# ---------------------------------------------------------------------------------------
# Streaming pass 2: mean spectra per pixel group (sections, off-tissue)
# ---------------------------------------------------------------------------------------
def stream_mean_spectra(reader, labels: np.ndarray, groups: list, grid: PpmGrid | None,
                        chunk: int, log=None, every=10000, max_pixels_per_group: int | None = None,
                        rng=None, profile: bool = False) -> tuple[np.ndarray, dict]:
    """Mean spectrum of each pixel group. labels[i] = group of spectrum i (or None/-1 to skip).

    Returns (axis, {group: mean_spectrum}). For continuous files the axis is the native one.
    max_pixels_per_group subsamples pixels evenly at random, to keep QC passes quick.
    """
    labels = np.asarray(labels, dtype=object)
    rng = np.random.default_rng(0) if rng is None else rng
    chosen = {}
    for g in groups:
        idx = np.flatnonzero(labels == g)
        if max_pixels_per_group and len(idx) > max_pixels_per_group:
            idx = np.sort(rng.choice(idx, max_pixels_per_group, replace=False))
        chosen[g] = idx
    if reader.shares_mz_axis:
        axis = reader.mz_axis()
        acc = {g: np.zeros(len(axis)) for g in groups}
        for g in groups:
            for _, block in reader.iter_blocks(indices=chosen[g], chunk=chunk):
                acc[g] += block.sum(axis=0, dtype=np.float64)
    else:
        assert grid is not None, "processed data needs a PpmGrid"
        axis = grid.centers
        acc = {g: np.zeros(grid.n) for g in groups}
        order = np.concatenate([chosen[g] for g in groups]) if groups else np.array([], int)
        lab_of = {int(i): g for g in groups for i in chosen[g]}
        for i, mz, it in reader.iter_spectra(indices=np.sort(order), log=log, every=every):
            (grid.add_profile if profile else grid.add)(acc[lab_of[i]], mz, it)
    means = {g: acc[g] / max(1, len(chosen[g])) for g in groups}
    return axis, means


# ---------------------------------------------------------------------------------------
# Peak picking on a mean spectrum
# ---------------------------------------------------------------------------------------
def noise_level(spec: np.ndarray, window: int) -> tuple[np.ndarray, np.ndarray]:
    """Baseline = rolling median; noise = 1.4826 x rolling median of |spec - baseline|.

    On sparse (binned centroid) spectra both can be 0, so the noise gets a floor at the
    5th percentile of the positive values. Otherwise everything would have infinite SNR.
    """
    window = max(3, int(window) | 1)
    base = ndimage.median_filter(spec, size=window, mode="nearest")
    mad = ndimage.median_filter(np.abs(spec - base), size=window, mode="nearest") * 1.4826
    pos = spec[spec > 0]
    floor = np.percentile(pos, 5) if len(pos) else 1e-12
    return base, np.maximum(mad, max(floor, 1e-12))


def peak_center(axis: np.ndarray, y: np.ndarray, i: int) -> float:
    """Sub-sample peak position from the apex point i and its two neighbours.

    For a Gaussian-ish peak, log(intensity) is a parabola, so fitting a parabola through the
    logs of the 3 top points gives the centre exactly, even with only ~4 points per FWHM.
    (A plain weighted mean over "points above half height" is biased by several ppm when
    the peak is sampled that coarsely, because which edge points get in is lopsided.)
    Falls back to a symmetric weighted mean if a neighbour is <= 0 (sparse binned data).
    """
    if 0 < i < len(axis) - 1 and y[i - 1] > 0 and y[i] > 0 and y[i + 1] > 0:
        a, b, c = np.log(y[i - 1]), np.log(y[i]), np.log(y[i + 1])
        denom = a - 2 * b + c
        if denom < 0:
            delta = 0.5 * (a - c) / denom               # offset in units of points, -0.5..0.5
            delta = float(np.clip(delta, -1, 1))
            step = axis[i + 1] - axis[i] if delta > 0 else axis[i] - axis[i - 1]
            return float(axis[i] + delta * step)
    lo, hi = max(i - 3, 0), min(i + 4, len(axis))
    w = np.clip(y[lo:hi], 0, None)
    return float(np.sum(axis[lo:hi] * w) / w.sum()) if w.sum() > 0 else float(axis[i])


def pick_peaks(axis: np.ndarray, spec: np.ndarray, snr: float, noise_window: int,
               smooth_bins: float = 0.0) -> dict:
    """Find peaks in a mean spectrum. Returns a dict of arrays (mz, height, snr, fwhm, ...)."""
    y = ndimage.gaussian_filter1d(spec, smooth_bins) if smooth_bins > 0 else spec
    base, noise = noise_level(y, noise_window)
    s = (y - base) / noise
    idx, _ = signal.find_peaks(y, height=base + snr * noise, prominence=snr * noise)
    if len(idx) == 0:
        return {k: np.array([]) for k in ("mz", "height", "snr", "fwhm", "fwhm_ppm", "lo", "hi", "index")}
    w, _, lft, rgt = signal.peak_widths(y, idx, rel_height=0.5)
    pos = np.arange(len(axis))
    lo_mz = np.interp(lft, pos, axis)
    hi_mz = np.interp(rgt, pos, axis)
    fwhm = hi_mz - lo_mz
    cmz = np.array([peak_center(axis, y - base, i) for i in idx])
    return {"mz": cmz, "height": y[idx] - base[idx], "snr": s[idx], "fwhm": fwhm,
            "fwhm_ppm": fwhm / cmz * 1e6, "lo": lo_mz, "hi": hi_mz, "index": idx}


def find_ion(axis: np.ndarray, spec: np.ndarray, target: float, search_ppm: float,
             noise_window: int) -> dict | None:
    """Look for one expected ion near `target`; return its measured m/z, SNR, width."""
    lo, hi = target * (1 - search_ppm * 1e-6), target * (1 + search_ppm * 1e-6)
    a, b = np.searchsorted(axis, lo), np.searchsorted(axis, hi)
    if b - a < 1:
        return None
    pad = max(noise_window, 50)
    a0, b0 = max(0, a - pad), min(len(axis), b + pad)
    seg_axis, seg = axis[a0:b0], spec[a0:b0]
    pk = pick_peaks(seg_axis, seg, snr=0.0, noise_window=min(noise_window, len(seg) // 2 * 2 - 1 or 3))
    if len(pk["mz"]) == 0:
        return None
    inside = (pk["mz"] >= lo) & (pk["mz"] <= hi)
    if not inside.any():
        return None
    j = np.flatnonzero(inside)[np.argmax(pk["height"][inside])]
    return {"mz": float(pk["mz"][j]), "ppm_error": float((pk["mz"][j] - target) / target * 1e6),
            "snr": float(pk["snr"][j]), "height": float(pk["height"][j]),
            "fwhm_ppm": float(pk["fwhm_ppm"][j])}


# ---------------------------------------------------------------------------------------
# Aligning peak lists across sections
# ---------------------------------------------------------------------------------------
def align_peaks(tables: list, align_ppm: float) -> list[list[int]]:
    """Group peaks from many sections into consensus features.

    tables: list of 1D arrays of (recalibrated) peak m/z, one per section.
    Returns clusters as lists of (section_index, peak_index) flattened into pairs.
    Greedy single pass over sorted m/z: a peak joins the current cluster if it is within
    align_ppm of the cluster's running median, otherwise it starts a new cluster.
    """
    pool = [(mz, s, k) for s, arr in enumerate(tables) for k, mz in enumerate(arr)]
    pool.sort()
    clusters: list[list] = []
    cur: list = []
    cur_mz: list = []
    for mz, s, k in pool:
        if cur and abs(mz - np.median(cur_mz)) / mz * 1e6 > align_ppm:
            clusters.append(cur)
            cur, cur_mz = [], []
        cur.append((s, k))
        cur_mz.append(mz)
    if cur:
        clusters.append(cur)
    return clusters


def clip_windows(centers: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Neighbouring windows must not overlap: cut both at the midpoint between centres."""
    order = np.argsort(centers)
    c, l, h = centers[order], lo[order].copy(), hi[order].copy()
    mid = (c[1:] + c[:-1]) / 2
    h[:-1] = np.minimum(h[:-1], mid)
    l[1:] = np.maximum(l[1:], mid)
    out_l, out_h = np.empty_like(l), np.empty_like(h)
    out_l[order], out_h[order] = l, h
    return out_l, out_h


def fit_resolution(mz: np.ndarray, fwhm: np.ndarray) -> dict:
    """Fit R = m/FWHM as R ~ a * m^b. b ~ 0 looks like TOF, -0.5 Orbitrap, -1 FT-ICR."""
    ok = (fwhm > 0) & np.isfinite(fwhm)
    if ok.sum() < 3:
        return {"n": int(ok.sum()), "R_median": np.nan, "exponent_b": np.nan}
    R = mz[ok] / fwhm[ok]
    b, loga = np.polyfit(np.log(mz[ok]), np.log(R), 1)
    hint = "TOF-like (R flat with m/z)" if b > -0.25 else ("Orbitrap-like (R ~ m^-0.5)" if b > -0.75
                                                         else "FT-ICR-like (R ~ m^-1)")
    return {"n": int(ok.sum()), "R_median": float(np.median(R)), "exponent_b": float(b),
            "R_at_400": float(np.exp(loga) * 400 ** b), "R_at_800": float(np.exp(loga) * 800 ** b),
            "analyzer_hint": hint}
