"""Small helpers to turn per-pixel values (rows of an AnnData) back into 2D images."""
from __future__ import annotations

import numpy as np
import pandas as pd


def section_box(obs: pd.DataFrame, slide: str, section: str, pad: int = 3):
    """Bounding box (r0, r1, c0, c1) of one section, 0-based rows/cols, padded."""
    sel = (obs["slide"] == slide) & (obs["section"] == section)
    y, x = obs.loc[sel, "y"].to_numpy() - 1, obs.loc[sel, "x"].to_numpy() - 1
    return max(0, y.min() - pad), y.max() + pad + 1, max(0, x.min() - pad), x.max() + pad + 1


def section_image(obs: pd.DataFrame, values: np.ndarray, slide: str, section: str, pad: int = 3,
                  fill=np.nan, box=None):
    """Image of one section (cropped to its box). Returns (img, mask, box, row_index_of_pixels)."""
    sel = ((obs["slide"] == slide) & (obs["section"] == section)).to_numpy()
    r0, r1, c0, c1 = box or section_box(obs, slide, section, pad)
    img = np.full((r1 - r0, c1 - c0), fill, dtype=float)
    mask = np.zeros(img.shape, bool)
    yy = obs["y"].to_numpy()[sel] - 1 - r0
    xx = obs["x"].to_numpy()[sel] - 1 - c0
    img[yy, xx] = values[sel]
    mask[yy, xx] = True
    return img, mask, (r0, r1, c0, c1), np.flatnonzero(sel), (yy, xx)


def scaled01(v: np.ndarray, lo_pct: float = 1, hi_pct: float = 99) -> np.ndarray:
    """Clip to percentiles and rescale to 0..1 (like the tutorial's min01, per lipid)."""
    v = np.asarray(v, float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        return np.zeros_like(v)
    lo, hi = np.percentile(finite, [lo_pct, hi_pct])
    return np.clip((v - lo) / (hi - lo if hi > lo else 1.0), 0, 1)
