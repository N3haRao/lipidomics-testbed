"""Feature QC before testing, following the Lipid Brain Atlas (Nature 2026) Methods.

Their rules, and what we do with each:
  [M+0] isotopologues only      isotope peaks are never tested here (05 flags them)
  biological: tissue > matrix   off/on-tissue signal ratio from 04 below max_off_on_ratio
  spatially coherent            Moran's I from 08 (explore.morans_i: I >= min_I in at least
                                min_section_fraction of the sections)
  annotated in a database       any name from 05 (class rules, curated, HMDB, second pass)
  one adduct per lipid          below

None of these looks at genotype, so the filter can't favour a result. The peaks are not
removed from the matrix (TIC normalization stays the same); the QC only decides which peaks
form the untargeted FDR family (differential.untargeted_peaks = lamanno) and which species
enter the feature-QC sensitivity version of the class tests.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def adduct_choice(var: pd.DataFrame, X, cols: np.ndarray, min_r: float = 0.5, max_pixels: int = 50000,
                  rng=None) -> pd.DataFrame:
    """One ion per lipid name.

    var: peaks table (index = peak ids, columns annotation, is_isotope); X: pixels x all peaks
    (array or h5py-backed), cols: positions in X of var's rows. Within a group of non-isotope
    peaks with the same name, the ion with the smallest fraction of zero pixels wins (ties: the
    larger total signal). Another ion is 'redundant' if its log intensities correlate with the
    winner's at r >= min_r, otherwise 'kept, not correlated' (probably a different compound).

    Returns a table indexed like var: adduct_group (the name), adduct_n, adduct_status
    ('single', 'chosen', 'redundant', 'kept, not correlated', '' for unnamed), adduct_r,
    adduct_keep (bool).
    """
    rng = rng or np.random.default_rng(0)
    name = var.get("annotation", pd.Series("", index=var.index)).fillna("").astype(str)
    iso = var.get("is_isotope", pd.Series(False, index=var.index)).fillna(False).astype(bool)
    out = pd.DataFrame({"adduct_group": np.where(name != "", name, ""), "adduct_n": 0,
                        "adduct_status": "", "adduct_r": np.nan, "adduct_keep": True}, index=var.index)
    named = (name != "") & ~iso
    groups = name[named].groupby(name[named]).groups
    n_pix = X.shape[0]
    rows = np.sort(rng.choice(n_pix, size=min(n_pix, max_pixels), replace=False))
    pos = dict(zip(var.index, cols))
    for nm, ids in groups.items():
        ids = list(ids)
        out.loc[ids, "adduct_n"] = len(ids)
        if len(ids) == 1:
            out.loc[ids, "adduct_status"] = "single"
            continue
        M = np.asarray(X[:, [pos[i] for i in ids]][rows], dtype=np.float64)
        zero = (M <= 0).mean(axis=0)
        total = M.sum(axis=0)
        L = np.log1p(M)
        ok = L.std(axis=0) > 0
        C = np.full((len(ids), len(ids)), np.nan)
        if ok.sum() > 1:
            C[np.ix_(ok, ok)] = np.corrcoef(L[:, ok].T)
        linked = np.nan_to_num(C, nan=0.0) >= min_r
        np.fill_diagonal(linked, True)
        # The lipid's ions are the largest set that correlate with one another: start from the
        # ion with the most correlated partners (ties: fewest zeros, then most signal), so an
        # unrelated ion with few zeros can't become the reference. Within that set the ion with
        # the fewest zero pixels (then the most signal) is kept.
        hub = np.lexsort((-total, zero, -linked.sum(axis=1)))[0]
        members = np.flatnonzero(linked[hub])
        best = members[np.lexsort((-total[members], zero[members]))[0]]
        for j, pid in enumerate(ids):
            r = 1.0 if j == best else float(C[best, j]) if np.isfinite(C[best, j]) else np.nan
            out.loc[pid, "adduct_r"] = r
            if j == best:
                out.loc[pid, "adduct_status"] = "chosen"
            elif j in members:
                out.loc[pid, ["adduct_status", "adduct_keep"]] = ["redundant", False]
            else:
                out.loc[pid, "adduct_status"] = "kept, not correlated"
    return out


def table(var: pd.DataFrame, structured: pd.Series, adducts: pd.DataFrame, max_off_on: float) -> pd.DataFrame:
    """Per-peak QC flags and the combined pass/fail (qc_pass)."""
    iso = var.get("is_isotope", pd.Series(False, index=var.index)).fillna(False).astype(bool)
    named = var.get("annotation", pd.Series("", index=var.index)).fillna("").astype(str) != ""
    ratio = var.get("off_on_ratio", pd.Series(np.nan, index=var.index)).astype(float)
    T = pd.DataFrame({
        "mz": var.get("mz"), "label": var.get("label"),
        "confidence_tier": var.get("confidence_tier", pd.Series("", index=var.index)),
        "monoisotopic": ~iso,
        "annotated": named,
        "spatially_structured": structured.reindex(var.index).fillna(False).astype(bool),
        "off_on_ratio": ratio,
        # no ratio measured (no off-tissue pixels) -> don't hold it against the peak
        "tissue_over_matrix": ratio.isna() | (ratio < max_off_on),
    }, index=var.index)
    T = T.join(adducts)
    T["qc_pass"] = (T["monoisotopic"] & T["annotated"] & T["spatially_structured"] & T["tissue_over_matrix"]
                    & T["adduct_keep"].fillna(True).astype(bool))
    return T
