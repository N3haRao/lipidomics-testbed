"""Turning peak m/z values into sum-composition names (SM 36:1;O2, HexCer 42:2;O2, ...).

The library is generated from reference/lipid_classes.yaml (anchor formula + CH2 / H2
steps), so every candidate has an exact formula. Optionally we keep only compositions whose
formula exists in LIPID MAPS (cached by 00_setup), which removes chemically odd ones.
For the metabolite mode we add HMDB (or the small curated list) on top.

A match = |ppm error| <= tolerance for some (species, adduct). One peak can match several
species: that's reported, never hidden. Flags we add:
  ambiguous      more than one composition fits within tolerance
  isotope        looks like the 13C (M+1) peak of a lower-mass peak in our list
  typeII         a peak one double bond up (-2.016 Da) is strong enough that ITS M+2 isotope
                 could sit inside this peak (they are only ~9 mDa apart)
  matrix         matches a DAN matrix ion
"""
from __future__ import annotations


import numpy as np
import pandas as pd
import yaml

from . import chem
from .config import REPO_ROOT

REF = REPO_ROOT / "reference"
SPHINGO = {"SM", "Cer", "HexCer", "Hex2Cer", "SHexCer", "GM3", "GM2", "GM1", "GD1", "LSM"}


def composition_name(cls: str, c: int, db: int, ox: str | None, ether: bool) -> str:
    e = "O-" if ether else ""
    o = f";{ox}" if ox else ""
    return f"{cls} {e}{c}:{db}{o}"


def build_lipid_library(lmsd: pd.DataFrame | None = None, log=None) -> pd.DataFrame:
    spec = yaml.safe_load((REF / "lipid_classes.yaml").read_text())
    rows = []
    for cls, d in spec["classes"].items():
        oxy = d.get("oxygen") or {None: 0}
        ethers = [False, True] if d.get("ether") else [False]
        for c in range(d["c_range"][0], d["c_range"][1] + 1):
            for db in range(d["db_range"][0], d["db_range"][1] + 1):
                if db > c // 2:     # chemically silly
                    continue
                for ox, d_o in oxy.items():
                    for eth in ethers:
                        f = chem.shift_formula(d["anchor"], c - d["anchor_c"], db - d["anchor_db"], d_o)
                        if eth:   # ether: one ester carbonyl O becomes CH2 -> formula -O +2H
                            cc = chem.parse_formula(f)
                            cc["O"] -= 1
                            cc["H"] += 2
                            f = chem.format_formula(cc)
                        rows.append({"name": composition_name(cls, c, db, ox, eth), "lipid_class": cls,
                                     "tag": d.get("tag"), "total_c": c, "total_db": db, "oxygen": ox,
                                     "ether": eth, "formula": f, "source": "class_rules"})
    for fx in spec.get("fixed", []):
        rows.append({"name": fx["name"], "lipid_class": fx["name"].split()[0], "tag": fx.get("tag"),
                     "total_c": None, "total_db": None, "oxygen": None, "ether": False,
                     "formula": fx["formula"], "source": "fixed"})
    lib = pd.DataFrame(rows)
    lib["neutral_mass"] = [chem.mono_mass(f) for f in lib["formula"]]
    lib["in_lmsd"] = np.nan
    if lmsd is not None and len(lmsd):
        known = set(lmsd["formula"])
        lib["in_lmsd"] = lib["formula"].isin(known)
        if log:
            log.info(f"lipid library: {len(lib)} generated compositions, {int(lib.in_lmsd.sum())} have a "
                     "formula present in LIPID MAPS")
    return lib


def build_metabolite_library(ctx, log) -> pd.DataFrame:
    which = ctx.get("annotation.metabolite_db", "curated")
    cur = pd.read_csv(REF / "metabolites_curated.csv").fillna({"tag": ""})
    cur = cur.rename(columns={"class": "lipid_class"})
    cur["source"] = "curated"
    if which == "hmdb":
        from .refdb import cached_table
        h = cached_table(ctx, "hmdb")
        if h is None:
            log.warning("annotation.metabolite_db = hmdb but no HMDB cache found "
                        "(run 00_setup.py --download). Falling back to the curated list.")
        else:
            h = h.rename(columns={"db_class": "lipid_class"})
            h["tag"] = ""
            h["source"] = "hmdb"
            # Curated entries keep their tags (GPD etc), HMDB adds breadth
            cur = pd.concat([cur, h[["name", "formula", "lipid_class", "tag", "source"]]], ignore_index=True)
    cur["neutral_mass"] = [chem.mono_mass(f) for f in cur["formula"]]
    for col in ("total_c", "total_db", "oxygen"):
        cur[col] = None
    cur["ether"] = False
    cur["in_lmsd"] = np.nan
    return cur


def ion_table(lib: pd.DataFrame, polarity: str, adducts_cfg: dict) -> pd.DataFrame:
    """Every (species, adduct) pair with its m/z. Class-specific adducts only for their class."""
    general = adducts_cfg[polarity]
    special = (adducts_cfg.get("class_specific") or {}).get(polarity) or {}
    excluded = (adducts_cfg.get("excluded") or {}).get(polarity) or {}
    parts = []
    for a in general:
        t = lib.copy()
        drop = [cls for cls, adds in excluded.items() if a in adds]
        t = t[~t["lipid_class"].isin(drop)]
        t["adduct"] = a
        t["adduct_kind"] = "default"
        parts.append(t)
    for cls, adds in special.items():
        sub = lib[lib["lipid_class"] == cls]
        for a in adds:
            t = sub.copy()
            t["adduct"] = a
            t["adduct_kind"] = "class_specific"
            parts.append(t)
    ions = pd.concat(parts, ignore_index=True)
    # Note: ion_mz handles the n*M, the pieces and the charge, see chem.py
    ions["ion_mz"] = [chem.ion_mz(f, a) for f, a in zip(ions["formula"], ions["adduct"])]
    return ions.sort_values("ion_mz").reset_index(drop=True)


def match_peaks(peaks: pd.DataFrame, ions: pd.DataFrame, tol_ppm: float, use_lmsd_filter: bool,
                max_report: int) -> pd.DataFrame:
    """For each peak (column 'mz'), list candidates within tol_ppm and pick a best one."""
    ions_mz = ions["ion_mz"].to_numpy()
    out = []
    for mz in peaks["mz"].to_numpy():
        lo, hi = mz * (1 - tol_ppm * 1e-6), mz * (1 + tol_ppm * 1e-6)
        a, b = np.searchsorted(ions_mz, lo), np.searchsorted(ions_mz, hi, side="right")
        cand = ions.iloc[a:b].copy()
        if use_lmsd_filter and "in_lmsd" in cand and cand["in_lmsd"].notna().any():
            # class-rule compositions not in LMSD are dropped; fixed/curated/HMDB entries stay
            cand = cand[(cand["in_lmsd"] != False) | (cand["source"] != "class_rules")]  # noqa: E712
        cand["ppm_error"] = (mz - cand["ion_mz"]) / cand["ion_mz"] * 1e6
        # Ranking: even chains first for glycerolipids/sphingolipids (odd chains are rare in
        # mammalian brain), then default adducts over class-specific ones, then smallest |ppm|.
        odd = cand["total_c"].fillna(0).astype(int) % 2 == 1
        cand["rank_key"] = ((cand["adduct_kind"] != "default").astype(int) * 5
                            + odd.astype(int) * 10 + cand["ppm_error"].abs() / max(tol_ppm, 1e-9))
        cand = cand.sort_values("rank_key")
        if len(cand):
            best = cand.iloc[0].copy()
            if not (isinstance(best["tag"], str) and best["tag"]):
                # e.g. an HMDB entry won the ranking but a curated entry with the SAME formula has a tag
                same = cand[(cand["formula"] == best["formula"]) & cand["tag"].apply(lambda t: isinstance(t, str) and t != "")]
                if len(same):
                    best["tag"] = same["tag"].iloc[0]
            names = list(dict.fromkeys(cand["name"] + " " + cand["adduct"]))
            comps = cand["name"].unique()
            tags = sorted({t for t in cand["tag"] if isinstance(t, str) and t})
            out.append({
                "annotation": best["name"], "adduct": best["adduct"], "formula": best["formula"],
                "lipid_class": best["lipid_class"], "tag": best["tag"] if isinstance(best["tag"], str) else "",
                "ppm_error": float(best["ppm_error"]), "n_candidates": len(comps),
                "ambiguous": len(comps) > 1, "candidate_tags": ";".join(tags),
                "candidates": " | ".join(names[:max_report]) + (" | ..." if len(names) > max_report else ""),
                "total_c": best["total_c"], "total_db": best["total_db"],
            })
        else:
            out.append({"annotation": "", "adduct": "", "formula": "", "lipid_class": "", "tag": "",
                        "ppm_error": np.nan, "n_candidates": 0, "ambiguous": False, "candidate_tags": "",
                        "candidates": "", "total_c": None, "total_db": None})
    return pd.DataFrame(out, index=peaks.index)


def isotope_flags(mz: np.ndarray, intensity: np.ndarray, tol_ppm: float) -> pd.DataFrame:
    """Flag M+1 isotopes and type II (M+2 of one-more-double-bond) overlaps among our own peaks."""
    order = np.argsort(mz)
    smz, sint = mz[order], intensity[order]
    is_iso = np.zeros(len(mz), bool)
    iso_of = np.full(len(mz), np.nan)
    t2_frac = np.zeros(len(mz))

    def nearest(target, extra_da=0.0):
        j = np.searchsorted(smz, target)
        best = None
        for k in (j - 1, j):
            if 0 <= k < len(smz) and abs(smz[k] - target) <= target * tol_ppm * 1e-6 + extra_da:
                if best is None or abs(smz[k] - target) < abs(smz[best] - target):
                    best = k
        return best

    for i in range(len(smz)):
        n_c = max(1, int(round(smz[i] / 14.5)))    # rough carbon count for a lipid-like ion
        k = nearest(smz[i] - chem.C13_SHIFT)
        if k is not None and sint[i] <= sint[k] * max(1.0, 1.5 * chem.m1_ratio(n_c)):
            is_iso[order[i]] = True
            iso_of[order[i]] = smz[k]
        # The M+2 of the one-more-double-bond species sits ~8.9 mDa BELOW our peak, so we look
        # for that species with a tolerance widened by that shift.
        k2 = nearest(smz[i] - chem.H2, extra_da=chem.H2 - 2 * chem.C13_SHIFT)
        if k2 is not None and sint[i] > 0:
            t2_frac[order[i]] = chem.m2_overlap_fraction(n_c) * sint[k2] / sint[i]
    t2_frac[is_iso] = 0.0     # an isotope peak's own "type II" partner is just its M+0
    return pd.DataFrame({"is_isotope": is_iso, "isotope_of_mz": iso_of, "typeII_fraction": t2_frac,
                         "typeII_flag": t2_frac > 0.10})


def hypothesis_direction(tag: str, hyp_cfg: dict) -> str:
    """'accumulated', 'depleted', 'both', 'secondary' or '' for a class tag."""
    if not tag:
        return ""
    acc = tag in (hyp_cfg.get("accumulated") or [])
    dep = tag in (hyp_cfg.get("depleted") or [])
    if acc and dep:
        return "both"
    if acc:
        return "accumulated"
    if dep:
        return "depleted"
    if tag in (hyp_cfg.get("secondary") or []):
        return "secondary"
    return ""
