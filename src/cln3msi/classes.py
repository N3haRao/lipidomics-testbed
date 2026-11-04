"""Class-level and ratio tests, and the slide consistency check (both run inside phase 09).

Class-level tests
-----------------
Single lipids are noisy and there are hundreds of them, so the FDR correction is strict.
Here we add up the species of one class (all sulfatides, all gangliosides, ...) per animal
and region, and test those totals with the same paired test as single peaks. Ratios
(sulfatide / HexCer, GPI / LPI, ...) are taken within the same animal and region, so the
TIC normalization, the matrix coating and most of a left/right trend cancel out.

The classes and ratios are set in config.yaml (class_tests), before looking at results.
All class and ratio tests of one mode form one BH family.

Slide consistency check
-----------------------
For every slide, log2(Mutant / Control) of every tested peak is that slide's "difference
profile". If the labels are right, the slides share the genotype effect, so their profiles
correlate positively (or at least not negatively). A slide whose profile is the mirror
image of the others (clearly negative correlation, the same in several regions) looks like
its two sections were swapped: wrong sample sheet, or the slide went into the instrument
the other way round. It's a QC warning. Never relabel a slide because of this; check the
sectioning and mounting records first, then use qc.sections.flip_slides if they confirm it.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import stats

CHAIN = re.compile(r"^(.+?)\s+(?:[OP]-)?\d+:\d+")


def candidate_class(name: str) -> str:
    """'LPC O-16:2 [M+H]+' -> 'LPC', 'BMP/PG 36:1 [M-H]-' -> 'BMP/PG', 'Taurine [M-H]-' -> 'Taurine'."""
    name = re.sub(r"\s*\[.*$", "", str(name)).strip()
    m = CHAIN.match(name)
    return m.group(1) if m else name


def usable_peaks(ann: pd.DataFrame, min_sections: int, ambiguity: str = "tags") -> pd.DataFrame:
    """Annotated, not an isotope, seen in enough sections, and not ambiguous between classes.

    ambiguity = "tags":    drop a peak only if its candidates carry more than one hypothesis tag
                           (e.g. SHexCer and BMP_PG). PC/PE-type overlaps are kept.
               "classes":  drop a peak if its candidates belong to more than one lipid class.
               "none":     keep everything.
    """
    a = ann.copy()
    a = a[a["annotation"].fillna("").astype(str) != ""]
    if "is_isotope" in a:
        a = a[~a["is_isotope"].fillna(False).astype(bool)]
    if "n_sections" in a and min_sections:
        a = a[a["n_sections"].fillna(0) >= min_sections]
    if ambiguity == "tags" and "candidate_tags" in a:
        a = a[a["candidate_tags"].fillna("").astype(str).map(lambda t: len({x for x in t.split(";") if x}) <= 1)]
    elif ambiguity == "classes" and "candidates" in a:
        def one_class(c):
            parts = [p for p in str(c).split(" | ") if p and p != "..."]
            return len({candidate_class(p) for p in parts}) <= 1
        a = a[a["candidates"].fillna("").map(one_class)]
    return a


def members(ann: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """Peaks of one class entry. spec keys: classes (lipid_class list), names (exact annotation
    list), suffix (annotation must end with it, e.g. ';O3')."""
    sel = pd.Series(False, index=ann.index)
    if spec.get("classes"):
        sel |= ann["lipid_class"].isin(spec["classes"])
    if spec.get("names"):
        sel |= ann["annotation"].isin(spec["names"])
    if spec.get("suffix"):
        sel &= ann["annotation"].astype(str).str.endswith(spec["suffix"])
    return ann[sel]


def totals(means_long: pd.DataFrame, ann: pd.DataFrame, class_cfg: dict, ratio_cfg: list,
           units: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per unit: animals x features table of class totals and ratios.

    Returns (values_long with columns unit, kind, feature, n_species, animal_id, value,
             members table with columns feature, n_species, species)."""
    peak_col = "peak_id" if "peak_id" in ann else ann.index.name or "index"
    a = ann.reset_index() if peak_col not in ann else ann
    rows, mem = [], []
    for unit in units:
        sub = means_long[means_long["unit"] == unit]
        if sub.empty:
            continue
        wide = sub.pivot_table(index="animal_id", columns="peak", values="mean_norm_intensity", aggfunc="sum")
        tot, nsp = {}, {}
        for name, spec in class_cfg.items():
            m = members(a, spec or {})
            ids = [p for p in m[peak_col] if p in wide.columns]
            if not ids:
                continue
            # an animal with too few pixels has NaN for every peak: keep it NaN, don't sum to 0
            tot[name] = wide[ids].sum(axis=1, min_count=1)
            nsp[name] = len(ids)
            if unit == units[0]:
                lab = m.set_index(peak_col).loc[ids]
                mem.append({"feature": name, "n_species": len(ids),
                            "species": "; ".join((lab["annotation"].astype(str) + " " +
                                                  lab.get("adduct", pd.Series("", index=lab.index)).fillna("").astype(str)).str.strip())})
        for name, v in tot.items():
            rows.append(pd.DataFrame({"unit": unit, "kind": "class", "feature": name, "n_species": nsp[name],
                                      "animal_id": v.index, "value": v.to_numpy()}))
        for num, den in ratio_cfg or []:
            if num in tot and den in tot:
                v = tot[num] / tot[den].replace(0, np.nan)
                rows.append(pd.DataFrame({"unit": unit, "kind": "ratio", "feature": f"{num} / {den}",
                                          "n_species": np.nan, "animal_id": v.index, "value": v.to_numpy()}))
    vals = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["unit", "kind", "feature", "n_species", "animal_id", "value"])
    return vals, pd.DataFrame(mem)


def per_pair(values: pd.Series, animals: pd.DataFrame, control: str, mutant: str) -> pd.Series:
    """values: animal -> value (log2 scale). Returns slide -> Mutant minus Control."""
    out = {}
    for s, g in animals.groupby("slide"):
        c, m = g.index[g.genotype == control], g.index[g.genotype == mutant]
        if len(c) == 1 and len(m) == 1 and c[0] in values.index and m[0] in values.index:
            out[s] = values[m[0]] - values[c[0]]
    return pd.Series(out, dtype=float)


def test_features(vals: pd.DataFrame, animals: pd.DataFrame, control: str, mutant: str,
                  min_pairs: int) -> pd.DataFrame:
    """Paired tests for every (unit, feature), plus the per-slide log2FC as pair_<slide> columns."""
    out = []
    for unit, g in vals.groupby("unit", sort=False):
        wide = g.pivot_table(index="animal_id", columns="feature", values="value")
        wide = wide.reindex(index=[i for i in wide.index if i in animals.index])
        log2 = np.log2(wide + stats.pseudocount(wide))
        res = stats.paired_tests(log2, animals, control, mutant, min_pairs=min_pairs)
        res = res.rename(columns={"peak": "feature"})
        res.insert(0, "unit", unit)
        meta = g.drop_duplicates("feature").set_index("feature")[["kind", "n_species"]]
        res = res.join(meta, on="feature")
        pp = pd.DataFrame({f: per_pair(log2[f].dropna(), animals, control, mutant) for f in log2.columns}).T
        pp.columns = [f"pair_{s}" for s in pp.columns]
        res = res.join(pp, on="feature")
        out.append(res)
    if not out:
        return pd.DataFrame()
    R = pd.concat(out, ignore_index=True)
    R["q_bh"] = stats.bh(R["p_ttest"])
    R["pct_change"] = 100 * (2 ** R["log2fc"] - 1)
    return R


# ---- slide consistency ------------------------------------------------------------------

def pair_profiles(means_long: pd.DataFrame, animals: pd.DataFrame, control: str, mutant: str,
                  unit: str) -> pd.DataFrame:
    """peaks x slides table of log2(Mutant/Control); peaks with a zero or missing mean are left out."""
    sub = means_long[means_long["unit"] == unit]
    w = sub.pivot_table(index="peak", columns="animal_id", values="mean_norm_intensity")
    w = w[(w > 0).all(axis=1)]
    lw = np.log2(w)
    return pd.DataFrame({s: per_pair_cols(lw, g, control, mutant) for s, g in animals.groupby("slide")}).dropna(axis=1, how="all")


def per_pair_cols(lw: pd.DataFrame, g: pd.DataFrame, control: str, mutant: str) -> pd.Series:
    c, m = g.index[g.genotype == control], g.index[g.genotype == mutant]
    if len(c) == 1 and len(m) == 1 and c[0] in lw and m[0] in lw:
        return lw[m[0]] - lw[c[0]]
    return pd.Series(np.nan, index=lw.index)


def consistency(means_long: pd.DataFrame, animals: pd.DataFrame, control: str, mutant: str,
                units: list[str], flag_below: float, min_peaks: int = 20) -> tuple[dict, pd.DataFrame, list[str]]:
    """Spearman correlation of the per-slide difference profiles, per unit.

    Returns (unit -> correlation matrix, table slide x unit of the median correlation with the
    other slides, list of flagged slides). A slide is flagged when its median correlation is
    below flag_below in at least half of the units checked."""
    mats, med = {}, {}
    for u in units:
        P = pair_profiles(means_long, animals, control, mutant, u)
        if P.shape[1] < 3 or len(P) < min_peaks:
            continue
        C = P.rank().corr()             # Spearman = Pearson on ranks
        mats[u] = C
        med[u] = pd.Series({s: float(np.median(C.loc[s].drop(s))) for s in C.index})
    if not med:
        return {}, pd.DataFrame(), []
    T = pd.DataFrame(med)
    T.index.name = "slide"
    below = (T < flag_below).sum(axis=1)
    T["units_below_threshold"] = below
    T["flagged"] = below >= np.ceil(len(med) / 2)
    return mats, T, list(T.index[T["flagged"]])


def sensitivity(vals: pd.DataFrame, animals: pd.DataFrame, control: str, mutant: str,
                flagged: list[str], min_pairs: int) -> pd.DataFrame:
    """Redo the class/ratio tests with the flagged slides (a) left out and (b) flipped.
    Exploratory only: picking the analysis after seeing the data is circular."""
    if not flagged:
        return pd.DataFrame()
    keys = ["unit", "feature"]
    cols = ["log2fc", "p_ttest", "p_wilcoxon", "n_pairs", "n_pairs_up", "q_bh"]
    keep = animals[~animals["slide"].isin(flagged)]
    drop = test_features(vals[vals.animal_id.isin(keep.index)], keep, control, mutant, min_pairs=2)
    flip = animals.copy()
    sel = flip["slide"].isin(flagged)
    flip.loc[sel, "genotype"] = flip.loc[sel, "genotype"].map({control: mutant, mutant: control})
    flp = test_features(vals, flip, control, mutant, min_pairs=min_pairs)
    tag = "+".join(flagged)
    out = drop[keys + cols].rename(columns={c: f"without_{tag}_{c}" for c in cols})
    return out.merge(flp[keys + cols].rename(columns={c: f"flipped_{tag}_{c}" for c in cols}), on=keys, how="outer")
