"""Brain evidence for our peak names, from the Lipid Brain Atlas reference tables.

Source: "The lipidomic architecture of the mouse brain" (Lipid Brain Atlas, La Manno lab, EPFL), Nature (2026),
doi:10.1038/s41586-026-11050-0, Supplementary Tables 1 and 2 (CC BY 4.0). The compact CSVs
in reference/lamanno2026/ are made by scripts/tools/build_lamanno_reference.py.

What it gives each name (species level, e.g. "HexCer 42:2;O2"):

  points, following the paper's annotation score (Methods, "Lipid annotation")
    msms        8   seen by LC-MS/MS in mouse brain (their tandem MS table or their
                    quantitative whole-brain HILIC-MS/MS)
    lcms        2   seen by LC-MS without MS/MS (their HILIC and reverse-phase LC-MS lists)
    fitzner     1   published brain LC-MS (Fitzner et al.)
    vandenbosch 1   published brain MALDI-MSI (Vandenbosch et al.)
    maldi       1   their own MALDI peak within maldi_match_ppm, annotated with this name
                    and score > 4 (positive mode only, so it can only help lip_pos)
  Their METASPACE term is left out: its FDR belongs to their data, not ours. LIPID MAPS
  alone scores 0, as in the paper.

  tier
    high            score > 4 (in practice: MS/MS, or both LC-MS lists plus a publication)
    medium          score 2 to 4
    low             score 1
    database only   score 0, but their tables do cover this lipid class
    not covered     no table covers the class (gangliosides, LSM, metabolites, ...), so
                    their data can't say anything either way
    "" (empty)      no name

  Ambiguous peaks: only candidates in the same ionization regime are compared (if default
  adducts such as [M-H]- compete with class-specific ones such as [M-CH3]- of a choline
  lipid, only the default ones count: abundance says nothing about how well PC ionizes in
  negative mode). If every remaining candidate belongs to a class their quantitative whole-brain
  LC-MS/MS measures, a candidate that holds more than molar_fraction (0.8) of the
  candidates' summed nmol% is picked, as in the paper. A candidate of a measured class that
  wasn't detected counts as 0. Otherwise nothing is decided.

BMP and PG have the same formula; our name "BMP/PG 36:1" is looked up as both "BMP 36:1"
and "PG 36:1" (evidence then says the formula exists in brain, not which of the two it is).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

POINTS = {"msms": 8, "lcms": 2, "fitzner": 1, "vandenbosch": 1, "maldi": 1}
SOURCES = list(POINTS)

SPHINGO = {"Cer", "SM", "HexCer", "Hex2Cer", "Hex3Cer", "SHexCer", "LSM", "SHexSph", "HexSph",
           "GM1", "GM2", "GM3", "GD1", "GD2", "GD3", "GT1", "CerP", "PE-Cer", "PI-Cer", "IPC"}
ALIASES = {"Hex1Cer": "HexCer", "GlcCer": "HexCer", "GalCer": "HexCer", "LacCer": "Hex2Cer",
           "LBPA": "BMP", "DHCer": "Cer", "DHSM": "SM", "TAG": "TG", "SE 27:1": "CE"}
OX_LETTER = {"m": 1, "d": 2, "t": 3}


def _clean(s: str) -> str:
    s = str(s).strip()
    s = s.split("&/or")[0].strip()
    s = re.sub(r"\s*\[.*$", "", s)                      # drop an adduct: "PC 34:1 [M+H]+"
    s = re.sub(r"([OP])-\s+", r"\1-", s)                 # "PC O- 36:2"
    s = re.sub(r"(\d+)\.(\d+)$", r"\1:\2", s)            # "PC 34.1"
    return s.strip()


def canonical(name) -> str | None:
    """Any of the naming styles in the reference tables -> 'CLASS [O-]C:DB[;On]'.

    Handles 'PE 18:1_22:6', 'Cer d18:1/24:1', 'PE(O-16:1/20:4)', 'PE(P-16:0/20:4)' (vinyl ether
    = O- with one more double bond), 'HexCer(t42:2)', 'SHexCer 36:1;3', 'Hex1Cer d18:1/24:1',
    'LacCer(d30:0)', 'LBPA(34:1)', 'BMP/PG 36:1', 'SM t38:5'. Oxidized ('(OH)') and unparsable
    names give None.
    """
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return None
    s = _clean(name)
    if not s or "(OH)" in s or "OH)" in s:
        return None
    m = re.match(r"^([A-Za-z][\w\-/]*?)\((.+)\)$", s) or re.match(r"^([A-Za-z][\w\-/]*)\s+(.+)$", s)
    if not m:
        return None
    cls, body = m.group(1).strip(), m.group(2).strip().replace(" ", "")
    if cls == "SE" and body.startswith("27:1"):
        cls, body = "CE", body.split("/", 1)[-1] if "/" in body else body
    cls = ALIASES.get(cls, cls)
    ox = None
    mo = re.search(r";(O?)(\d*)$", body)
    if mo:
        ox = int(mo.group(2)) if mo.group(2) else 1
        body = body[:mo.start()]
    ether, plasm = False, False
    if body[:2] in ("O-", "P-"):
        ether, plasm = True, body[0] == "P"
        body = body[2:]
    if body[:1] in OX_LETTER and re.match(r"^[mdt]\d", body):
        ox = OX_LETTER[body[0]]
        body = body[1:]
    chains = re.split(r"[_/]", body)
    try:
        parts = [tuple(int(x) for x in re.match(r"^(\d+):(\d+)", c).groups()) for c in chains if c]
    except AttributeError:
        return None
    if not parts:
        return None
    # 'CL 72:8_20:4': a total plus one identified chain; the first number is already the sum
    if cls in ("CL", "TG") and len(parts) > 1 and parts[0][0] >= 40:
        parts = parts[:1]
    c, db = sum(p[0] for p in parts), sum(p[1] for p in parts)
    if plasm:
        db += 1
    if cls in SPHINGO and ox is None:
        ox = 2
    if ox is not None and cls not in SPHINGO:
        return None                                       # oxidized glycerophospholipid ('PC 38:6;O')
    out = f"{cls} {'O-' if ether else ''}{c}:{db}"
    if ox is not None and cls in SPHINGO:
        out += f";O{ox}"
    return out


def class_of(canon: str | None) -> str:
    if not canon:
        return ""
    return canon.split(" ", 1)[0]


def lookup_names(name) -> list[str]:
    """Canonical names to look up for one of OUR names (BMP/PG -> both)."""
    c = canonical(name)
    if not c:
        return []
    cls, rest = c.split(" ", 1)
    if cls == "BMP/PG":
        return [f"BMP {rest}", f"PG {rest}"]
    return [c]


def split_multi(annotation: str) -> list[str]:
    """Their MALDI 'Annotation' can hold several names: 'PC 34:1 PC 36:4 PE 37:1'."""
    s = re.sub(r"_db$", "", str(annotation).strip())
    toks = re.findall(r"[A-Za-z][\w\-/]*\s+(?:O-)?\d+:\d+(?:;O\d?)?", s)
    return toks or ([s] if s else [])


@dataclass
class Reference:
    species: dict = field(default_factory=dict)       # source -> set of canonical names
    covered: set = field(default_factory=set)         # classes any table covers
    quant: dict = field(default_factory=dict)         # canonical name -> mean nmol%
    quant_classes: set = field(default_factory=set)   # classes the quantitative assay measures
    maldi: pd.DataFrame | None = None                 # mz, names (list), score

    def score(self, name, mz: float | None = None, maldi_ppm: float = 5.0, maldi_min_score: float = 4.0,
              polarity: str = "positive") -> dict:
        names = lookup_names(name)
        hit = {s: False for s in SOURCES}
        for n in names:
            for s in ("msms", "lcms", "fitzner", "vandenbosch"):
                hit[s] |= n in self.species.get(s, set())
        if self.maldi is not None and mz is not None and np.isfinite(mz) and polarity == "positive" and names:
            M = self.maldi
            near = M[(np.abs(M["mz"] - mz) / mz * 1e6 <= maldi_ppm) & (M["score"] > maldi_min_score)]
            hit["maldi"] = any(any(n in canon for n in names) for canon in near["names"])
        score = sum(POINTS[s] for s, v in hit.items() if v)
        cls = class_of(names[0]) if names else ""
        covered = any(class_of(n) in self.covered for n in names)
        if not names:
            # a name that isn't a lipid sum composition (metabolites, LIPID MAPS natural products)
            tier = "not covered" if str(name).strip() else ""
        elif score > 4:
            tier = "high"
        elif score >= 2:
            tier = "medium"
        elif score >= 1:
            tier = "low"
        elif covered:
            tier = "database only"
        else:
            tier = "not covered"
        return {"evidence_score": score, "confidence_tier": tier,
                "evidence_sources": ",".join(s for s, v in hit.items() if v),
                "evidence_note": ("BMP/PG looked up as both; evidence is for the formula" if cls in ("BMP", "PG")
                                  and str(name).startswith("BMP/PG") else "")}

    def resolve(self, candidate_names: list[str], fraction: float = 0.8) -> str | None:
        """The candidate holding > fraction of the summed whole-brain nmol%, or None."""
        cands = list(dict.fromkeys(n for n in candidate_names if n))
        if len(cands) < 2:
            return None
        ab = []
        for n in cands:
            ln = lookup_names(n)
            if not ln or any(class_of(x) not in self.quant_classes for x in ln):
                return None                       # a class the assay doesn't measure: can't decide
            ab.append(sum(self.quant.get(x, 0.0) for x in ln))
        ab = np.asarray(ab, float)
        if ab.sum() <= 0:
            return None
        share = ab / ab.sum()
        j = int(np.argmax(share))
        return cands[j] if share[j] > fraction else None


def load(ref_dir: Path) -> Reference | None:
    """Read reference/lamanno2026/*.csv; None if the folder isn't there."""
    ref_dir = Path(ref_dir)
    sp = ref_dir / "species_sources.csv"
    if not sp.exists():
        return None
    S = pd.read_csv(sp)
    R = Reference()
    for src, g in S.groupby("source"):
        R.species[src] = set(g["species"])
    R.covered = set(S["lipid_class"])
    q = ref_dir / "quant_whole_brain.csv"
    if q.exists():
        Q = pd.read_csv(q)
        R.quant = dict(zip(Q["species"], Q["mean_nmol_pct"]))
        R.quant_classes = set(Q["lipid_class"])
    mp = ref_dir / "maldi_peaks.csv"
    if mp.exists():
        M = pd.read_csv(mp)
        M["names"] = [[c for c in (canonical(t) for t in split_multi(a)) if c] for a in M["annotation"].fillna("")]
        R.maldi = M[["mz", "names", "score"]]
    return R
