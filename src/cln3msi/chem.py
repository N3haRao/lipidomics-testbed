"""Masses from formulas, so no m/z value in this pipeline is ever typed in by hand.

    mono_mass("C42H82NO8P")            -> 759.5778  (PC 34:1, neutral)
    ion_mz("C42H82NO8P", "[M+K]+")     -> 798.5410

Adducts are written the usual way: [M+H]+, [M-H]-, [M+Na]+, [M-2H]2-, [M+H-H2O]+, [2M-H]-.
The parser reads the multiplier of M, the list of +/- formula pieces, and the charge.
"""
from __future__ import annotations

import re
from collections import Counter

# Monoisotopic masses of the most abundant isotope (IUPAC / NIST values)
ELEMENTS = {
    "H": 1.00782503207, "C": 12.0, "N": 14.0030740048, "O": 15.99491461956,
    "P": 30.97376163, "S": 31.97207100, "Na": 22.9897692809, "K": 38.96370668,
    "Cl": 34.96885268, "F": 18.99840322, "Li": 7.01600455,
}
ELECTRON = 0.00054857990946
PROTON = ELEMENTS["H"] - ELECTRON                  # 1.007276
C13_SHIFT = 1.0033548378                           # 13C minus 12C, used for isotope checks
H2 = 2 * ELEMENTS["H"]                             # one double bond = minus H2

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Counter:
    """'C42H82NO8P' -> Counter({'H': 82, 'C': 42, 'O': 8, 'N': 1, 'P': 1})."""
    formula = formula.strip()
    counts: Counter = Counter()
    pos = 0
    for m in _TOKEN.finditer(formula):
        if m.start() != pos:
            raise ValueError(f"can't parse formula '{formula}' near '{formula[pos:]}'")
        el, n = m.group(1), m.group(2)
        if el not in ELEMENTS:
            raise ValueError(f"unknown element '{el}' in '{formula}'")
        counts[el] += int(n) if n else 1
        pos = m.end()
    if pos != len(formula):
        raise ValueError(f"can't parse formula '{formula}' near '{formula[pos:]}'")
    return counts


def format_formula(counts: Counter) -> str:
    """Hill order: C, H, then alphabetical. Zero counts are dropped."""
    keys = [k for k in ("C", "H") if counts.get(k, 0)] + sorted(k for k in counts if k not in ("C", "H") and counts[k])
    return "".join(f"{k}{counts[k] if counts[k] != 1 else ''}" for k in keys)


def mono_mass(formula: str | Counter) -> float:
    counts = parse_formula(formula) if isinstance(formula, str) else formula
    if any(v < 0 for v in counts.values()):
        raise ValueError(f"negative atom count in {dict(counts)}")
    return float(sum(ELEMENTS[el] * n for el, n in counts.items()))


_ADDUCT = re.compile(r"^\[(\d*)M((?:[+-][0-9]*[A-Za-z0-9]+)*)\](\d*)([+-])$")
_PIECE = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)")


def parse_adduct(adduct: str) -> tuple[int, list[tuple[int, str]], int]:
    """'[M-2H]2-' -> (1, [(-2, 'H')], -2). Returns (n_M, [(signed count, formula)], charge)."""
    a = adduct.replace(" ", "")
    m = _ADDUCT.match(a)
    if not m:
        raise ValueError(f"can't parse adduct '{adduct}'")
    n_m = int(m.group(1)) if m.group(1) else 1
    pieces = []
    for sign, n, f in _PIECE.findall(m.group(2)):
        pieces.append(((1 if sign == "+" else -1) * (int(n) if n else 1), f))
    z = int(m.group(3)) if m.group(3) else 1
    z = z if m.group(4) == "+" else -z
    return n_m, pieces, z


def ion_mz(formula: str, adduct: str) -> float:
    """m/z of an ion: (n*M + pieces - z*electron) / |z|."""
    n_m, pieces, z = parse_adduct(adduct)
    mass = n_m * mono_mass(formula)
    for count, f in pieces:
        mass += count * mono_mass(f)
    mass -= z * ELECTRON          # a + charge means electrons were lost
    return mass / abs(z)


def ion_formula(formula: str, adduct: str) -> str:
    n_m, pieces, _ = parse_adduct(adduct)
    c = Counter({k: v * n_m for k, v in parse_formula(formula).items()})
    for count, f in pieces:
        for el, n in parse_formula(f).items():
            c[el] += count * n
    return format_formula(c)


def adduct_polarity(adduct: str) -> str:
    return "positive" if parse_adduct(adduct)[2] > 0 else "negative"


def ppm(observed, reference):
    return 1e6 * (observed - reference) / reference


def shift_formula(formula: str, d_carbons: int, d_double_bonds: int, d_oxygen: int = 0) -> str:
    """Change chain length / unsaturation of a sum composition.

    Adding a CH2 adds one carbon and two H; each extra double bond removes two H.
    So PC 34:1 -> PC 36:2 is shift_formula('C42H82NO8P', +2, +1).
    """
    c = parse_formula(formula)
    c["C"] += d_carbons
    c["H"] += 2 * d_carbons - 2 * d_double_bonds
    c["O"] += d_oxygen
    return format_formula(c)


def n_carbons(formula: str) -> int:
    return parse_formula(formula).get("C", 0)


def m2_overlap_fraction(n_c: int) -> float:
    """Rough relative abundance of the M+2 (two 13C) isotope vs the monoisotopic peak.

    Binomial with 1.07% 13C: P(2)/P(0) ~ C(n,2) * p^2 / (1-p)^2. Good enough for a flag.
    """
    p = 0.0107
    return n_c * (n_c - 1) / 2 * p * p / (1 - p) ** 2


def m1_ratio(n_c: int) -> float:
    p = 0.0107
    return n_c * p / (1 - p)
