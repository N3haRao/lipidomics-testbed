"""Sample sheet: load it, check it, and match imzML files to animals.

The design rules we check here (they get asserted at the start of every phase that uses
the sheet, so a typo in the CSV can't silently mess up the statistics):
  * every slide has exactly one Control and one Mutant
  * every slide has exactly one Top and one Bottom
  * the slide letter is uppercase in the sheet and lowercase in the file names
  * each file named in the sheet appears for exactly the two animals of its slide
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

REQUIRED_COLS = ["animal_id", "genotype", "sex", "slide", "position_on_slide", "run_order"]


class SampleSheetError(RuntimeError):
    pass


def load_sheet(ctx) -> pd.DataFrame:
    if not ctx.sample_sheet.exists():
        raise SampleSheetError(f"sample sheet not found: {ctx.sample_sheet}")
    df = pd.read_csv(ctx.sample_sheet, dtype=str).apply(lambda c: c.str.strip())
    validate_sheet(df, ctx)
    df["run_order"] = df["run_order"].astype(int)
    return df


def validate_sheet(df: pd.DataFrame, ctx) -> None:
    problems = []
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SampleSheetError(f"sample sheet is missing columns {missing}")
    g = ctx.cfg["project"]["genotypes"]
    allowed = {g["control"], g["mutant"]}
    bad = set(df["genotype"]) - allowed
    if bad:
        problems.append(f"unknown genotype values {bad} (allowed {allowed}, see project.genotypes)")
    if df["animal_id"].duplicated().any():
        problems.append(f"duplicated animal_id: {df.loc[df.animal_id.duplicated(), 'animal_id'].tolist()}")
    if not df["slide"].str.fullmatch(r"[A-Z]").all():
        problems.append("slide column should be single uppercase letters")
    for slide, sub in df.groupby("slide"):
        geno = sorted(sub["genotype"])
        if geno != sorted([g["control"], g["mutant"]]):
            problems.append(f"slide {slide}: needs exactly one {g['control']} and one {g['mutant']}, has {geno}")
        pos = sorted(sub["position_on_slide"])
        if pos != ["Bottom", "Top"]:
            problems.append(f"slide {slide}: needs one Top and one Bottom, has {pos}")
        if sub["run_order"].nunique() != 1 or sub["sex"].nunique() != 1:
            problems.append(f"slide {slide}: the two animals have different run_order or sex, check the sheet")
    # file columns: imzml_<mode>
    for col in [c for c in df.columns if c.startswith("imzml_")]:
        for slide, sub in df.groupby("slide"):
            files = sub[col].dropna().unique()
            if len(files) != 1:
                problems.append(f"{col}, slide {slide}: the two animals should share one file, found {list(files)}")
                continue
            m = re.match(ctx.cfg["modes"]["file_regex"], files[0])
            if not m:
                problems.append(f"{col}, slide {slide}: '{files[0]}' does not match modes.file_regex")
            elif m.group("slide").upper() != slide:
                problems.append(f"{col}, slide {slide}: file '{files[0]}' is for slide {m.group('slide').upper()}")
        if df[col].dropna().duplicated().sum() != len(df[col].dropna()) // 2:
            problems.append(f"{col}: each file should be listed exactly twice (one Control, one Mutant)")
    if problems:
        raise SampleSheetError("sample sheet problems:\n  - " + "\n  - ".join(problems))


def parse_filename(name: str, regex: str) -> dict | None:
    """slide_e_lip_neg.imzML -> {'slide': 'E', 'mode': 'lip_neg'}; None if it doesn't match."""
    m = re.match(regex, name)
    if not m:
        return None
    return {"slide": m.group("slide").upper(), "mode": f"{m.group('analyte')}_{m.group('pol')}"}


def files_for_mode(ctx, mode: str, sheet: pd.DataFrame | None = None) -> dict[str, Path]:
    """slide letter -> imzML path for one mode. Uses the sheet's file column if it has one."""
    sheet = load_sheet(ctx) if sheet is None else sheet
    col = f"imzml_{mode}"
    folder = ctx.mode_dir(mode)
    out = {}
    if col in sheet.columns:
        for slide, sub in sheet.groupby("slide"):
            out[slide] = folder / sub[col].iloc[0]
    else:   # fall back to globbing with the naming convention
        for p in sorted(folder.glob("*.imzML")):
            info = parse_filename(p.name, ctx.cfg["modes"]["file_regex"])
            if info and info["mode"] == mode:
                out[info["slide"]] = p
    return out


def section_table(sheet: pd.DataFrame) -> pd.DataFrame:
    """One row per (slide, position): who is on top, who is on the bottom."""
    cols = ["slide", "position_on_slide", "animal_id", "genotype", "sex", "run_order"]
    return sheet[cols].sort_values(["slide", "position_on_slide"], ascending=[True, False]).reset_index(drop=True)
