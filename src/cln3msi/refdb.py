"""Reference databases: download once (00_setup --download), parse, cache as parquet.

LIPID MAPS LMSD and HMDB both ship as zipped SDF files. An SDF record looks like

    <molblock...>
    > <FORMULA>
    C42H82NO8P

    > <COMMON_NAME>
    ...
    $$$$

We don't need structures, only names and formulas, so the parser just reads the "> <FIELD>"
blocks and skips the molblocks. Field names differ between databases and versions, so we
look for a few likely names and log which ones were found. Masses are recomputed from the
formula with chem.py, so we don't depend on how a database rounded them.

Nothing here uploads anything. The only network access is the plain HTTPS GET in
download(), and only when you run 00_setup.py with --download.
"""
from __future__ import annotations

import io
import shutil
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import chem

NAME_FIELDS = ["COMMON_NAME", "NAME", "GENERIC_NAME", "SYSTEMATIC_NAME", "ABBREVIATION"]
FORMULA_FIELDS = ["FORMULA", "MOLECULAR_FORMULA", "CHEMICAL_FORMULA"]
ID_FIELDS = ["LM_ID", "HMDB_ID", "DATABASE_ID", "ID"]
CLASS_FIELDS = ["MAIN_CLASS", "SUB_CLASS", "CATEGORY", "SUPER_CLASS", "CLASS"]


def download(url: str, dest: Path, timeout: int, log) -> Path:
    """Plain download to dest (skipped if it is already there)."""
    if dest.exists() and dest.stat().st_size > 0:
        log.info(f"already cached: {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    log.info(f"downloading {url}\n   -> {dest}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "cln3msi-setup/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 20)
    tmp.rename(dest)
    log.info(f"downloaded {dest.stat().st_size / 1e6:.1f} MB")
    return dest


def iter_sdf_records(fh):
    """Yield {field: value} dicts from an SDF text stream (molblocks ignored)."""
    rec: dict = {}
    field = None
    buf: list[str] = []
    for raw in fh:
        line = raw.rstrip("\r\n")
        if line.startswith("$$$$"):
            if field is not None:
                rec[field] = "\n".join(buf).strip()
            if rec:
                yield rec
            rec, field, buf = {}, None, []
            continue
        if line.startswith(">") and "<" in line and ">" in line[1:]:
            if field is not None:
                rec[field] = "\n".join(buf).strip()
            field = line[line.index("<") + 1: line.rindex(">")].strip().upper()
            buf = []
            continue
        if field is not None:
            if line == "":
                rec[field] = "\n".join(buf).strip()
                field, buf = None, []
            else:
                buf.append(line)


def _open_sdf(path: Path):
    """Return a text stream for an .sdf, or the first .sdf inside a .zip."""
    if path.suffix.lower() == ".zip":
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if n.lower().endswith(".sdf")]
        if not names:
            raise ValueError(f"no .sdf inside {path}: {z.namelist()[:10]}")
        return io.TextIOWrapper(z.open(names[0]), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def sdf_to_table(path: Path, log) -> pd.DataFrame:
    """Parse names + formulas from an SDF (zipped or not) into a DataFrame."""
    rows = []
    seen_fields: set = set()
    with _open_sdf(path) as fh:
        for rec in iter_sdf_records(fh):
            seen_fields.update(rec)
            f = next((rec[k] for k in FORMULA_FIELDS if rec.get(k)), None)
            if not f:
                continue
            rows.append({
                "id": next((rec[k] for k in ID_FIELDS if rec.get(k)), ""),
                "name": next((rec[k] for k in NAME_FIELDS if rec.get(k)), ""),
                "abbreviation": rec.get("ABBREVIATION", ""),
                "db_class": next((rec[k] for k in CLASS_FIELDS if rec.get(k)), ""),
                "formula": f.replace(" ", ""),
            })
    log.info(f"{path.name}: {len(rows)} records with a formula. Fields seen: {sorted(seen_fields)[:40]}")
    df = pd.DataFrame(rows)
    # Recompute mass from formula. Charged or odd formulas (with + / - or isotopes) are dropped.
    masses = []
    for f in df["formula"]:
        try:
            masses.append(chem.mono_mass(f))
        except ValueError:
            masses.append(float("nan"))
    df["neutral_mass"] = masses
    bad = df["neutral_mass"].isna().sum()
    if bad:
        log.info(f"  {bad} records had formulas we can't parse (charges, isotopes, R groups); dropped")
    return df.dropna(subset=["neutral_mass"]).reset_index(drop=True)


def swisslipids_to_table(path: Path, log) -> pd.DataFrame:
    """SwissLipids lipids.tsv(.gz) -> one row per "Species" (sum composition, e.g. PC(34:1)).

    SwissLipids gives the formula of the form present at pH 7.3 plus its charge; the neutral
    formula is that formula with as many H added as the charge is negative (and removed if
    positive). We check it against the file's own neutral exact mass. The finer levels
    (molecular / structural / isomeric subspecies, ~770,000 rows) have the same formulas as
    their species, so for MS1 annotation the ~10,000 species are all we need.
    """
    import re
    cols = ["Lipid ID", "Level", "Name", "Abbreviation*", "Formula (pH7.3)", "Charge (pH7.3)",
            "Exact Mass (neutral form)", "LIPID MAPS"]
    d = pd.read_csv(path, sep="\t", usecols=cols, dtype=str, encoding="latin-1", quoting=3)
    d = d[d["Level"] == "Species"].dropna(subset=["Formula (pH7.3)"])
    rows = []
    for _, r in d.iterrows():
        try:
            c = chem.parse_formula(r["Formula (pH7.3)"])
            c["H"] += -int(float(r["Charge (pH7.3)"] or 0))
            f = chem.format_formula(c)
        except Exception:
            continue
        ab = str(r["Abbreviation*"] or r["Name"]).split(" | ")[0].strip()
        m = re.match(r"^([A-Za-z0-9\-]+)\((.*)\)$", ab)
        name = f"{m.group(1)} {m.group(2)}" if m else ab
        rows.append({"id": r["Lipid ID"], "name": name, "abbreviation": ab,
                     "db_class": m.group(1) if m else "", "formula": f,
                     "mass_file": float(r["Exact Mass (neutral form)"]) if r["Exact Mass (neutral form)"] else np.nan})
    df = pd.DataFrame(rows)
    df["neutral_mass"] = [chem.mono_mass(f) for f in df["formula"]]
    dev = (df["neutral_mass"] - df["mass_file"]).abs()
    log.info(f"{path.name}: {len(df)} SwissLipids species; neutral mass vs the file's: median "
             f"|diff| {dev.median():.5f} Da, {int((dev > 0.01).sum())} entries off by > 0.01 Da (dropped)")
    return df[~(dev > 0.01)].drop(columns="mass_file").reset_index(drop=True)


def cached_table(ctx, which: str) -> pd.DataFrame | None:
    """Load the parsed cache for 'lipidmaps', 'hmdb' or 'swisslipids', or None if 00_setup hasn't made it."""
    p = ctx.cache_root / which / f"{which}.parquet"
    return pd.read_parquet(p) if p.exists() else None
