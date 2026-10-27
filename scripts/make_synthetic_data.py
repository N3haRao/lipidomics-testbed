"""Make a small fake dataset with the same names and layout as the real one.

    python scripts/make_synthetic_data.py --out synthetic_data

What you get (all small, ~400 MB total, a minute or two to write):

  synthetic_data/
    Brain - Pilar - Spatial (synthetic)/        <- use as data_root
      imzML files for biostats/
        Lipid Negative/       slide_{e,f,g,i,j,k}_lip_neg.imzML/.ibd   continuous + PROFILE
        Lipid Positive/       slide_..._lip_pos.imzML/.ibd             processed + CENTROID
        Metabolite Negative/  slide_..._met_neg.imzML/.ibd             processed + PROFILE
      export_All_Lip_Neg.csv   a fake SCiLS-style export (semicolon CSV with # comments)
      export_All_Lip_Pos.xlsx  same idea as an Excel file
      optical_slide_E.tif      a fake optical image for slide E (tests the TIFF path)
    synthetic_atlas.npz        a toy "Allen" atlas (annotation volume + structure tree)
    lmsd_fake.zip, hmdb_fake.zip   tiny SDF files so 00_setup can be tested offline
    truth.csv                  the planted effects, so you can check the pipeline finds them

The three modes use the three storage/spectrum combinations on purpose, so every branch of
the pipeline (memmap fast path, binning, centroid vs profile) gets exercised.

What's planted in the biology:
  * two sections per slide, Control on top, Mutant (smaller, ~8% "atrophy") on the bottom
  * white matter (corpus callosum-like band) is rich in sulfatides and HexCer
  * Mutant: Cer, BMP/PG, HexCer, lyso-SM and glycerophosphodiesters go UP,
            sulfatide SHexCer 42:2;O2 goes DOWN in white matter
  * every slide has its own intensity scale (slide effect) and m/z offset (2 to 4 ppm)
  * a mild top-to-bottom intensity gradient on every slide, stronger on slide G
"""
from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cln3msi import chem  # noqa: E402

SLIDES = ["E", "F", "G", "I", "J", "K"]
MZ_OFFSET_PPM = {"E": 3.0, "F": 2.0, "G": 4.0, "I": 1.5, "J": 3.5, "K": 2.5}
SLIDE_SCALE = {"E": 1.0, "F": 0.75, "G": 1.3, "I": 0.9, "J": 1.15, "K": 0.8}
GRADIENT = {"E": 0.08, "F": 0.05, "G": 0.30, "I": 0.06, "J": 0.08, "K": 0.04}
GRID_H, GRID_W = 56, 40            # MSI pixels (rows, cols), 50 um each
ATLAS_RES_UM = 25

# Synthetic structure tree with real Allen acronyms (ids are made up for the toy atlas)
STRUCTURES = [
    (997, "root", "root", [997]),
    (8, "grey", "Basic cell groups and regions", [997, 8]),
    (315, "Isocortex", "Isocortex", [997, 8, 315]),
    (1089, "HPF", "Hippocampal formation", [997, 8, 1089]),
    (549, "TH", "Thalamus", [997, 8, 549]),
    (1097, "HY", "Hypothalamus", [997, 8, 1097]),
    (1009, "fiber tracts", "fiber tracts", [997, 1009]),
    (776, "cc", "corpus callosum", [997, 1009, 776]),
    (484682516, "ccb", "corpus callosum, body", [997, 1009, 776, 484682516]),
    (484682520, "ccs", "corpus callosum, splenium", [997, 1009, 776, 484682520]),
    (1000, "alv", "alveus", [997, 1009, 1000]),
]
ID = {s[1]: s[0] for s in STRUCTURES}


# ---------------------------------------------------------------------------------------
# Toy atlas
# ---------------------------------------------------------------------------------------
def atlas_slice(ap_frac: float, dv: int = 56, ml: int = 80) -> np.ndarray:
    """A coronal toy brain at 25 um. ap_frac in [0,1] changes hippocampus size a bit."""
    yy, xx = np.mgrid[0:dv, 0:ml].astype(float)
    cy, cx = dv * 0.5, ml * 0.5
    ry, rx = dv * 0.42, ml * 0.45
    r = ((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2
    lab = np.zeros((dv, ml), np.uint32)
    brain = r <= 1
    lab[brain] = ID["HY"]
    # cortex = outer dorsal rim
    lab[brain & (r > 0.55) & (yy < cy + 0.15 * ry)] = ID["Isocortex"]
    # corpus callosum = an arc under the cortex
    arc = brain & (r > 0.40) & (r <= 0.55) & (yy < cy + 0.05 * ry)
    lab[arc & (np.abs(xx - cx) < rx * 0.35)] = ID["ccb"]
    lab[arc & (np.abs(xx - cx) >= rx * 0.35)] = ID["ccs"]
    # hippocampus = two blobs under the arc, a little bigger further back
    hr = 0.16 + 0.06 * ap_frac
    for sx in (-1, 1):
        h = ((yy - (cy - 0.05 * ry)) / (ry * hr * 1.1)) ** 2 + ((xx - (cx + sx * 0.3 * rx)) / (rx * hr * 1.6)) ** 2
        lab[brain & (h <= 1) & (lab == ID["HY"])] = ID["HPF"]
        alv = brain & (h > 1) & (h <= 1.5) & (lab == ID["HY"]) & (yy < cy)
        lab[alv] = ID["alv"]
    # thalamus = central block
    th = ((yy - (cy + 0.25 * ry)) / (ry * 0.28)) ** 2 + ((xx - cx) / (rx * 0.30)) ** 2 <= 1
    lab[brain & th & (lab == ID["HY"])] = ID["TH"]
    return lab


def make_atlas(path: Path, n_ap: int = 40):
    ann = np.stack([atlas_slice(i / (n_ap - 1)) for i in range(n_ap)])     # (AP, DV, ML)
    ref = np.zeros(ann.shape, np.float32)
    ref[ann > 0] = 0.4
    wm = np.isin(ann, [ID["ccb"], ID["ccs"], ID["alv"]])
    ref[wm] = 1.0
    structs = [{"id": i, "acronym": a, "name": n, "structure_id_path": p} for i, a, n, p in STRUCTURES]
    np.savez_compressed(path, annotation=ann, reference=ref, resolution_um=ATLAS_RES_UM,
                        structures=json.dumps(structs))
    return ann


# ---------------------------------------------------------------------------------------
# Section geometry on the slide
# ---------------------------------------------------------------------------------------
def place_section(atlas_lab: np.ndarray, center_rc, scale: float, angle_deg: float):
    """Downsample the 25 um atlas labels to 50 um, scale/rotate and drop them on the MSI grid."""
    small = atlas_lab[::2, ::2]
    out = np.zeros((GRID_H, GRID_W), np.uint32)
    rr, cc = np.mgrid[0:GRID_H, 0:GRID_W].astype(float)
    th = np.deg2rad(angle_deg)
    dy, dx = (rr - center_rc[0]) / scale, (cc - center_rc[1]) / scale
    sy = np.cos(th) * dy - np.sin(th) * dx + small.shape[0] / 2
    sx = np.sin(th) * dy + np.cos(th) * dx + small.shape[1] / 2
    iy, ix = np.round(sy).astype(int), np.round(sx).astype(int)
    ok = (iy >= 0) & (iy < small.shape[0]) & (ix >= 0) & (ix < small.shape[1])
    out[ok] = small[iy[ok], ix[ok]]
    return out


# ---------------------------------------------------------------------------------------
# Ion lists (formula, adduct) with abundance per region and genotype effects
# ---------------------------------------------------------------------------------------
# region weights: WM (fiber tracts), CTX, HPF, TH, other grey
def ions_for(mode: str, matrix_formula: str):
    W = dict
    if mode == "lip_neg":
        ions = [
            ("PE 38:4", "C43H78NO8P", "[M-H]-", 400, W(wm=0.6, ctx=1.0, hpf=1.1, th=1.0, oth=1.0), 1.0, None),
            ("PE 40:6", "C45H78NO8P", "[M-H]-", 600, W(wm=0.5, ctx=1.2, hpf=1.0, th=1.0, oth=1.0), 1.0, None),
            ("PS 40:6", "C46H78NO10P", "[M-H]-", 300, W(wm=0.7, ctx=1.1, hpf=1.0, th=1.2, oth=1.0), 1.0, None),
            ("PI 38:4", "C47H83O13P", "[M-H]-", 800, W(wm=0.4, ctx=1.2, hpf=1.3, th=1.0, oth=1.0), 1.0, None),
            ("SHexCer 42:2;O2", "C48H91NO11S", "[M-H]-", 900, W(wm=3.0, ctx=0.4, hpf=0.5, th=1.0, oth=0.6), 0.55, "wm"),
            ("SHexCer 42:1;O2", "C48H93NO11S", "[M-H]-", 250, W(wm=3.0, ctx=0.4, hpf=0.5, th=1.0, oth=0.6), 1.0, None),
            ("SHexCer 42:2;O3", "C48H91NO12S", "[M-H]-", 250, W(wm=3.0, ctx=0.3, hpf=0.4, th=0.9, oth=0.5), 1.0, None),
            ("HexCer 42:2;O2", "C48H91NO8", "[M-H]-", 150, W(wm=3.0, ctx=0.4, hpf=0.5, th=1.0, oth=0.6), 1.6, None),
            ("Cer 42:2;O2", "C42H81NO3", "[M-H]-", 100, W(wm=1.2, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 2.0, None),
            ("BMP/PG 34:1", "C40H77O10P", "[M-H]-", 90, W(wm=0.8, ctx=1.0, hpf=1.2, th=1.0, oth=1.0), 1.8, None),
            ("unknown 700.12", None, 700.1234, 150, W(wm=0.3, ctx=0.5, hpf=3.0, th=0.6, oth=0.5), 1.0, None),
        ]
        matrix = [("DAN 4M-H", matrix_formula, "[4M-H]-", 150)]
        return ions, matrix
    if mode == "lip_pos":
        ions = [
            ("PC 34:1 H", "C42H82NO8P", "[M+H]+", 500, W(wm=0.6, ctx=1.2, hpf=1.1, th=1.0, oth=1.0), 1.0, None),
            ("PC 34:1 Na", "C42H82NO8P", "[M+Na]+", 300, W(wm=0.6, ctx=1.2, hpf=1.1, th=1.0, oth=1.0), 1.0, None),
            ("PC 34:1 K", "C42H82NO8P", "[M+K]+", 900, W(wm=0.6, ctx=1.2, hpf=1.1, th=1.0, oth=1.0), 1.0, None),
            ("PC 32:0 K", "C40H80NO8P", "[M+K]+", 500, W(wm=0.5, ctx=1.3, hpf=1.0, th=1.0, oth=1.0), 1.0, None),
            ("PC 36:1 K", "C44H86NO8P", "[M+K]+", 400, W(wm=1.5, ctx=0.8, hpf=0.9, th=1.0, oth=1.0), 1.0, None),
            ("PC 38:6 K", "C46H80NO8P", "[M+K]+", 300, W(wm=0.5, ctx=1.0, hpf=1.4, th=1.0, oth=1.0), 1.0, None),
            ("SM 36:1;O2 H", "C41H83N2O6P", "[M+H]+", 250, W(wm=1.0, ctx=1.0, hpf=1.0, th=1.2, oth=1.0), 1.5, None),
            ("HexCer 42:2;O2 Na", "C48H91NO8", "[M+Na]+", 200, W(wm=3.0, ctx=0.4, hpf=0.5, th=1.0, oth=0.6), 1.6, None),
            ("LSM 18:1;O2 H", "C23H49N2O5P", "[M+H]+", 40, W(wm=1.0, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 2.5, None),
        ]
        matrix = [("DAN 2M+H", matrix_formula, "[2M+H]+", 300), ("DAN 3M+H", matrix_formula, "[3M+H]+", 150)]
        return ions, matrix
    if mode == "met_neg":
        ions = [
            ("Taurine", "C2H7NO3S", "[M-H]-", 700, W(wm=0.6, ctx=1.0, hpf=1.2, th=1.0, oth=1.0), 1.0, None),
            ("Glutamate", "C5H9NO4", "[M-H]-", 500, W(wm=0.5, ctx=1.2, hpf=1.1, th=1.0, oth=1.0), 1.0, None),
            ("N-acetylaspartate", "C6H9NO5", "[M-H]-", 600, W(wm=0.8, ctx=1.0, hpf=1.0, th=1.2, oth=1.0), 1.0, None),
            ("Glutathione", "C10H17N3O6S", "[M-H]-", 200, W(wm=1.0, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 1.0, None),
            ("AMP", "C10H14N5O7P", "[M-H]-", 150, W(wm=0.7, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 1.0, None),
            ("Glycerophosphoinositol", "C9H19O11P", "[M-H]-", 80, W(wm=1.0, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 1.7, None),
            ("Glycerophosphoethanolamine", "C5H14NO6P", "[M-H]-", 80, W(wm=1.0, ctx=1.0, hpf=1.0, th=1.0, oth=1.0), 1.6, None),
            ("Cholesterol sulfate", "C27H46O4S", "[M-H]-", 120, W(wm=2.0, ctx=0.7, hpf=0.8, th=1.0, oth=1.0), 1.0, None),
        ]
        matrix = [("DAN M-H", matrix_formula, "[M-H]-", 500), ("DAN 2M-H", matrix_formula, "[2M-H]-", 250)]
        return ions, matrix
    raise ValueError(mode)


MODE_SETUP = {
    # mode: (storage, spectrum type, polarity, (mz_lo, mz_hi) of the profile axis, R)
    "lip_neg": ("continuous", "profile", "negative", (600.0, 950.0), 4200.0),
    "lip_pos": ("processed", "centroid", "positive", (300.0, 1000.0), None),
    "met_neg": ("processed", "profile", "negative", (100.0, 500.0), 6000.0),
}

REGION_OF = {ID["Isocortex"]: "ctx", ID["HPF"]: "hpf", ID["TH"]: "th", ID["HY"]: "oth",
             ID["ccb"]: "wm", ID["ccs"]: "wm", ID["alv"]: "wm"}


def ion_list_with_isotopes(ions, matrix):
    """Expand each ion into M, M+1, M+2 with rough binomial ratios."""
    out = []
    for name, f, add, base, reg, fold, where in ions:
        mz = chem.ion_mz(f, add) if f else float(add)
        n_c = chem.n_carbons(f) if f else int(mz / 14.5)
        for k, ratio in ((0, 1.0), (1, chem.m1_ratio(n_c)), (2, chem.m2_overlap_fraction(n_c))):
            out.append(dict(name=name, iso=k, mz=mz + k * chem.C13_SHIFT, base=base * ratio, reg=reg,
                            fold=fold, where=where, matrix=False))
    for name, f, add, base in matrix:
        out.append(dict(name=name, iso=0, mz=chem.ion_mz(f, add), base=base, reg=None, fold=1.0,
                        where=None, matrix=True))
    return out


def pixel_intensities(lst, region_img, sec_img, geno_of_sec, rng, slide, r):
    """Per pixel (flattened rows) intensity of each ion in lst. Returns (n_pix, n_ion)."""
    n = region_img.size
    reg_flat = region_img.ravel()
    sec_flat = sec_img.ravel()
    rows = np.repeat(np.arange(GRID_H), GRID_W)
    grad = 1 + GRADIENT[slide] * (rows / GRID_H - 0.5)          # technical top-to-bottom trend
    animal = {1: rng.normal(1, 0.08), 2: rng.normal(1, 0.08)}   # animal-to-animal noise
    out = np.zeros((n, len(lst)), np.float32)
    on = sec_flat > 0
    for j, ion in enumerate(lst):
        if ion["matrix"]:
            # Matrix ions: strong off tissue, suppressed on tissue, follows the gradient too
            val = np.where(on, 0.3, 1.0) * ion["base"]
        else:
            w = np.array([ion["reg"].get(REGION_OF.get(int(t), "oth"), 1.0) if t else 0.0 for t in reg_flat])
            val = ion["base"] * w
            for s in (1, 2):
                sel = sec_flat == s
                if geno_of_sec[s] == "Mutant" and ion["fold"] != 1.0:
                    if ion["where"] == "wm":
                        is_wm = np.array([REGION_OF.get(int(t)) == "wm" for t in reg_flat[sel]])
                        val[sel] = np.where(is_wm, val[sel] * ion["fold"], val[sel])
                    else:
                        val[sel] *= ion["fold"]
                val[sel] *= animal[s]
        noise = rng.gamma(20.0, 1 / 20.0, n)      # pixel-to-pixel noise, CV ~22%
        out[:, j] = val * noise * grad * SLIDE_SCALE[slide]
    return out


def write_slide(path: Path, mode: str, slide: str, region_img, sec_img, geno_of_sec, rng, matrix_formula):
    from pyimzml.ImzMLWriter import ImzMLWriter
    storage, stype, pol, (lo, hi), R = MODE_SETUP[mode]
    ions, matrix = ions_for(mode, matrix_formula)
    lst = ion_list_with_isotopes(ions, matrix)
    I = pixel_intensities(lst, region_img, sec_img, geno_of_sec, rng, slide, region_img)
    true_mz = np.array([d["mz"] for d in lst])
    off = MZ_OFFSET_PPM[slide] + (0.5 if mode == "lip_pos" else 0.0)
    obs_mz = true_mz * (1 + off * 1e-6)

    if storage == "continuous":
        step = np.log1p(60e-6)
        axis = lo * np.exp(np.arange(int(np.log(hi / lo) / step)) * step)
    with ImzMLWriter(str(path), mode=storage, spec_type=stype, polarity=pol,
                     mz_dtype=np.float64, intensity_dtype=np.float32) as w:
        k = 0
        for r in range(GRID_H):
            for c in range(GRID_W):
                inten = I[k]
                jitter = obs_mz * (1 + rng.normal(0, 1.0, len(obs_mz)) * 1e-6)
                keep = (jitter > lo) & (jitter < hi)
                if stype == "centroid":
                    m = jitter[keep]
                    it = inten[keep]
                    order = np.argsort(m)
                    m, it = m[order], it[order]
                    # a few random noise centroids, like a real peak list
                    nm = np.sort(rng.uniform(lo, hi, 15))
                    m = np.concatenate([m, nm])
                    it = np.concatenate([it, rng.gamma(1.0, 2.0, 15)]).astype(np.float32)
                    o = np.argsort(m)
                    m, it = m[o], it[o]
                    sel = it > 0.5
                    w.addSpectrum(m[sel], it[sel], (c + 1, r + 1, 1))
                else:
                    sig = jitter / R / 2.3548
                    if storage == "continuous":
                        y = np.abs(rng.normal(0, 0.4, len(axis))).astype(np.float32)
                        for mu, s, a in zip(jitter[keep], sig[keep], inten[keep]):
                            i0, i1 = np.searchsorted(axis, [mu - 5 * s, mu + 5 * s])
                            y[i0:i1] += a * np.exp(-0.5 * ((axis[i0:i1] - mu) / s) ** 2)
                        w.addSpectrum(axis, y, (c + 1, r + 1, 1))
                    else:   # processed profile: only points near peaks are stored
                        ms, ys = [], []
                        for mu, s, a in zip(jitter[keep], sig[keep], inten[keep]):
                            x = mu + np.arange(-16, 17) * (s / 3.5)
                            ms.append(x)
                            ys.append(a * np.exp(-0.5 * ((x - mu) / s) ** 2) + np.abs(rng.normal(0, 0.3, len(x))))
                        m = np.concatenate(ms)
                        y = np.concatenate(ys).astype(np.float32)
                        o = np.argsort(m)
                        w.addSpectrum(m[o], y[o], (c + 1, r + 1, 1))
                k += 1
    return lst


def fake_sdf(records: list[dict]) -> str:
    out = []
    for rec in records:
        out.append("\n  fake\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END")
        for k, v in rec.items():
            out.append(f"> <{k}>\n{v}\n")
        out.append("$$$$")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="synthetic_data")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--sheet", default=str(Path(__file__).resolve().parents[1] / "metadata" / "sample_sheet.csv"))
    ap.add_argument("--matrix-formula", default="C10H10N2")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    out = Path(args.out).resolve()
    data_root = out / "Brain - Pilar - Spatial (synthetic)"
    folders = {"lip_neg": "Lipid Negative", "lip_pos": "Lipid Positive", "met_neg": "Metabolite Negative"}
    for f in folders.values():
        (data_root / "imzML files for biostats" / f).mkdir(parents=True, exist_ok=True)

    sheet = pd.read_csv(args.sheet)
    atlas = make_atlas(out / "synthetic_atlas.npz")
    ap_index = 20                                   # all sections cut at this toy AP level
    truth = []
    for slide in SLIDES:
        sub = sheet[sheet.slide == slide]
        geno = {1: sub.loc[sub.position_on_slide == "Top", "genotype"].iloc[0],
                2: sub.loc[sub.position_on_slide == "Bottom", "genotype"].iloc[0]}
        top = place_section(atlas[ap_index], (14, 20), 1.0 if geno[1] == "Control" else 0.92, rng.uniform(-4, 4))
        bot = place_section(atlas[ap_index + rng.integers(-2, 3)], (42, 20),
                            0.92 if geno[2] == "Mutant" else 1.0, rng.uniform(-4, 4))
        region = np.where(top > 0, top, bot)
        sec = np.where(top > 0, 1, np.where(bot > 0, 2, 0))
        for mode, folder in folders.items():
            p = data_root / "imzML files for biostats" / folder / f"slide_{slide.lower()}_{mode}.imzML"
            lst = write_slide(p, mode, slide, region, sec, geno, rng, args.matrix_formula)
            print(f"wrote {p.name} ({p.with_suffix('.ibd').stat().st_size / 1e6:.1f} MB)")
            if slide == "E":
                for d in lst:
                    if d["iso"] == 0:
                        truth.append({"mode": mode, "name": d["name"], "mz": round(d["mz"], 5),
                                      "mutant_fold": d["fold"], "where": d["where"] or "everywhere",
                                      "matrix": d["matrix"]})
        if slide == "E":
            # Fake optical image, 10 um/px, same field of view as the MSI grid
            import tifffile
            big = np.kron(np.where(np.isin(region, [ID["ccb"], ID["ccs"], ID["alv"]]), 200,
                                   np.where(region > 0, 120, 30)), np.ones((5, 5))).astype(np.uint8)
            big = ndimage.gaussian_filter(big, 2) + rng.integers(0, 10, big.shape).astype(np.uint8)
            tifffile.imwrite(data_root / "optical_slide_E.tif", big)

    pd.DataFrame(truth).to_csv(out / "truth.csv", index=False)

    # Fake SCiLS exports: we don't know the real layout, so make two different ones
    feats = [d for d in truth if d["mode"] == "lip_neg"]
    with open(data_root / "export_All_Lip_Neg.csv", "w") as fh:
        fh.write("# Exported from fake SCiLS Lab\n# Region: All\n")
        fh.write("Spot index;x;y;" + ";".join(f"{d['mz'] * (1 + 3e-6):.4f}" for d in feats) + "\n")
        for i in range(30):
            fh.write(f"{i};{i % 10};{i // 10};" + ";".join(f"{rng.gamma(2, 50):.2f}" for _ in feats) + "\n")
    fp = pd.DataFrame({"m/z": [d["mz"] for d in truth if d["mode"] == "lip_pos"],
                       "Interval Width (+/- Da)": 0.01, "Name": [d["name"] for d in truth if d["mode"] == "lip_pos"]})
    fp.to_excel(data_root / "export_All_Lip_Pos.xlsx", index=False)

    # Fake LMSD / HMDB as zipped SDF (a handful of formulas is enough to test parsing)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from cln3msi.annotate import build_lipid_library
    lib = build_lipid_library()
    keep = lib.sample(n=min(3000, len(lib)), random_state=1)
    planted = {d["name"].split(" H")[0].split(" K")[0].split(" Na")[0] for d in truth}
    keep = pd.concat([keep, lib[lib["name"].isin(planted)]]).drop_duplicates("formula")
    recs = [{"LM_ID": f"LMFAKE{i:06d}", "COMMON_NAME": n, "FORMULA": f, "MAIN_CLASS": c}
            for i, (n, f, c) in enumerate(zip(keep["name"], keep["formula"], keep["lipid_class"]))]
    with zipfile.ZipFile(out / "lmsd_fake.zip", "w") as z:
        z.writestr("structures.sdf", fake_sdf(recs))
    mets = pd.read_csv(Path(__file__).resolve().parents[1] / "reference" / "metabolites_curated.csv")
    recs = [{"DATABASE_ID": f"HMDBFAKE{i:05d}", "GENERIC_NAME": n, "FORMULA": f}
            for i, (n, f) in enumerate(zip(mets["name"], mets["formula"]))]
    recs.append({"DATABASE_ID": "HMDBFAKE99999", "GENERIC_NAME": "Charged thing", "FORMULA": "C5H14NO+"})
    with zipfile.ZipFile(out / "hmdb_fake.zip", "w") as z:
        z.writestr("structures.sdf", fake_sdf(recs))
    print(f"\nsynthetic data_root: {data_root}")


if __name__ == "__main__":
    main()
