"""Small unit tests: masses, spectra helpers, tissue splitting, config and sample-sheet checks.

    python -m pytest -q tests
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import chem, spectra, tissue  # noqa: E402
from cln3msi.config import ConfigError, load  # noqa: E402
from cln3msi.imzml_io import detect_profile  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


# ---- masses: reference values from LIPID MAPS / literature --------------------------------
@pytest.mark.parametrize("formula,adduct,expected", [
    ("C42H82NO8P", "[M+H]+", 760.5851),     # PC 34:1
    ("C42H82NO8P", "[M+K]+", 798.5410),
    ("C42H82NO8P", "[M+Na]+", 782.5670),
    ("C47H83O13P", "[M-H]-", 885.5499),     # PI 38:4
    ("C48H91NO11S", "[M-H]-", 888.6240),    # SHexCer 42:2;O2 (ST 24:1)
    ("C73H131N3O31", "[M-2H]2-", 771.9311), # GM1 36:1;O2
    ("C10H10N2", "[M-H]-", 157.0771),       # DAN
])
def test_ion_mz(formula, adduct, expected):
    assert chem.ion_mz(formula, adduct) == pytest.approx(expected, abs=2e-4)


def test_shift_formula():
    # PC 34:1 -> PC 36:2
    assert chem.shift_formula("C42H82NO8P", 2, 1) == "C44H84NO8P"


def test_bad_formula():
    with pytest.raises(ValueError):
        chem.mono_mass("C5H14NO+")


# ---- spectra helpers ----------------------------------------------------------------------
def gaussian(axis, mu, fwhm, a=100.0):
    s = fwhm / 2.3548
    return a * np.exp(-0.5 * ((axis - mu) / s) ** 2)


def test_peak_center_coarse_sampling():
    # only ~4 points per FWHM, like the synthetic data; centre must be within 0.5 ppm
    axis = 800 * np.exp(np.arange(-200, 200) * np.log1p(60e-6))
    mu = 800.0123
    y = gaussian(axis, mu, mu / 4200)
    pk = spectra.pick_peaks(axis, y + 0.01, snr=3, noise_window=51)
    assert len(pk["mz"]) == 1
    assert abs(pk["mz"][0] - mu) / mu * 1e6 < 0.5


def test_window_sums_match_loop():
    rng = np.random.default_rng(0)
    mz = np.sort(rng.uniform(100, 200, 500))
    it = rng.random(500)
    lo, hi = np.array([110, 150.5]), np.array([120, 151])
    got = spectra.window_sums(mz, it, lo, hi)
    want = [it[(mz >= a) & (mz <= b)].sum() for a, b in zip(lo, hi)]
    assert np.allclose(got, want)


def test_ppmgrid_profile_interpolation_keeps_peak():
    mz = np.linspace(499.9, 500.1, 41)
    it = gaussian(mz, 500.0, 0.05)
    g = spectra.PpmGrid(400, 600, 2.0)
    acc = np.zeros(g.n)
    g.add_profile(acc, mz, it)
    j = np.argmax(acc)
    assert abs(g.centers[j] - 500.0) / 500 * 1e6 < 2.5


def test_detect_profile_vs_centroid():
    axis = np.linspace(500, 510, 4000)
    prof = sum(gaussian(axis, m, 0.02) for m in (501, 503, 507)) + 0.01
    assert detect_profile(axis, prof)["guess"] == "profile"
    rng = np.random.default_rng(1)
    cmz = np.sort(rng.uniform(500, 510, 300))
    assert detect_profile(cmz, rng.gamma(1.0, 10.0, 300))["guess"] == "centroid"


def test_align_peaks():
    t = [np.array([500.0, 600.0]), np.array([500.0005, 600.001]), np.array([700.0])]
    cl = spectra.align_peaks(t, align_ppm=3)
    assert sorted(len(c) for c in cl) == [1, 2, 2]


# ---- tissue / sections --------------------------------------------------------------------
def two_blobs(gap=6, touching=False):
    img = np.full((60, 40), 1.0)
    yy, xx = np.mgrid[0:60, 0:40]
    top = ((yy - 15) / 11) ** 2 + ((xx - 20) / 16) ** 2 <= 1
    cy = 26 if touching else 15 + 22 + gap
    bot = ((yy - cy) / 11) ** 2 + ((xx - 20) / 16) ** 2 <= 1
    img[top | bot] = 100.0
    return img, top, bot


def qc_cfg(**sections):
    return {"tissue": {"smooth_sigma_px": 0, "min_object_px": 10, "close_radius_px": 0, "erode_edge_px": 1},
            "sections": {"n_expected": 2, "top_is_low_row_index": True, "flip_slides": [], "manual_split_row": {},
                         **sections}}


def test_split_and_assign():
    img, top, bot = two_blobs()
    seg = tissue.segment_slide(img, "E", qc_cfg())
    assert seg["n_sections"] == 2
    assert (seg["positions"][top] == tissue.TOP).mean() > 0.95
    assert (seg["positions"][bot] == tissue.BOTTOM).mean() > 0.95


def test_flip_slide():
    img, top, _ = two_blobs()
    seg = tissue.segment_slide(img, "E", qc_cfg(flip_slides=["E"]))
    assert (seg["positions"][top] == tissue.BOTTOM).mean() > 0.95


def test_touching_sections_are_cut():
    img, _, _ = two_blobs(touching=True)
    seg = tissue.segment_slide(img, "E", qc_cfg())
    assert seg["n_sections"] == 2


def test_manual_split():
    img, _, _ = two_blobs()
    seg = tissue.segment_slide(img, "F", qc_cfg(manual_split_row={"F": 30}))
    assert "manual" in seg["how"]


# ---- config / sample sheet ----------------------------------------------------------------
def test_output_inside_data_root_is_refused(tmp_path):
    data = tmp_path / "data root"
    data.mkdir()
    cfg = (REPO / "config.yaml").read_text()
    cfg = cfg.replace('"/Users/neharao/Desktop/Brain - Pilar - Spatial"', f'"{data}"')
    cfg = cfg.replace('"/Users/neharao/cln3_msi_outputs"', f'"{data / "outputs"}"')
    p = tmp_path / "c.yaml"
    p.write_text(cfg.replace('sample_sheet: "metadata/sample_sheet.csv"', f'sample_sheet: "{REPO / "metadata" / "sample_sheet.csv"}"'))
    with pytest.raises(ConfigError):
        load(p)


def test_sample_sheet_validation_catches_two_controls():
    from cln3msi import samples
    ctx = load(REPO / "config.test.yaml")
    df = pd.read_csv(REPO / "metadata" / "sample_sheet.csv", dtype=str)
    samples.validate_sheet(df, ctx)            # the real sheet passes
    bad = df.copy()
    bad.loc[1, "genotype"] = "Control"         # slide E now has two Controls
    with pytest.raises(samples.SampleSheetError):
        samples.validate_sheet(bad, ctx)
