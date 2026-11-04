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
            "sections": {"n_expected": 2, "split_axis": "auto", "top_is_low_index": True, "flip_slides": [], "manual_split": {},
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


def test_side_by_side_sections():
    # the real exports: wider than tall, two sections left and right
    img, top, bot = two_blobs()
    img, top, bot = img.T.copy(), top.T, bot.T
    seg = tissue.segment_slide(img, "E", qc_cfg())
    assert seg["axis"] == "cols" and seg["n_sections"] == 2
    assert (seg["positions"][top] == tissue.TOP).mean() > 0.95      # left one = Top by default
    seg = tissue.segment_slide(img, "E", qc_cfg(top_is_low_index=False))
    assert (seg["positions"][top] == tissue.BOTTOM).mean() > 0.95


def test_touching_side_by_side_are_cut():
    img, _, _ = two_blobs(touching=True)
    seg = tissue.segment_slide(img.T.copy(), "E", qc_cfg())
    assert seg["axis"] == "cols" and seg["n_sections"] == 2


def test_old_config_keys_still_work():
    img, top, _ = two_blobs()
    cfg = qc_cfg()
    cfg["sections"] = {"n_expected": 2, "top_is_low_row_index": False, "flip_slides": [], "manual_split_row": {}}
    seg = tissue.segment_slide(img, "E", cfg)
    assert (seg["positions"][top] == tissue.BOTTOM).mean() > 0.95


def test_manual_split():
    img, _, _ = two_blobs()
    seg = tissue.segment_slide(img, "F", qc_cfg(manual_split={"F": 30}))
    assert "manual" in seg["how"]


# ---- config / sample sheet ----------------------------------------------------------------
def test_output_inside_data_root_is_refused(tmp_path):
    data = tmp_path / "data root"
    data.mkdir()
    cfg = (REPO / "config.yaml").read_text()
    cfg = cfg.replace('data_root: "/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL"', f'data_root: "{data}"')
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


def test_provenance_record(tmp_path):
    from cln3msi import provenance
    cfgp = tmp_path / "c.yaml"
    cfgp.write_text(f'inherit: "{REPO / "config.test.yaml"}"\npaths: {{output_root: "{tmp_path / "out"}"}}\n')
    ctx = load(cfgp)
    fp = provenance.code_fingerprint()
    assert len(fp) == 12 and fp == provenance.code_fingerprint()      # deterministic
    rec = provenance.write(ctx, "09_differential", ["lip_neg", "met_neg"])
    files = sorted(p.name for p in (tmp_path / "out" / "metadata" / "provenance").iterdir())
    assert "09_differential_lip_neg.json" in files and "09_differential_met_neg_config.yaml" in files
    assert rec["code_fingerprint"] == fp and rec["packages"].get("numpy")
    assert {r["mode"] for r in provenance.read_all(ctx)} == {"lip_neg", "met_neg"}


def test_sign_test_reported():
    from cln3msi import stats as st
    animals = pd.DataFrame({"slide": list("AABBCCDDEEFF"), "genotype": ["Control", "Mutant"] * 6,
                            "sex": ["F"] * 12}, index=[f"a{i}" for i in range(12)])
    log2 = pd.DataFrame({"p1": [1.0, 1.5, 1.0, 1.6, 1.0, 1.4, 1.0, 1.7, 1.0, 1.5, 1.0, 1.45]}, index=animals.index)
    r = st.paired_tests(log2, animals, "Control", "Mutant", min_pairs=5).iloc[0]
    assert abs(r["p_sign"] - 0.03125) < 1e-9 and r["n_pairs_up"] == 6   # 6 of 6 up: 2 x 0.5^6


def raster_hugging_tissue():
    """Like the real lip_neg exports: the raster barely extends past the tissue (a 2-pixel
    glass margin), and inside each section the grey-matter core is much dimmer than the
    bright white-matter rim. Plain Otsu then cuts grey matter away as "background"."""
    from scipy import ndimage
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:50, 0:90]
    secs, cores = [], []
    for cx in (22, 67):
        secs.append(((yy - 25) / 22) ** 2 + ((xx - cx) / 18) ** 2 <= 1)
        cores.append(((yy - 25) / 17) ** 2 + ((xx - cx) / 13) ** 2 <= 1)
    tissue_px = secs[0] | secs[1]
    core = cores[0] | cores[1]
    acquired = ndimage.binary_dilation(tissue_px, iterations=2)
    img = np.full(xx.shape, np.nan)
    img[acquired] = 5.0            # glass
    img[tissue_px] = 400.0         # white-matter-rich rim
    img[core] = 40.0               # dim grey-matter core, half of each section
    gap = tissue_px & (np.abs(yy - 25) <= 3)                # the rim is broken: dim tissue reaches the edge
    img[gap] = 40.0
    img[acquired] *= rng.gamma(30, 1 / 30, acquired.sum())
    return img, tissue_px


def test_tissue_mask_when_raster_hugs_tissue():
    img, truth = raster_hugging_tissue()
    m_otsu, _, _ = tissue.tissue_mask(img, 0, 10, 0, "otsu")
    m_auto, _, how = tissue.tissue_mask(img, 0, 10, 0, "auto")
    recall_otsu = (m_otsu & truth).sum() / truth.sum()
    recall_auto = (m_auto & truth).sum() / truth.sum()
    assert recall_otsu < 0.9, f"test image doesn't reproduce the problem (otsu recall {recall_otsu:.2f})"
    assert recall_auto > 0.95, (recall_auto, how)
    assert (m_auto & ~truth).sum() / truth.sum() < 0.15      # and it doesn't swallow much glass


def test_auto_keeps_otsu_with_plenty_of_glass():
    img, top, bot = two_blobs()
    m, _, how = tissue.tissue_mask(img, 0, 10, 0, "auto")
    assert how == "otsu"
    assert ((m & (top | bot)).sum() / (top | bot).sum()) > 0.95


def test_sm_in_negative_mode_is_m_minus_ch3_not_odd_chain():
    # the real-data trap: "SM 35:1 [M-H]-" has exactly the formula of SM 36:1 [M-CH3]-
    import yaml
    from cln3msi import annotate
    cfg = yaml.safe_load((REPO / "config.yaml").read_text())["annotation"]["adducts"]
    lib = annotate.build_lipid_library()
    ions = annotate.ion_table(lib, "negative", cfg)
    assert not ((ions.lipid_class == "SM") & (ions.adduct == "[M-H]-")).any()
    peaks = pd.DataFrame({"mz": [chem.ion_mz("C41H83N2O6P", "[M-CH3]-")]})
    hit = annotate.match_peaks(peaks, ions, 2.0, False, 5).iloc[0]
    assert hit["annotation"] == "SM 36:1;O2" and hit["adduct"] == "[M-CH3]-", hit.to_dict()
    pos = annotate.ion_table(lib, "positive", cfg)
    assert not (pos.lipid_class == "SHexCer").any()


def test_paired_tests_need_min_pairs():
    from cln3msi import stats as st
    animals = pd.DataFrame({"slide": list("AABBCC"), "genotype": ["Control", "Mutant"] * 3,
                            "sex": ["F"] * 6}, index=[f"a{i}" for i in range(6)])
    log2 = pd.DataFrame({"p1": [1.0, 2.0, 1.1, 2.1, 0.9, 1.9]}, index=animals.index)
    r = st.paired_tests(log2, animals, "Control", "Mutant", min_pairs=5).iloc[0]
    assert np.isnan(r["p_ttest"]) and r["n_pairs"] == 3 and abs(r["log2fc"] - 1.0) < 1e-9
    r = st.paired_tests(log2, animals, "Control", "Mutant", min_pairs=3).iloc[0]
    assert r["p_ttest"] < 0.01


def _fake_means(swap_slide=None, seed=0):
    """6 slides x (Control, Mutant), 200 peaks; the Mutant has a shared effect profile.
    swap_slide gets its two animals swapped, like a sample-sheet mix-up."""
    rng = np.random.default_rng(seed)
    slides = list("EFGIJK")
    effect = rng.normal(0, 0.4, 200)
    effect[0] = -1.0                      # p0 (a sulfatide in the class test) clearly down
    base = rng.uniform(5, 9, 200)
    rows, animals = [], []
    for s in slides:
        for g in ("Control", "Mutant"):
            aid = f"{s}_{g}"
            animals.append({"animal_id": aid, "slide": s, "genotype": g, "sex": "F" if s in "EFG" else "M"})
            truly_mut = (g == "Mutant") != (s == swap_slide)
            v = base + rng.normal(0, 0.1, 200) + (effect if truly_mut else 0)
            for unit in ("whole_section", "white_matter"):
                for k in range(200):
                    rows.append({"unit": unit, "animal_id": aid, "peak": f"p{k}",
                                 "mean_norm_intensity": 2 ** v[k], "n_pixels": 500})
    return pd.DataFrame(rows), pd.DataFrame(animals).set_index("animal_id")


def test_pair_check_flags_a_swapped_slide_and_only_that():
    from cln3msi import classes
    ML, animals = _fake_means(swap_slide="J")
    _, T, flagged = classes.consistency(ML, animals, "Control", "Mutant", ["whole_section", "white_matter"], -0.3)
    assert flagged == ["J"], T
    ML, animals = _fake_means(swap_slide=None)
    _, T, flagged = classes.consistency(ML, animals, "Control", "Mutant", ["whole_section", "white_matter"], -0.3)
    assert flagged == [], T


def test_class_totals_ratios_and_sensitivity():
    from cln3msi import classes
    ML, animals = _fake_means(swap_slide="J")
    ann = pd.DataFrame({"peak_id": [f"p{k}" for k in range(6)],
                        "annotation": ["SHexCer 42:2;O2", "SHexCer 42:1;O3", "HexCer 42:2;O2", "PE 38:4",
                                       "PE 40:6", "PE 40:6"],
                        "adduct": ["[M-H]-"] * 6, "lipid_class": ["SHexCer", "SHexCer", "HexCer", "PE", "PE", "PE"],
                        "is_isotope": [False] * 5 + [True], "n_sections": [12] * 6,
                        "candidate_tags": ["SHexCer", "SHexCer;BMP_PG", "HexCer", "", "", ""]})
    ok = classes.usable_peaks(ann, 10, "tags")
    assert list(ok["peak_id"]) == ["p0", "p2", "p3", "p4"]       # isotope and cross-tag ambiguous dropped
    cls = {"Sulfatide": {"classes": ["SHexCer"]}, "HexCer": {"classes": ["HexCer"]}, "PE": {"classes": ["PE"]},
           "Sulfatide 2-OH": {"classes": ["SHexCer"], "suffix": ";O3"}}
    vals, mem = classes.totals(ML, ok, cls, [["Sulfatide", "HexCer"], ["Sulfatide 2-OH", "PE"]], ["whole_section"])
    assert set(vals["feature"]) == {"Sulfatide", "HexCer", "PE", "Sulfatide / HexCer"}   # empty class skipped
    R = classes.test_features(vals, animals, "Control", "Mutant", 5)
    assert {"pair_E", "pair_J", "q_bh", "pct_change"} <= set(R.columns)
    r = R[R.feature == "Sulfatide"].iloc[0]
    assert np.sign(r["pair_J"]) == -np.sign(r["pair_E"])          # the swapped slide points the other way
    X = classes.sensitivity(vals, animals, "Control", "Mutant", ["J"], 5).set_index("feature")
    assert X.at["Sulfatide", "flipped_J_n_pairs_up"] in (0, 6)    # flipped: all six pairs agree
    assert X.at["Sulfatide", "without_J_n_pairs"] == 5
    assert classes.candidate_class("LPC O-16:2 [M+H]+") == "LPC"
    assert classes.candidate_class("BMP/PG 36:1 [M-H]-") == "BMP/PG"
    assert classes.candidate_class("Glycerol 3-phosphate [M-H]-") == "Glycerol 3-phosphate"
