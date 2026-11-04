"""After the synthetic run: did the pipeline find what make_synthetic_data.py planted?

    python tests/check_synthetic_results.py --config config.test.yaml

Exit code 0 = all checks passed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cln3msi.config import load  # noqa: E402

TRUE_OFFSETS = {"E": 3.0, "F": 2.0, "G": 4.0, "I": 1.5, "J": 3.5, "K": 2.5}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ctx = load(ap.parse_args().config)
    R = ctx.output_root / "results"
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))

    # 03: per-slide m/z offsets recovered (lip_neg, met_neg have the plain offsets)
    ma = json.loads((R / "03_qc_all" / "lip_neg" / "mass_accuracy.json").read_text())
    err = max(abs(ma["per_slide_offset_ppm"][s] - v) for s, v in TRUE_OFFSETS.items())
    check("03 per-slide ppm offsets within 0.5 ppm", err < 0.5, f"max error {err:.2f} ppm")
    g = pd.read_csv(R / "03_qc_all" / "lip_neg" / "gradient.csv").set_index("slide")
    check("03 slide G has the strongest gradient", g["matrix_off_pred_log2_bottom_vs_top"].idxmax() == "G",
          g["matrix_off_pred_log2_bottom_vs_top"].round(3).to_dict())
    sec = pd.read_csv(R / "03_qc_all" / "lip_neg" / "sections.csv")
    check("03 Top sections are the Controls", (sec.loc[sec.position == "Top", "genotype"] == "Control").all())
    area = sec.pivot_table(index="slide", columns="genotype", values="area_mm2")
    check("03 Mutant sections smaller (planted 8% atrophy)", (area["Mutant"] < area["Control"]).all())

    # 05: planted species annotated
    ann = pd.read_csv(ctx.output_root / "processed" / "lip_neg" / "annotations.csv")
    for name in ("Cer 42:2;O2", "SHexCer 42:2;O2", "HexCer 42:2;O2", "BMP/PG 34:1", "PI 38:4"):
        check(f"05 annotated {name}", (ann["annotation"] == name).any())
    check("05 M+1 isotopes flagged", ann["is_isotope"].sum() >= 8, f"{int(ann['is_isotope'].sum())} flagged")

    # 07: registration
    reg = pd.read_csv(R / "07_register" / "lip_neg" / "registration.csv")
    check("07 all sections registered with Dice >= 0.8", (reg["dice"] >= 0.8).all(), f"min {reg['dice'].min():.2f}")

    # 09: planted genotype effects
    def hit(mode, label_start, unit, sign):
        d = pd.read_csv(R / "09_differential" / mode / "differential.csv")
        r = d[(d["unit"] == unit) & d["label"].astype(str).str.startswith(label_start)]
        if r.empty:
            return False, "not tested"
        r = r.iloc[0]
        ok = (r["log2fc"] * sign > 0.2) and (r["q_targeted"] < 0.05 if r["targeted"] else r["q_all"] < 0.05)
        return ok, f"log2FC {r['log2fc']:+.2f}, p {r['p_ttest']:.1e}"

    for mode, lab, unit, sign in [("lip_neg", "Cer 42:2;O2", "whole_section", 1),
                                  ("lip_neg", "BMP/PG 34:1", "whole_section", 1),
                                  ("lip_neg", "HexCer 42:2;O2", "white_matter", 1),
                                  ("lip_neg", "SHexCer 42:2;O2", "white_matter", -1),
                                  ("lip_pos", "LSM 18:1;O2", "whole_section", 1),
                                  ("met_neg", "Glycerophosphoinositol", "whole_section", 1)]:
        ok, det = hit(mode, lab, unit, sign)
        check(f"09 {mode} {lab} {'up' if sign > 0 else 'down'} in {unit}", ok, det)
    d = pd.read_csv(R / "09_differential" / "lip_neg" / "differential.csv")
    r = d[(d["unit"] == "whole_section") & (d["label"] == "PI 38:4 [M-H]-")].iloc[0]
    check("09 unchanged PI 38:4 stays small", abs(r["log2fc"]) < 0.25, f"log2FC {r['log2fc']:+.2f}")

    # 09 class tests and slide consistency (labels are correct in the synthetic data)
    for mode in ("lip_neg", "lip_pos", "met_neg"):
        pc = R / "09_differential" / mode / "pair_consistency.csv"
        ok = pc.exists() and not pd.read_csv(pc)["flagged"].astype(bool).any()
        check(f"09 {mode} slide consistency: no slide flagged", ok, "" if pc.exists() else "not written")
    ct = {m: pd.read_csv(R / "09_differential" / m / "class_tests.csv").set_index(["unit", "feature"])
          for m in ("lip_neg", "met_neg")}
    for mode, unit, feat, sign in [("lip_neg", "whole_section", "Cer", 1), ("lip_neg", "white_matter", "HexCer", 1),
                                   ("lip_neg", "white_matter", "Sulfatide non-OH", -1),
                                   ("met_neg", "whole_section", "GPI", 1)]:
        r = ct[mode].loc[(unit, feat)]
        check(f"09 class test {mode} {feat} {'up' if sign > 0 else 'down'} in {unit}",
              r["log2fc"] * sign > 0.2 and r["q_bh"] < 0.05, f"log2FC {r['log2fc']:+.2f}, q {r['q_bh']:.1e}")

    # provenance: every phase recorded, all with the same code and config
    from cln3msi import provenance
    recs = [r for r in provenance.read_all(ctx) if r.get("phase") != "10_report"]
    fps = {r["code_fingerprint"] for r in recs}
    cfgs = {r["config_sha256"] for r in recs}
    check("provenance recorded for phases 00 to 09, one code version and one config",
          len({r["phase"] for r in recs}) >= 10 and len(fps) == 1 and len(cfgs) == 1,
          f"{len(recs)} records, {len(fps)} code version(s), {len(cfgs)} config(s)")

    report = ctx.output_root / "report" / "report.html"
    check("10 report written", report.exists() and report.stat().st_size > 100_000)

    w = max(len(c[0]) for c in checks)
    for name, ok, det in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {name:<{w}}  {det}")
    n_fail = sum(not ok for _, ok, _ in checks)
    print(f"\n{len(checks) - n_fail}/{len(checks)} checks passed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
