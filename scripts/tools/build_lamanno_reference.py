"""Build reference/lamanno2026/ from the Lipid Brain Atlas supplementary tables.

    python scripts/tools/build_lamanno_reference.py \
        --supp1 LaManno_SuppTable1.xlsx --supp2 LaManno_SuppTable2.xlsx

Paper: "The lipidomic architecture of the mouse brain", Nature (2026),
doi:10.1038/s41586-026-11050-0 (CC BY 4.0). Run once; the outputs are committed, so the
pipeline never needs the Excel files. Writes:

  species_sources.csv     species, lipid_class, source (msms / lcms / fitzner / vandenbosch)
  quant_whole_brain.csv   species, lipid_class, mean_nmol_pct, n_entries (8 mice, summed over
                          the chain isomers of one sum composition)
  maldi_peaks.csv         their MALDI peaks (positive mode, DHB): mz, annotation, score, status
  validation.csv          our scoring applied to their own peaks, next to their score
  README.md

Which tab feeds which source (Supplementary Table 1):
  QualitativeTANDEM_Goslin   in-house LC-MS/MS                      -> msms
  QuantitativeLCMS           targeted HILIC-ESI-MS/MS, whole brain   -> msms (and the quant table)
  HILIC_column               HILIC LC-MS list                        -> lcms
  ReversePhase_column        reverse-phase LC-MS (Skyline)           -> lcms
  Fitzner_Goslin             Fitzner et al.                          -> fitzner
  AnalyticalChemistry_PosIonMode  Vandenbosch et al. MALDI-MSI       -> vandenbosch
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cln3msi import evidence  # noqa: E402

README = """# Lipid Brain Atlas reference (mouse brain lipids with LC-MS/MS evidence)

Made by `scripts/tools/build_lamanno_reference.py` from Supplementary Tables 1 and 2 of

> The lipidomic architecture of the mouse brain. Nature (2026). doi:10.1038/s41586-026-11050-0
> (Lipid Brain Atlas, La Manno lab, EPFL; code: github.com/lamanno-epfl/lipidbrainatlas)

The article and its supplementary information are licensed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/). These files are derived from them: names were
converted to sum composition (e.g. `Cer d18:1/24:1` -> `Cer 42:2;O2`, `PE(P-16:0/20:4)` ->
`PE O-36:5`), oxidized species were dropped, and the quantitative table was averaged over
the 8 mice.

Used by phase 05 (`annotation.brain_evidence` in config.yaml) to score our names the way the
paper scored theirs; see `src/cln3msi/evidence.py`. `validation.csv` applies our scoring to
their own MALDI peaks: {agree} of {n} names with one candidate get the same call (score > 4
or not) as theirs. Differences are expected where their METASPACE term or their own MALDI
data (both left out here) tipped the balance.

Classes per source:

{classes}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--supp1", required=True)
    ap.add_argument("--supp2", required=True)
    ap.add_argument("--out", default=str(ROOT / "reference" / "lamanno2026"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    s1 = pd.read_excel(a.supp1, sheet_name=None)
    rows, failed = [], {}

    def add(names, source, tab):
        for n in pd.Series(names).dropna().astype(str):
            c = evidence.canonical(n)
            if c:
                rows.append({"species": c, "lipid_class": evidence.class_of(c), "source": source})
            else:
                failed.setdefault(tab, []).append(n)

    add(s1["QualitativeTANDEM_Goslin"]["Species Name"], "msms", "QualitativeTANDEM_Goslin")
    Q = s1["QuantitativeLCMS"].iloc[1:]           # row 0 holds the sample IDs
    add(Q["LIPID_ID"], "msms", "QuantitativeLCMS")
    add(s1["HILIC_column"]["Lipid"], "lcms", "HILIC_column")
    rp = s1["ReversePhase_column"]["PrecursorName"].astype(str)
    add(rp[~rp.str.contains(r"\(ST\)")], "lcms", "ReversePhase_column")      # (ST) = internal standards
    add(s1["Fitzner_Goslin"]["Species Name"], "fitzner", "Fitzner_Goslin")
    add(s1["AnalyticalChemistry_PosIonMode"]["Species Level ID"], "vandenbosch", "AnalyticalChemistry_PosIonMode")
    S = pd.DataFrame(rows).drop_duplicates().sort_values(["source", "lipid_class", "species"])
    S = S[~S["lipid_class"].isin(["0", "ST"])]
    S.to_csv(out / "species_sources.csv", index=False)

    # quantitative whole-brain table: mean over mice, summed over isomers of one composition
    vals = Q.drop(columns=["LIPID_ID", "Units"]).apply(pd.to_numeric, errors="coerce")
    q = pd.DataFrame({"species": [evidence.canonical(x) for x in Q["LIPID_ID"]], "v": vals.mean(axis=1).to_numpy()})
    q = q.dropna(subset=["species"])
    q = q.groupby("species").agg(mean_nmol_pct=("v", "sum"), n_entries=("v", "size")).reset_index()
    q.insert(1, "lipid_class", [evidence.class_of(x) for x in q["species"]])
    q.to_csv(out / "quant_whole_brain.csv", index=False)

    b = pd.read_excel(a.supp2, sheet_name="detailed annotations all lipids")
    M = b[["Measured MSI m/z", "Annotation", "Score", "Status"]].copy()
    M.columns = ["mz", "annotation", "score", "status"]
    M.to_csv(out / "maldi_peaks.csv", index=False)

    # validation: our score (no METASPACE, no their-own-MALDI term) vs their Score, one-name peaks
    ref = evidence.load(out)
    ref.maldi = None
    V = M[M["annotation"].notna()].copy()
    V["names"] = [evidence.split_multi(x) for x in V["annotation"]]
    V = V[V["names"].map(len) == 1]
    V["name"] = [n[0] for n in V["names"]]
    sc = [ref.score(n) for n in V["name"]]
    V["our_score"] = [d["evidence_score"] for d in sc]
    V["our_tier"] = [d["confidence_tier"] for d in sc]
    V["their_high"] = V["score"] > 4
    V["our_high"] = V["our_score"] > 4
    V = V.drop(columns="names")
    V.to_csv(out / "validation.csv", index=False)
    agree = int((V["their_high"] == V["our_high"]).sum())

    cov = S.groupby("source")["lipid_class"].apply(lambda s: ", ".join(sorted(set(s)))).to_dict()
    (out / "README.md").write_text(README.format(
        agree=agree, n=len(V),
        classes="\n".join(f"- {k}: {v}" for k, v in cov.items())))
    print(f"species_sources: {len(S)} rows ({S.groupby('source').size().to_dict()})")
    print(f"quant: {len(q)} species in {q['lipid_class'].nunique()} classes")
    print(f"maldi peaks: {len(M)}; validation: {agree}/{len(V)} same call as theirs "
          f"(theirs high {int(V['their_high'].sum())}, ours high {int(V['our_high'].sum())})")
    for tab, names in failed.items():
        print(f"  {tab}: {len(names)} names not converted, e.g. {names[:6]}")


if __name__ == "__main__":
    main()
