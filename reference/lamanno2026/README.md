# Lipid Brain Atlas reference (mouse brain lipids with LC-MS/MS evidence)

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
their own MALDI peaks: 721 of 824 names with one candidate get the same call (score > 4
or not) as theirs. Differences are expected where their METASPACE term or their own MALDI
data (both left out here) tipped the balance.

Classes per source:

- fitzner: CE, CL, Cer, DG, HexCer, LPA, LPC, LPE, LPG, LPI, LPS, PA, PC, PE, PG, PI, PS, SM, TG
- lcms: BMP, CE, CL, Cer, CerP, DG, DGCC, GA1, Hex2Cer, HexCer, LPC, LPE, LPI, PA, PC, PE, PE-Cer, PG, PGP, PI, PI-Cer, PIP, PS, SHexCer, SM, TG
- msms: CE, CL, Cer, DG, FA, Hex2Cer, HexCer, LPC, LPE, LPG, LPI, LPS, PC, PE, PG, PI, PS, SM, TG
- vandenbosch: Cer, Hex2Cer, HexCer, LPC, LPE, PC, PE, SM
