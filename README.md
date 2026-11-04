# CLN3 KO brain MALDI-MSI lipidomics pipeline

Analysis code for Pilar Andrade's (Boeynaems lab, NRI) MALDI-MSI data of the CLN3 knockout
mouse brain: which lipids accumulate or get depleted in the KO, and where, with a focus on
demyelination.

The structure follows my CAJAL tutorial repo (`spatial-lipidomics-project`, notebooks 01-09):
matrix building, annotation, normalization, anatomy, Moran's I, NMF, territories, differential
testing. The difference is that every step here is a script with a command-line interface, so
the same code runs on the laptop and (later) as a SLURM job. A few notebooks are there just to
look at the figures.

> Unpublished collaborator data. Nothing in this repo uploads anything. The only network access
> is the one-time download in `00_setup.py --download` (LIPID MAPS, HMDB, Allen atlas). The
> imaging data are not in this repository.

## Documentation

| document | for |
|---|---|
| [`docs/METHODS.md`](docs/METHODS.md) | study design, every processing step, and which statistical tests are used and why |
| [`docs/REPLICATION.md`](docs/REPLICATION.md) | rerunning the analysis from the raw files, checking a rerun matches, which phase to rerun after a config change, and a guide to the result tables |
| [`docs/ADA_RUN_GUIDE.md`](docs/ADA_RUN_GUIDE.md) | the same on the Ada cluster, every command from the ssh login on |
| this README | repository layout, configuration keys, phase overview, troubleshooting |

## Quick start

```bash
bash bin/setup_env.sh && conda activate cln3-msi     # environment (once)
bash bin/run_synthetic_test.sh                        # check the install on fake data: 27/27 checks
# set paths.data_root and paths.output_root in config.yaml, then:
python scripts/00_setup.py --config config.yaml --download
bash bin/run_all.sh config.yaml                       # phases 01 to 10
open "<output_root>/report/report.html"
```

Every run records the code fingerprint, config hash and package versions in
`<output_root>/metadata/provenance/`, and the report shows them, so a result can always be
traced to the exact code and settings that produced it.

---

## Repository layout

```
config.yaml              every path and parameter (Mac defaults); the only file you normally edit
config.test.yaml         the synthetic test config (inherits config.yaml, overrides a few keys)
environment.yml          conda env, Python 3.11 (pins in requirements.txt)
requirements-optional.txt   harmonypy, jupyterlab (uMAIA: separate env, see below)
metadata/sample_sheet.csv   the 12 animals + the imzML file name for each mode
reference/               lipid class rules, adducts, calibrants, curated metabolites (all as formulas)
src/cln3msi/             the shared code the scripts import
scripts/00_setup.py ... 10_report.py, 04b_scils_crosscheck.py, make_synthetic_data.py
notebooks/               review notebooks (they only read outputs, they never recompute)
bin/                     setup_env.sh, run_phase.sh, run_all.sh, run_synthetic_test.sh
tests/                   unit tests + the synthetic end-to-end check
```

Outputs all go under `paths.output_root` (never inside the data folder; the config loader refuses):

```
<output_root>/
  metadata/    files.csv (01), header_<file>.txt (every header cvParam)
  processed/<mode>/  qc/*.npz (03), matrix_raw.h5ad + peaks.csv (04), annotations.csv (05),
                     matrix_norm.h5ad (06), regions.parquet (07), territories.parquet (08)
  results/<phase>/<mode>/   tables (csv/json)
  figures/<phase>/<mode>/   png
  summaries/<phase>_<mode>.md   <- short text summary of each phase, with warnings
  metadata/provenance/          code fingerprint, config and package versions of every phase run
  logs/<phase>_<mode>_<time>.log  (includes memory use)
  cache/       LIPID MAPS, HMDB, Allen atlas
  report/report.html
```

---

## 1. Setup (Mac, once)

Needs [Miniforge](https://github.com/conda-forge/miniforge) (gives you `mamba`/`conda`).

```bash
cd ~/path/to/lipidomics-testbed
bash bin/setup_env.sh
conda activate cln3-msi
```

That creates the `cln3-msi` env (Python 3.11, pinned packages), installs this repo with
`pip install -e .`, tries the optional extras (a failure there is fine), and registers a
Jupyter kernel called `cln3-msi`.

**Optional: uMAIA.** It needs JAX and numpy < 2, which clash with this env (same as in the
CAJAL tutorial). If you want it, make a separate env following the uMAIA instructions
(https://github.com/lamanno-epfl/uMAIA), run `pip install -e .` from this repo inside it, set
`normalization.umaia.enabled: true`, and run phase 06 in that env. If uMAIA can't be imported,
06 says so in the log and carries on with TIC/RMS/median only. (The uMAIA route has not been
tested; the analysis as run does not use it.)

## 2. Test on fake data first (no real data or internet needed)

```bash
conda activate cln3-msi
bash bin/run_synthetic_test.sh
```

This writes ~600 MB of fake imzML/ibd files to `synthetic_data/` (same folder names with spaces,
same file names, 6 slides x 3 modes; each mode uses a different storage/spectrum type so every
code branch runs), runs phases 00 to 10 with `config.test.yaml`, runs the unit tests, and checks
that the planted effects were found (27 checks, e.g. "Cer up in KO", "sulfatide down in white
matter", "slide G has the strongest gradient"). Takes ~5 minutes. `FRESH=1 bash bin/run_synthetic_test.sh`
regenerates the fake data.

## 3. Configure for the real data

Open `config.yaml`. The things to check:

| key | what |
|---|---|
| `paths.data_root` | `"/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL"` (keep the quotes and the upper case, there are spaces) |
| `paths.output_root` | where outputs go. Must not be inside `data_root` |
| `modes.run` | which modes to process: `[lip_neg, lip_pos, met_neg]` |
| `qc.slide`, `qc.mode` | the slide/mode phase 02 looks at first (default E, lip_neg) |
| `qc.sections.*` | split direction (`split_axis`: auto, rows or cols), which side is the sheet's Top (`top_is_low_index`, `flip_slides`), or a forced split line per slide (`manual_split`) |
| `annotation.ppm_tolerance` | leave `null` to use the value measured in 02/03 |
| `normalization.for_stats` | `tic` (decided up front; don't change it after seeing 09) |
| `registration.optical_tiffs` | per-slide TIFF path, relative to `data_root`, once Pilar sends them |
| `scils_crosscheck.enabled` | `true` to compare with the SCiLS exports |
| `class_tests.*` | which lipid classes are summed and which ratios are tested, per mode (set before looking at results) |
| `pair_check.*` | slide consistency check: regions used, flag threshold, and whether to add exploratory left-out / flipped columns for a flagged slide |

### What's in the data folder, and what the pipeline uses

Only `imzML files for biostats/` is read (36 files: 6 slides x 3 modes, `.imzML` + `.ibd` each,
all present). Everything else is ignored by the preliminary analysis:

| item | what it is | used? |
|---|---|---|
| `imzML files for biostats/` | SCiLS export of slides E-K, the main input | **yes** |
| `raw lip neg/`, `raw lip pos/`, `raw met neg/` | one SCiLS project per slide (`.slx` + `.sbd`); the optical scans live in here | no (SCiLS-only format). Source for the TIFF exports, if we ask |
| `*.slx/*.sbd` at the top and in `Female Data - Pilar analysis/` | Pilar's combined SCiLS projects | no |
| `*.srd` (e.g. `CC-Ctrl-All.srd`) | SCiLS region definitions (hand-drawn regions such as corpus callosum) | no; could serve as manual ROIs later if exported as spot lists |
| `*.mca` | SCiLS feature (m/z interval) lists | no |
| `export_All_*.csv/.xlsx` | SCiLS intensity exports | only by the optional 04b cross-check (off) |
| `Results_CLN3_Analysis*/`, `PC 38.1 - depleted.png`, `.Rhistory`, `.RDataTmp` | Pilar's results and R session | no |
| `09SEP26/` | a later, separate **cerebellum** slide Z (raw Bruker `.d` folders, method files, its own TIFF) | no. Only tells us the instrument is probably a Bruker timsTOF (fleX) |

## 4. Run

Every phase is `python scripts/<phase>.py --config config.yaml [--mode lip_neg]`. Without
`--mode`, it loops over `modes.run`. The same thing with the helper: `bash bin/run_phase.sh 04 config.yaml --mode lip_neg`.

```bash
conda activate cln3-msi
cd ~/path/to/lipidomics-testbed

python scripts/00_setup.py --config config.yaml --download     # once: validate + download/caches
python scripts/01_inventory.py --config config.yaml            # headers only, a few minutes
python scripts/02_qc_one_slide.py --config config.yaml         # slide E, lip_neg; look at the figures!
python scripts/03_qc_all.py --config config.yaml               # all slides: the slow one
python scripts/04_peaks_matrix.py --config config.yaml
python scripts/04b_scils_crosscheck.py --config config.yaml    # optional
python scripts/05_annotate.py --config config.yaml
python scripts/06_normalize.py --config config.yaml
python scripts/07_register.py --config config.yaml
python scripts/08_explore.py --config config.yaml
python scripts/09_differential.py --config config.yaml
python scripts/10_report.py --config config.yaml
```

Or everything in order, stopping at the first error:

```bash
bash bin/run_all.sh config.yaml
FROM=05 bash bin/run_all.sh config.yaml     # restart from a phase
```

If a download in 00 is blocked (LIPID MAPS sometimes wants a click-through), download the file
in the browser and point to it:

```bash
python scripts/00_setup.py --config config.yaml \
  --lipidmaps-file "$HOME/Downloads/LMSD.sdf.zip" --hmdb-file "$HOME/Downloads/structures.zip"
```

**Suggested first pass on new data:** run 00, 01 and 02, and check `summaries/01_inventory_*.md`,
`summaries/02_qc_lip_neg.md` and the `figures/02_qc/lip_neg/*.png` (especially
`sections_E.png` and `mean_spectrum_E.png`). Then we tune before the long phases.

---

## What each phase does and writes

| phase | does | main outputs |
|---|---|---|
| **00_setup** | validates config + sample sheet (one Control and one Mutant per slide, Top/Bottom, file names), checks the data folder, caches LIPID MAPS / HMDB as parquet, downloads the Allen CCF via brainglobe, lists optional packages | `cache/`, `summaries/00_setup.md` |
| **01_inventory** | reads every imzML header with lxml (no intensities): spectrum type, storage, polarity, instrument/software cvParams, pixel grid, #spectra, m/z range, data types. Checks profile vs centroid **from the data** (30 spectra) and flags mismatches | `metadata/files.csv`, `metadata/header_*.txt`, acquired-pixel maps |
| **02_qc_one_slide** | 2 streaming passes over one slide: TIC / non-matrix TIC / DAN images, tissue mask, Top/Bottom split + labelled figure, mean spectra (Top, Bottom, off-tissue), calibrant mass error (m/z from formulas), DAN ions, resolving power from single profile spectra (fit says TOF-, Orbitrap- or FT-ICR-like), gradient check, suggested ppm tolerance | `results/02_qc/<mode>/`, `figures/02_qc/<mode>/` |
| **03_qc_all** | 02 on every slide + cross-slide mean spectra, per-slide m/z offsets, pooled tolerance, section areas, the top-vs-bottom gradient table | `results/03_qc_all/<mode>/{mass_accuracy.json,gradient.csv,sections.csv}`, `sections_all.png` |
| **04_peaks_matrix** | peak picking on each section's mean spectrum (= centroiding if profile), per-slide recalibration, alignment across the 12 sections, keep peaks in >= 50% of sections, integrate each peak in every tissue pixel, presence filter, DAN peaks moved out of X | `processed/<mode>/matrix_raw.h5ad`, `peaks.csv` |
| **04b_scils_crosscheck** | inspects a SCiLS CSV/XLSX export (prints what it finds), compares its m/z list to ours | `results/04b_scils_crosscheck/<mode>/` |
| **05_annotate** | sum-composition library generated from class rules (exact formulas), filtered by LIPID MAPS; HMDB + curated list for metabolites; mode adducts; flags ambiguous, isotope (M+1/M+2) and type II overlaps; hypothesis tags | `processed/<mode>/annotations.csv`, `hypothesis_peaks.csv` |
| **06_normalize** | TIC, RMS, median (and uMAIA if available); slide effect before/after (partial R2 of slide on section means, PCA); X = TIC for stats; any batch-corrected embedding only in obsm for plots | `processed/<mode>/matrix_norm.h5ad`, `slide_effect_pca.png` |
| **07_register** | per section: guide image (TIFF, else sulfatide/HexCer image, else PC1), AP search around Bregma -1.58 to -2.18, similarity -> affine -> B-spline (SimpleITK, mutual information), warped Allen labels, region masks, Dice QC, lipid WM/GM fallback | `processed/<mode>/regions.parquet`, `overlays.png`, `registration.csv` |
| **08_explore** | Moran's I per section with permutations, seeded NMF on Control, Leiden territories on Control transferred to all pixels by kNN | `territories.parquet`, `morans_i.csv`, `nmf_programs.csv` |
| **09_differential** | per-animal means per region / territory, log2, paired over slides: paired t (primary, with 95% CI), Wilcoxon; BH within region (targeted and all); position flag; exploratory sex. Then, from the same means: **class totals and ratios** (`class_tests`, one BH family per mode) and the **slide consistency check** (`pair_check`: a slide whose KO/Control difference is the mirror image of the others is flagged as a possible label swap) | `results/09_differential/<mode>/{differential,class_tests,pair_consistency}.csv`, volcanoes, forest plots, per-slide dot plot, correlation heatmap |
| **10_report** | self-contained HTML with QC, overlays, tables, top-hit ion images (Control vs KO per slide, shared scale), volcanoes, caveats, all summaries | `report/report.html` |

### How the "don't assume" rules are handled

* **Profile vs centroid, continuous vs processed:** detected per file in 01 (header *and* data),
  used everywhere after. Continuous uncompressed files are read through a read-only numpy
  memmap; processed files spectrum by spectrum. Processed data are binned on a ppm grid
  (`qc.bin_ppm`, never coarser than the data's own spacing for profile data).
* **Instrument / resolution:** whatever the header says is logged; resolving power is measured.
* **Mass accuracy:** measured on calibrants (formulas in `reference/calibrants.yaml`); 04
  recalibrates each slide by its median offset; 05 takes its tolerance from that.
* **Memory:** nothing loads a whole slide. Peak RAM is logged in every log file. The biggest
  object is the pixels x peaks matrix in 04/06 (e.g. 250k pixels x 1000 peaks = 1 GB).

### Design limits the code respects

1. **Paired by slide**, n = 6 pairs; slide is the pairing factor.
2. **Genotype = position** (Control always on top). 03 measures the top-to-bottom trend in the
   off-tissue DAN signal and TIC; 09 flags hits within `qc.gradient.flag_ratio` x that trend
   (`could_be_position`). Not corrected, only reported.
3. **Sex = run** (females runs 1-3, males run 4). Sex results are labelled exploratory.
4. **Animal = unit.** All tests are on per-animal means, never on pixels.
5. **Labels are checked, never changed by the code.** 09 compares the slides' KO/Control
   difference profiles; a slide that runs opposite to the others is flagged as a possible swap.
   Only the sectioning/mounting records can confirm it; then set `qc.sections.flip_slides`.

---

## HPC (Ada, SLURM)

**Step-by-step guide with every command, including the ssh login: [`docs/ADA_RUN_GUIDE.md`](docs/ADA_RUN_GUIDE.md).**

Everything is driven from the Mac terminal by the scripts in `bin/hpc/`:

| script | runs on | does |
|---|---|---|
| `settings.sh` | (edit once) | your Ada login, work folder, partition/account, local data path |
| `upload_to_hpc.sh` | Mac | rsync the code + `imzML files for biostats` (~40 GB, resumable) |
| `start_on_hpc.sh` | Mac | ssh in, run `hpc_setup.sh`, then `submit_all.sh` |
| `hpc_setup.sh` | Ada login node | conda (module or Miniforge), env, unit tests, phase 00 + downloads |
| `submit_all.sh` | Ada login node | submits 01 → 02, and per mode 03 → 09 as dependent jobs, then 10 |
| `job_phase.sbatch` | Ada compute node | one phase (4 cores, 32 GB by default) |
| `status.sh` / `fetch_results.sh` | Mac | queue + progress / copy summaries, figures, report back |

```bash
cd ~/path/to/lipidomics-testbed
nano bin/hpc/settings.sh                 # already set for nehar@ada.nri.bcm.edu; change HPC_WORK / partition if needed
bash bin/hpc/upload_to_hpc.sh            # code + data
bash bin/hpc/start_on_hpc.sh             # setup + submit everything
bash bin/hpc/status.sh                   # any time
bash bin/hpc/fetch_results.sh            # when jobs are done -> ~/cln3_msi_outputs_from_ada
```

On Ada the layout is `HPC_WORK/{repo,data,outputs,envs,slurm_logs}`, and `config.hpc.yaml` points
at it through `${CLN3_WORK}` (outputs) and `${CLN3_DATA}` (data, `HPC_DATA` in `settings.sh`,
default `HPC_WORK/data`). A second user can read the data from another account's folder and
write to their own: see "Running as another user" in `docs/ADA_RUN_GUIDE.md`. Edit code and config on the Mac and re-run
`upload_to_hpc.sh --code` (the upload mirrors the repo, so edits made directly on Ada get
overwritten). To rerun from a phase: `FROM=05 bash bin/hpc/start_on_hpc.sh`.

## Troubleshooting

* `CONFIG / DATA ERROR` with an expected layout printed: `paths.data_root` is wrong or a mode
  folder is missing. Paths with spaces must be quoted in YAML.
* Wrong Top/Bottom: look at `figures/03_qc_all/<mode>/sections_all.png`, then use
  `qc.sections.top_is_low_index: false` (all slides) or `flip_slides: [J]`.
* Only one section found on a slide: set `qc.sections.manual_split: {E: <row or column>}`.
* Registration bad on some sections: check `overlays.png`; sections below `registration.min_dice`
  automatically get only the lipid WM/GM masks. A TIFF usually fixes it.
* Any other error: the newest log in `<output_root>/logs/` has the full traceback and the
  settings the phase used.

## Notebooks

`notebooks/02_qc_review.ipynb`, `07_registration_review.ipynb`, `08_territories_review.ipynb`,
`09_hits_review.ipynb`. Open with the `cln3-msi` kernel; set `CONFIG` in the first cell.
