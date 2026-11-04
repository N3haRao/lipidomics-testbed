# Replicating the analysis

Step-by-step instructions to rerun the whole analysis from the raw imzML files, either on a
Mac or Linux machine or on a SLURM cluster, and to check that a rerun matches. For what each
step does and why, see [`METHODS.md`](METHODS.md). For the Ada cluster specifically, with
every command including the ssh login, see [`ADA_RUN_GUIDE.md`](ADA_RUN_GUIDE.md).

Contents

1. [What you need](#1-what-you-need)
2. [Get the code](#2-get-the-code)
3. [Create the environment](#3-create-the-environment)
4. [Check the installation on synthetic data](#4-check-the-installation-on-synthetic-data)
5. [Put the data in place and point the config at it](#5-put-the-data-in-place-and-point-the-config-at-it)
6. [Run the pipeline](#6-run-the-pipeline)
7. [Check the outputs](#7-check-the-outputs)
8. [Reproduce a previous run exactly](#8-reproduce-a-previous-run-exactly)
9. [Change a decision and rerun only what is needed](#9-change-a-decision-and-rerun-only-what-is-needed)
10. [Guide to the result tables](#10-guide-to-the-result-tables)

---

## 1. What you need

| | |
|---|---|
| Data | the folder `imzML files for biostats/` from the shared drive (36 files: 6 slides × 3 modes, one `.imzML` and one `.ibd` each, about 40 GB). It is **not** in this repository: it is unpublished collaborator data. |
| Sample sheet | `metadata/sample_sheet.csv` (in the repository) |
| Software | [Miniforge](https://github.com/conda-forge/miniforge) (conda/mamba). Everything else is installed from `environment.yml` / `requirements.txt`. |
| Internet | once, for phase 00 to download LIPID MAPS, HMDB and the Allen atlas. After that the pipeline runs offline. |
| Machine | a laptop is enough for testing. For the full dataset, the cluster is faster: each phase job requests 4 cores and 32 GB (`bin/hpc/job_phase.sbatch`), and the three modes run in parallel. Peak memory is written to every log file. |
| Disk | about 40 GB for the data plus some tens of GB for outputs and caches |

The data folder must keep this layout (folder names with spaces are fine):

```
<data_root>/
  imzML files for biostats/
    Lipid Negative/       slide_e_lip_neg.imzML, slide_e_lip_neg.ibd, ... slide_k_...
    Lipid Positive/       slide_e_lip_pos.imzML, ...
    Metabolite Negative/  slide_e_met_neg.imzML, ...
```

The pipeline never writes into the data folder; the config loader refuses any output path
inside it.

## 2. Get the code

With git:

```bash
git clone https://github.com/N3haRao/lipidomics-testbed.git
cd lipidomics-testbed
git checkout <branch or commit to reproduce>
```

Or without git: on the GitHub page pick the branch, then Code > Download ZIP, and unzip. Both
ways work. Every run records a **code fingerprint** (a hash of the code itself, Section 8), so
a ZIP copy can be matched to a run just as well as a git commit.

## 3. Create the environment

```bash
bash bin/setup_env.sh
conda activate cln3-msi
```

This creates the `cln3-msi` environment (Python 3.11, versions pinned in `requirements.txt`),
installs this repository in editable mode, and registers a Jupyter kernel. Equivalent manual
route:

```bash
mamba env create -f environment.yml
conda activate cln3-msi
pip install -e .
```

## 4. Check the installation on synthetic data

```bash
bash bin/run_synthetic_test.sh
```

This generates a small fake dataset with known effects planted in it (same folder and file
layout as the real data, 6 slides × 3 modes, about 600 MB), runs every phase on it with
`config.test.yaml`, runs the unit tests, and checks that the planted effects are found.
It takes about 5 to 10 minutes and needs no internet. The last lines should read:

```
31 passed
...
27/27 checks passed
```

The checks cover recovery of the per-slide mass offsets and the planted position gradient,
correct Control/Mutant assignment, annotation of the planted lipids, registration quality,
the planted genotype effects (single peaks and class totals), that the slide consistency
check flags no slide when labels are correct, and that every phase recorded its provenance
with one code version and one config. A unit test also plants a swapped slide and
confirms that only that slide is flagged.

## 5. Put the data in place and point the config at it

Open `config.yaml` and set the `paths` block:

```yaml
paths:
  data_root: "/path/to/BRAIN - PILAR - SPATIAL"   # the folder that contains "imzML files for biostats"
  output_root: "/path/to/outputs"                 # anywhere outside data_root
```

Nothing else needs changing to reproduce the analysis as run. All other settings are the
decisions described in `METHODS.md`. For a cluster, `config.hpc.yaml` inherits `config.yaml`
and only overrides the paths, using `${CLN3_DATA}` for the data and `${CLN3_WORK}` for the
outputs (set from `HPC_DATA` and `HPC_WORK` in `bin/hpc/settings.sh`).

On Ada, the data are at `/home/nehar/cln3_msi/data/imzML files for biostats/`. Another user
can run the pipeline against that folder from their own account (read access needed, outputs
in their own home); see "Running as another user" in [`ADA_RUN_GUIDE.md`](ADA_RUN_GUIDE.md).

## 6. Run the pipeline

### On a laptop or workstation

```bash
conda activate cln3-msi
python scripts/00_setup.py --config config.yaml --download   # once: validate, download and cache references
bash bin/run_all.sh config.yaml                              # phases 01 to 10 in order, stops at the first error
```

`FROM=05 bash bin/run_all.sh config.yaml` restarts from a phase. Single phases:
`python scripts/<NN_name>.py --config config.yaml [--mode lip_neg]`.

If a download in phase 00 is blocked, download the files in a browser and pass them:
`python scripts/00_setup.py --config config.yaml --lipidmaps-file LMSD.sdf.zip --hmdb-file structures.zip`.

### On a SLURM cluster

From the Mac (or any machine with the data), after editing `bin/hpc/settings.sh` (login,
host, work folder, partition):

```bash
bash bin/hpc/upload_to_hpc.sh      # code + data, resumable
bash bin/hpc/start_on_hpc.sh       # environment, unit tests, phase 00, then submits all jobs
bash bin/hpc/status.sh             # progress
bash bin/hpc/fetch_results.sh      # copy summaries, figures and the report back
```

Jobs are chained with SLURM dependencies: 01 then 02, then for each mode 03 → 04 → 05 → 06 →
07 → 08 → 09 (the lipid negative chain first, because the other modes reuse its registration),
then 10. `FROM=09 bash bin/hpc/start_on_hpc.sh` resubmits from a phase.

### Phases

| phase | does | reads | writes (under `output_root`) |
|---|---|---|---|
| 00_setup | validates config and sample sheet; caches LIPID MAPS, HMDB, Allen atlas | config, sample sheet | `cache/` |
| 01_inventory | imzML headers and data type checks | imzML headers | `metadata/files.csv` |
| 02_qc_one_slide | detailed QC of one slide (for tuning) | one slide | `results/02_qc/` |
| 03_qc_all | QC of every slide, mass accuracy, section split, position gradient | all slides | `results/03_qc_all/<mode>/` |
| 04_peaks_matrix | peak list, recalibration, pixels × peaks matrix | 03, all slides | `processed/<mode>/matrix_raw.h5ad` |
| 05_annotate | sum-composition annotation and flags | 04 | `processed/<mode>/annotations.csv` |
| 06_normalize | TIC/RMS/median, slide-effect comparison | 04, 05 | `processed/<mode>/matrix_norm.h5ad` |
| 07_register | Allen atlas registration, region masks | 06 | `processed/<mode>/regions.parquet` |
| 08_explore | Moran's I, NMF programs, territories | 06, 07 | `processed/<mode>/territories.parquet` |
| 09_differential | per-animal paired tests, class tests, slide check | 06, 07, 08 | `results/09_differential/<mode>/` |
| 10_report | one self-contained HTML report | everything | `report/report.html` |

## 7. Check the outputs

* `report/report.html` opens offline in a browser. Read the **Slide consistency check** and
  **Reproducibility** sections near the top first.
* `summaries/<phase>_<mode>.md` is a short text summary of each phase, with warnings.
* `logs/` has one log per phase run, including memory use, the code fingerprint and the config
  hash.

Key quality numbers to compare with a previous run: Dice per section
(`results/07_register/<mode>/registration.csv`), the slide partial R² per normalization method
(`results/06_normalize/<mode>/method_comparison.csv`), the per-slide gradient
(`results/03_qc_all/<mode>/gradient.csv`), and the number of peaks kept and annotated
(`summaries/04_peaks_matrix_<mode>.md`, `summaries/05_annotate_<mode>.md`).

## 8. Reproduce a previous run exactly

Every phase writes `metadata/provenance/<phase>_<mode>.json` and the fully merged config it
used (`<phase>_<mode>_config.yaml`). The JSON has:

| field | meaning |
|---|---|
| `code_fingerprint` | hash of every file in `src/`, `scripts/`, `reference/` and `bin/`. Same fingerprint = same code, whether from git or a ZIP. |
| `git_commit` | the commit, if run from a git checkout ("dirty" if code was edited) |
| `config_sha256` | hash of the merged config |
| `sample_sheet_sha256` | hash of the sample sheet |
| `python`, `platform`, `packages` | interpreter and package versions |

To reproduce a run: use code with the same `code_fingerprint`, run with the saved config
(change only `paths`), and the same package versions. Random steps use `runtime.seed`, so the
numbers repeat. The report's Reproducibility section warns if the phases of one output folder
ran with different code or configs.

Reference databases are downloaded once and cached under `cache/`. Keep that folder with the
outputs if you need to repeat the annotation with exactly the same database versions; a fresh
download of LIPID MAPS or HMDB at a later date may differ slightly.

## 9. Change a decision and rerun only what is needed

| change in `config.yaml` | rerun from |
|---|---|
| `qc.sections.flip_slides`, `top_is_low_index`, `manual_split` (which section is which animal) | 03 |
| `qc.tissue.*` (tissue mask, edge trim) | 03 |
| `peaks.*` (peak picking, cap) | 04 |
| `annotation.*`, `hypotheses` | 05 |
| `normalization.*` | 06 |
| `registration.*` (including the list of atlas regions) | 07 |
| `explore.*` (Moran's I, NMF, territories) | 08 |
| `differential.*`, `class_tests.*`, `pair_check.*` | 09 |
| `report.*` | 10 |

Example, if the sectioning records confirm that slide J's sections are swapped:

```yaml
qc:
  sections:
    flip_slides: ["J"]
```

then `FROM=03 bash bin/run_all.sh config.yaml` (or `FROM=03 bash bin/hpc/start_on_hpc.sh`).

To add atlas regions, list them under `registration.regions` with Allen acronyms (all
descendants are included), for example `striatum: ["STR"]`, and rerun from 07.

## 10. Guide to the result tables

### `results/09_differential/<mode>/differential.csv` (one row per peak and region)

| column | meaning |
|---|---|
| `unit_kind`, `unit` | region (atlas region or whole section) or territory |
| `peak`, `mz`, `label`, `annotation`, `adduct`, `lipid_class` | the peak and its annotation |
| `hypothesis`, `targeted` | hypothesis tag; targeted peaks get their own FDR family |
| `n_pairs` | slides with both animals measured |
| `log2fc`, `ci_low`, `ci_high` | mean log2(Mutant/Control) over slides, 95% CI |
| `p_ttest` | paired t-test (primary) |
| `p_wilcoxon`, `p_sign` | exact Wilcoxon signed-rank and sign (binomial) tests, robustness checks |
| `n_pairs_up` | pairs with a positive difference |
| `q_targeted`, `q_all` | Benjamini-Hochberg q within region: targeted family, all peaks |
| `significant` | q < `fdr_alpha` in the peak's family |
| `matches_hypothesis` | direction agrees with the hypothesis tag |
| `could_be_position` | smaller than 2 × the position gradient, same sign (warning) |
| `log2fc_Female`, `log2fc_Male`, `p_sex_difference_exploratory` | effect per sex and the exploratory sex × genotype test |
| `ambiguous`, `typeII_flag`, `candidates` | annotation caveats |

### `results/09_differential/<mode>/class_tests.csv` (one row per class or ratio and region)

The same test columns as above, plus `feature`, `kind` (class or ratio), `n_species`,
`pct_change`, `q_bh` (one family per mode) and `pair_<slide>` (each slide's own log2
difference). If a slide was flagged by the consistency check, extra columns
`without_<slide>_*` and `flipped_<slide>_*` give the exploratory sensitivity analysis.
`class_members.csv` lists the species in each class, and `class_values_per_animal.csv` the
per-animal values that were tested.

### `results/09_differential/<mode>/pair_consistency.csv`

Each slide's median Spearman correlation with the other slides' difference profiles, per
region, and whether it was flagged. The full correlation matrices are in
`pair_correlation_<region>.csv`.

### `results/09_differential/<mode>/per_animal_means.csv`

The per-animal, per-region mean normalized intensity of every tested peak, with the number of
pixels. Everything in 09 is computed from these values, so they are the starting point for any
re-analysis with different statistics.
