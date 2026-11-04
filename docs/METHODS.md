# Methods: study design, processing and statistics

This document describes what the pipeline does, step by step, and why each statistical
choice was made. It contains no results. For how to run it, see
[`REPLICATION.md`](REPLICATION.md). Every parameter named here lives in `config.yaml`, and the
value used in a given run is saved in `<output_root>/metadata/provenance/`.

Contents

1. [Study design](#1-study-design)
2. [Data](#2-data)
3. [Processing workflow](#3-processing-workflow)
4. [Statistical analysis](#4-statistical-analysis)
5. [Quality checks for the design's weak points](#5-quality-checks-for-the-designs-weak-points)
6. [What the analysis cannot do](#6-what-the-analysis-cannot-do)
7. [Correspondence with the CAJAL tutorial](#7-correspondence-with-the-cajal-tutorial)
8. [Software](#8-software)

---

## 1. Study design

| | |
|---|---|
| Animals | 12 mice: 6 CLN3 knockout ("Mutant"), 6 controls |
| Tissue | one coronal brain section per animal, 10 µm, around Bregma −1.6 to −2.2 mm |
| Slides | 6 glass slides (E, F, G, I, J, K), each with one control and one knockout section |
| Position | the control is always the "Top" section, which in the exported images is the left one |
| Sex | slides E, F, G female; slides I, J, K male (both animals on a slide share sex) |
| Acquisition runs | E, F, G in runs 1, 2, 3 (one slide each); I, J, K together in run 4 |
| Ion modes | lipids negative (m/z 300 to 2,000), lipids positive (m/z 300 to 2,000), metabolites negative (m/z 50 to 750) |
| Matrix, pixel size | 1,5-diaminonaphthalene (DAN), 50 µm |

The sample sheet (`metadata/sample_sheet.csv`) lists every animal with its genotype, sex,
slide, position and run, and the imzML file for each mode. Phase 00 refuses to run if a slide
does not have exactly one control and one knockout, or if a position is missing.

**Unit of replication.** The animal. Each animal contributes one section, so there are 6
independent observations per genotype, arranged as **6 pairs** (one pair per slide). Pixels
are measurements within an animal, not replicates.

**Consequences of the layout.** Three features of the design shape the analysis:

* Two animals on one slide share the matrix coating, the instrument run and the day. This
  makes the slide a natural block, and comparing the two animals within a slide removes
  everything they share.
* Genotype is completely confounded with position on the slide (Section 5.1).
* Sex is completely confounded with acquisition run (Section 5.2).

## 2. Data

The input is the SCiLS Lab export in imzML format: one `.imzML` (XML metadata) and one `.ibd`
(binary spectra) per slide and mode, 36 files in total. The files are never modified;
the configuration loader refuses any write path inside the data folder.

Phase 01 reads every header without loading spectra and records the spectrum type, storage
mode, polarity, instrument parameters, pixel grid and m/z range. Whether spectra are profile
or centroid is also tested **from the data** (30 spectra per file), because headers can be
wrong. For this dataset the files are *processed* (each spectrum has its own m/z array) and
*centroided*. The header names a Bruker flex-series TOF instrument.

Spectra are read lazily, one batch at a time (`runtime.spectra_chunk`), so no slide is ever
held in memory in full.

## 3. Processing workflow

Each phase is a script (`scripts/NN_name.py`) that reads the outputs of earlier phases and
writes tables, figures and a short summary. Phases 03 to 09 run separately for each mode.

### 3.1 Inventory and quality control (phases 01 to 03)

* **Tissue detection.** A tissue mask is computed from the non-matrix total ion current (TIC
  minus the DAN matrix ions), smoothed, thresholded and cleaned (`qc.tissue`). The threshold
  method "auto" uses Otsu's method, but switches to the lower class boundary of a three-class
  Otsu split when that fills holes inside the sections (this happens when the imaged area
  barely extends past the tissue, so that Otsu separates white from grey matter instead of
  tissue from glass). One pixel at every tissue edge is excluded from statistics
  (`erode_edge_px`), because edge pixels mix tissue and glass.
* **Section split.** The two largest tissue pieces are assigned to the two animals. The
  direction of the split (sections above each other or side by side) is detected from the
  data (`split_axis: auto`); for this dataset the sections are side by side. The left
  section is the sample sheet's "Top" (`top_is_low_index: true`). A slide can be reversed
  individually with `flip_slides`, and a split line can be forced with `manual_split`.
* **Mass accuracy.** Abundant, well-known ions (for example PI 38:4, PE 40:6 and the main
  sulfatides in lipid negative mode, PC 34:1 adducts in positive mode, taurine and glutamate in
  the metabolite mode), defined by formula in `reference/calibrants.yaml`, are located on each
  slide's mean spectrum, and their ppm error is measured. The per-slide median error is used to recalibrate that slide in
  phase 04, and the spread of errors sets the annotation tolerance in phase 05.
* **Position gradient.** See Section 5.1.

### 3.2 Peak list and data matrix (phase 04)

1. A mean spectrum is computed for each of the 12 sections. Because the data are centroided
   and processed, the centroids are binned on a ppm grid (`qc.bin_ppm` = 2 ppm) before
   averaging.
2. Peaks are picked on each section's mean spectrum (signal-to-noise ≥ `peaks.snr` = 5, noise
   from a rolling median). The peak centre is refined by fitting a parabola to the log
   intensities of the top three points.
3. Each slide's peak positions are shifted by its own calibration offset, then peaks from
   all sections are grouped when they lie within the alignment tolerance.
4. A consensus peak is kept if it was picked in at least half of the sections
   (`min_section_fraction` = 0.5) and is present in at least 1% of tissue pixels
   (`min_pixel_fraction`). At most `max_peaks` = 2,000 peaks are kept, the strongest first.
5. Every kept peak is integrated in every tissue pixel (sum of intensities inside the peak
   window), giving a pixels × peaks matrix per mode. DAN matrix peaks are moved out of the
   matrix into separate columns, so they do not enter normalization or statistics.

### 3.3 Annotation (phase 05)

Lipids are annotated at **sum-composition level** (for example "HexCer 42:2;O2", meaning a
hexosylceramide with 42 carbons and 2 double bonds in total, and two oxygens on the sphingoid
base). A library of candidate compositions is generated from class rules in
`reference/lipid_classes.yaml` (an anchor formula plus CH2 and H2 increments), so every
candidate has an exact elemental formula. Compositions whose formula is absent from LIPID MAPS
are removed. For the metabolite mode, HMDB (or the curated list in
`reference/metabolites_curated.csv`) is added.

A peak matches a candidate if the ppm error of some (candidate, adduct) is within the
tolerance from phase 03. Adducts follow ion chemistry:

* negative mode: [M−H]⁻ by default; gangliosides also [M−2H]²⁻; choline lipids (PC, SM, LPC,
  LSM) only as [M+Cl]⁻ or [M−CH3]⁻, because they have no acidic proton. This rule matters:
  "SM 35:1 [M−H]⁻" has exactly the formula of SM 36:1 [M−CH3]⁻, so allowing [M−H]⁻ for SM
  would invent odd-chain sphingomyelins;
* positive mode: [M+H]⁺, [M+Na]⁺, [M+K]⁺; sulfatides are excluded (they are anions).

Every peak also gets flags: **ambiguous** (more than one composition fits), **isotope** (the
peak is the ¹³C M+1 of a lower-mass peak; isotope peaks are never tested), and **type II
overlap** (the M+2 isotope of a species with one more double bond could fall inside the
peak). Peaks of the classes the study is about get a **hypothesis tag** (`hypotheses` in the
config): accumulated (Cer, HexCer, Hex2Cer, gangliosides, LSM), depleted (SM, sulfatide), and
secondary (BMP/PG, glycerophosphodiesters, sterols).

### 3.4 Normalization (phase 06)

Each pixel is divided by a size factor: TIC (sum of all kept non-matrix peaks), RMS, or median.
The methods are compared by how much slide-to-slide variation remains: for every peak, the 12
section means are log-transformed, and the partial R² of slide (after genotype) is computed.
**TIC is used for all statistics, and this was fixed before any genotype result was
examined** (`normalization.for_stats`), so that the method cannot be chosen to favour a result.
Batch-corrected embeddings (for example Harmony) are used only for plots, never for tests.

### 3.5 Registration to the Allen atlas (phase 07)

Each section is aligned to the Allen Mouse Brain Common Coordinate Framework (25 µm, via
brainglobe) with SimpleITK:

1. **Guide image.** The strongest annotated sulfatide ion, which outlines white matter
   sharply. An optical scan can be used instead when one is available.
2. **Search.** Coronal atlas planes are tried across the expected anteroposterior range in
   0.1 mm steps, each from four starting rotations (0°, 90°, 180°, 270°).
3. **Fit.** Similarity, then affine transform, with Mattes mutual information as the metric.
   The best plane and rotation (by normalized cross-correlation) is then refined with a
   B-spline (1 mm control-point grid). The B-spline result is kept only if it does not reduce
   the overlap.
4. **Quality.** The Dice coefficient between the tissue mask and the warped atlas outline
   (0 = no overlap, 1 = perfect). Sections below `min_dice` = 0.8 do not get atlas regions.

Warped atlas labels give each pixel an Allen structure. The regions used are white matter
("fiber tracts"), corpus callosum, isocortex, hippocampal formation and thalamus, each with all
descendant structures. The three modes were exported on the same pixel grid, so the
registration is done once (lipid negative mode, which has the sulfatide guide) and the same
regions are applied to the other modes after checking that the pixels match.

### 3.6 Spatial exploration (phase 08)

This phase describes spatial structure. It does not test genotype.

* **Moran's I** per peak and section, on a 6-nearest-neighbour graph of pixel positions,
  with a permutation p-value (199 permutations). A peak counts as spatially structured if
  I ≥ 0.3 in at least half of the sections.
* **Lipid programs.** Non-negative matrix factorization (NMF) of the structured peaks,
  trained on control pixels only, with the number of programs chosen at the elbow of the
  reconstruction-error curve. Before NMF, each peak is divided by its per-slide median so
  that the programs follow anatomy rather than slide-to-slide intensity differences
  (exploration only).
* **Lipid territories.** Leiden clustering of the program activities of a control
  subsample, then transferred to every pixel (remaining controls and all knockouts) by a
  k-nearest-neighbour vote. Training on controls only means the knockout cannot shape the
  territories. Territories are used as additional regions in phase 09.

## 4. Statistical analysis

### 4.1 Per-animal summaries

For every region (whole section, each atlas region, each territory) and every peak:

1. take the TIC-normalized intensity of the animal's pixels in that region, excluding edge
   pixels;
2. average them, giving **one value per animal** (at least `min_pixels_per_unit` = 20 pixels
   are required, otherwise the value is missing);
3. take log2 (a pseudocount of half the smallest positive mean is added only for peaks where
   some animal's mean is exactly 0).

Averaging first is what makes the animal the unit. Testing pixels directly would treat
thousands of measurements from one brain as independent and give meaninglessly small p-values
(pseudoreplication).

### 4.2 The primary test: paired t-test

For each slide *s*, the difference d<sub>s</sub> = log2(Mutant<sub>s</sub>) − log2(Control<sub>s</sub>).
With 6 slides there are 6 differences. The primary test is a one-sample t-test of the mean
difference against 0 (equivalently, a paired t-test), with 5 degrees of freedom.

* **log2 fold change** = mean of the d<sub>s</sub>. 0 means no change, +1 doubling, −1 halving;
  −0.5 is about 29% lower.
* **95% confidence interval** = mean ± t<sub>0.975, n−1</sub> × SD / √n.
* A test is run only if at least `min_pairs` = 5 slides have both animals measured.

**Why paired.** The slide is a block: the two animals on a slide share the coating, the run
and the day, and in this dataset the slide explains most of the variance between sections
(phase 06 reports the partial R²). Subtracting within the slide removes that shared component,
so the test compares genotypes against the within-slide noise only. An unpaired test (6 vs 6)
would include the slide variance in its error term and lose power. The paired t-test is
mathematically identical to the linear model `log2 ~ genotype + slide` with slide as a fixed
block. When the slide effect happens to be small, an unpaired test can give a slightly lower
p-value, but the paired test is the one that matches how the samples were processed.

**Why a t-test on log values.** Intensities are positive and multiplicative (a lipid that is
"20% higher" is higher by a ratio), so log differences are roughly symmetric and the t-test's
assumptions are reasonable. The t-test uses the size of the changes, not only their sign,
which matters with 6 pairs.

**Why not a mixed model.** The model `log2 ~ genotype + (1 | slide)` is the textbook
alternative. With only two animals per slide, its slide-variance estimate collapsed and it
reported implausibly small p-values (far below the paired t-test on the same numbers). It is
off by default (`run_mixedlm: false`).

### 4.3 Robustness checks reported next to every test

* **Exact Wilcoxon signed-rank test** on the 6 differences. It ranks the differences, so it
  does not assume normality and shows whether one extreme pair drives a t-test result. With
  6 pairs, the smallest possible two-sided p-value is 2/64 = 0.031, so it cannot survive the
  multiple-testing correction on its own; it is a check, not a test for calling hits.
* **Sign test (binomial)**: how many of the pairs go up, compared with 50:50. It ignores the
  size of the changes and has the same floor of 0.031.
* **Number of pairs in the direction of the mean** (`n_pairs_up`) and every pair's own
  difference (`pair_<slide>` columns in the class tests).

### 4.4 Multiple testing

p-values are corrected with the Benjamini-Hochberg false discovery rate procedure. The
corrected value is the q-value, and q < 0.05 (`fdr_alpha`) is called significant. Families,
all set before looking at results:

* **single peaks, targeted**: within each region, only the peaks with a hypothesis tag;
* **single peaks, untargeted**: within each region, every tested peak. Tested peaks are those
  that are annotated or spatially structured (`untargeted_peaks: annotated_or_structured`), so
  that hundreds of unstructured noise peaks do not dilute the correction;
* **class-level and ratio tests**: one family per mode (Section 4.5).

### 4.5 Class-level and ratio tests

Single lipids are numerous and noisy. The class tests add up the species of one class per
animal and region (for example all sulfatides, all gangliosides, the three
glycerophosphodiesters), and also test ratios within the same animal and region (for example
sulfatide / (PE+PI+PS), GPI / LPI, HexCer / Cer). The paired t-test, the robustness checks and a
single BH family per mode are applied as above.

Ratios have a specific advantage here: two lipids measured in the same pixels share the
pixel's matrix coating, ionization conditions and normalization factor, so those cancel in the
ratio. The list of classes and ratios is in `config.yaml` (`class_tests`) and was fixed before
the genotype results were examined. Species with ambiguous class assignments, isotope peaks and
species detected in fewer than `min_sections` = 10 sections are excluded; the members of every
class are listed in `class_members.csv`.

### 4.6 Sex

Both animals on a slide are the same sex, so **sex cancels in the within-slide difference**:
the paired analysis is already adjusted for it, and sex cannot be added as a separate
covariate. What can be asked is whether the genotype effect differs between the sexes (the sex
× genotype interaction). This is tested as a Welch t-test of the 3 female pair differences
against the 3 male ones (`p_sex_difference_exploratory`, with `log2fc_Female` and
`log2fc_Male`). Because all females were run in runs 1 to 3 and all males in run 4, this is
equally a run × genotype comparison, and it is reported as exploratory only.

## 5. Quality checks for the design's weak points

### 5.1 Genotype and position: the gradient check

The control is always on the left. A technical trend across the slide (uneven matrix coating,
instrument drift during the raster) would therefore look exactly like a genotype effect, and no
statistical model can remove it from this design. Two measures limit the risk:

1. **Gradient check (phase 03).** Off-tissue pixels contain matrix but no tissue, so any trend
   in their signal along the split axis is purely technical. A robust (Theil-Sen) slope of the
   log off-tissue matrix signal against position is fitted per slide, and converted into the
   log2 ratio that this trend alone would produce between the two section centres. In phase 09,
   a result with the same sign and less than `flag_ratio` = 2 times the median predicted ratio
   is flagged `could_be_position`. This is a warning, not a correction.
2. **Ratios** (Section 4.5) cancel a trend that scales all signals in a pixel.

### 5.2 Sex and run

Sex and run vary together (Section 4.6). No model can separate them; sex results are
descriptive and labelled exploratory.

### 5.3 Slide consistency check (possible label swaps)

For each slide, the log2(Mutant/Control) of every tested peak forms a "difference profile".
If labels are correct, the slides share the genotype effect, so their profiles correlate
positively (or are at least not negatively correlated). Phase 09 computes the Spearman
correlation of the profiles between every pair of slides, in whole section, white matter and
cortex (`pair_check`). A slide is flagged if its median correlation with the other slides is
below −0.3 in at least half of these regions. A purely positional artefact would push every
slide in the same direction, so a sign reversal confined to one slide is what a swap of its two
sections would look like.

The pipeline **never relabels a slide** on the basis of this check: choosing labels from the
results would be circular. It reports the flag and, as exploratory extra columns, the class
tests with the flagged slide left out and flipped. Only independent records (sectioning or
mounting notes, the slide's loading orientation, the region names in the SCiLS project) can
justify setting `qc.sections.flip_slides`.

## 6. What the analysis cannot do

* **Identify isomers.** Annotation is by accurate mass (MS1) only. BMP and PG, galactosyl- and
  glucosylceramide, odd-chain and ether lipids, and some adducts of one class and protonated
  ions of another cannot be distinguished. On-tissue MS/MS or LC-MS/MS of adjacent sections
  would be needed.
* **Quantify absolutely.** MALDI intensities are relative and depend on local ion suppression.
* **Separate genotype from position, or sex from run** (Sections 5.1 and 5.2).
* **Detect small effects with 6 pairs.** Effects need to be consistent in nearly every pair to
  survive correction.

## 7. Correspondence with the CAJAL tutorial

| CAJAL tutorial notebook | Phase | Difference |
|---|---|---|
| 01 Mass spectra and data | 01 to 04 | header and data checks, per-slide QC, recalibration, one peak list across 12 sections |
| 02 Annotation | 05 | class-rule library filtered by LIPID MAPS, ion-chemistry adduct rules, isotope and overlap flags |
| 03 Normalization (uMAIA) | 06 | TIC/RMS/median compared; uMAIA optional in a separate environment |
| 04 Anatomy (Allen) | 07 | automatic plane search and B-spline registration, Dice QC |
| 05 Embedding, 06 Clustering and transfer | 08 | Moran's I, seeded NMF on controls, Leiden territories transferred by kNN |
| 07 Clustering and differential | 09 | per-animal paired tests instead of pixel-level comparisons; class tests; slide consistency check |
| 08 Multimodal (XGBoost, SHAP, GO) | not used | not appropriate with 6 pairs |

## 8. Software

Python 3.11 with the versions pinned in `requirements.txt` (numpy, scipy, pandas, anndata,
pyimzML, lxml, scikit-image, scikit-learn, statsmodels, SimpleITK, brainglobe-atlasapi, igraph,
leidenalg, matplotlib). Reference data: LIPID MAPS Structure Database and HMDB (downloaded and
cached by phase 00), the Allen Mouse Brain CCF at 25 µm (`allen_mouse_25um` via brainglobe).
Random steps (permutations, subsampling, NMF initialization) use `runtime.seed` = 0, so a rerun
of the same code and config gives the same numbers.
