# EPICflow: Meffil-Based Methylation QC and Normalization Pipeline

## Table of Contents

1. [Overview](#overview)
2. [Design Decisions](#design-decisions)
3. [Repository Structure](#repository-structure)
4. [Installation](#installation)
5. [Input Requirements](#input-requirements)
6. [Pipeline Parameters](#pipeline-parameters)
7. [Pipeline Steps](#pipeline-steps)
8. [Running the Pipeline](#running-the-pipeline)
9. [Output Files](#output-files)
10. [Interpreting Results](#interpreting-results)
11. [Cross-Study Identity Checking](#cross-study-identity-checking)
12. [Troubleshooting](#troubleshooting)
13. [Citations](#citations)

---

## Overview

EPICflow processes raw Illumina EPIC/EPICv2 array IDAT files through quality
control, functional normalization, and genotype-based identity checking. It is
designed for the H3Africa ecosystem: large cohorts, African population diversity,
HPC environments (CHPC Lengau/PBS), and integration with H3Africa genotyping
array data processed through
[H3Aflow](https://github.com/Glo-Bimoko/H3Aflow).

The pipeline is written in Nextflow DSL2 and uses meffil as the primary
analysis engine. All R and Python scripts live in `bin/` and are called by the
Nextflow processes.

---

## Design Decisions

This section records every significant architectural decision made during
development so that future maintainers understand why things are done the way
they are.

### Why meffil instead of minfi?

Meffil was chosen over minfi for four reasons specific to this cohort context.
First, meffil performs functional normalization using control probes to model
technical variation, which is substantially better than the quantile
normalization approach in minfi for multi-plate studies. Second, meffil
produces rich per-sample QC objects (stored as `.rds` files) that carry
SNP betas, sex predictions, detection p-values, and control probe metrics
together — making downstream steps like identity checking straightforward
without re-reading IDATs. Third, meffil's HTML QC reports are automatically
generated and are appropriate for sharing with wet-lab collaborators for
plate-level review. Fourth, meffil is faster and more memory-efficient than
minfi at scale.

### Why plate-based processing (12 chips = 96 samples per plate)?

Processing samples in plate units (12 chips, 8 samples per chip, 96 samples
per plate) rather than all at once or chip by chip reflects how Illumina arrays
are physically processed in the laboratory. Chip-to-chip variation within a
plate is a real source of batch effect. By running meffil QC and normalization
at the plate level first, chip effects are visible and correctable within each
batch unit. Cross-plate batch effects are then handled by ComBat in the
combined normalization step. This two-level structure gives the best separation
of within-plate and cross-plate technical variation.

### Why combined normalization after per-plate QC?

Meffil functional normalization requires a reasonable number of samples to
estimate the PCs of technical variation from control probes. Normalizing per
plate first (96 samples) corrects within-plate effects. Then, all plates are
combined and normalized together to correct cross-plate effects. ComBat is
applied automatically if a statistically significant plate effect is detected
(p < 0.05 by ANOVA on the first PC). This matches meffil best practices from
the [Meffil wiki](https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets).

### Why genotype concordance using the 65 rs-probe SNPs?

The EPIC/EPICv2 array contains 65 rs-prefixed probes that target non-CpG SNP
positions. Because these probes measure allelic signal rather than methylation,
bisulfite treatment does not alter their readout. The signal clusters tightly
at beta ≈ 0.0 (AA homozygote), 0.5 (AB heterozygote), and 1.0 (BB homozygote).
These 65 probes act as a genotype fingerprint embedded within every methylation
array, enabling sample identity verification without any additional genotyping.
Within-study all-vs-all concordance catches duplicate samples and sample swaps
between any two participants in the methylation cohort.

### Why allele-call concordance rather than Pearson correlation?

Two methods were evaluated. Pearson correlation on raw betas versus H3A integer
dosages (0/1/2) would in theory retain more information at borderline
heterozygous calls. However, the H3Africa genotyping dataset uses raw PLINK
binary files from H3Aflow (not imputed dosages), so the H3A side only ever
has integer values 0, 1, or 2. The continuous information advantage of Pearson
only materializes when the comparator side also has continuous values. With
integer H3A dosages, Pearson adds computational overhead for no meaningful
sensitivity gain over simple allele-call concordance. Allele-call concordance
was therefore chosen as the sole method.

### Why is allele flip correction not applied?

Allele flip correction is needed for imputed genotype data, where the imputation
software can arbitrarily assign which allele is the "effect allele," producing
negative correlations for genuine matches. The H3A dataset in this pipeline
comes from raw PLINK `.bed/.bim/.fam` files produced by H3Aflow from raw IDATs.
PLINK encodes alleles consistently from the `.bim` ALT column. `match_h3a_snps.r`
handles allele orientation during variant extraction. Flip issues do not arise
with this data provenance, so no flip correction is applied.

### Why positional matching as the primary strategy for H3A SNP lookup?

The H3Africa array `.bim` file uses positional variant IDs in column 2, not
rsIDs. The format is `seq-h3a_37_CHR_POS_REF_ALT` or `h3a_37_CHR_POS_REF_ALT`,
with legacy Affymetrix `kgp` IDs for a small subset. Because rsID-based direct
matching would find almost nothing (the `.bim` does not use `rs` prefixes for
most variants), `match_h3a_snps.r` uses GRCh37/hg19 chromosomal position from
the meffil manifest as its primary matching strategy. A supplementary direct
name match handles the `kgp` entries. This is the correct approach for all
H3Africa array datasets regardless of the specific array version.

### Why pre-QC identity checking?

A sample that was swapped at the laboratory stage may fail QC *because* of the
swap. For example, if the wrong DNA was loaded into the wrong well, the sex
prediction check will flag it, and the bisulfite conversion controls may be
consistent with the wrong sample. If identity checking only runs on QC-passed
samples, swapped samples that failed QC will never be caught. The pipeline
therefore runs identity checking twice: once immediately after `PLATE_QC` on
all samples including failures (`IDENTITY_CHECK_PRE_QC`), and once after
`COMBINED_NORMALIZE` on QC-passed samples only (`H3A_CONCORDANCE`).

### Why is the id_map a CSV with `h3a_id,epic_id` column order?

The H3Africa cohort is larger than the EPIC methylation cohort. Not every
genotyped participant was processed for methylation. The id_map records only
those participants enrolled in both studies. Column order is `h3a_id,epic_id`
because the H3A side is the reference (the larger dataset); H3A IDs that have
no corresponding EPIC entry are the "unenrolled" population against which
`outside_dataset` detection works. Both IDs come from the `Sample_ID` column
used by both H3Aflow and EPICflow samplesheets, so in practice the two columns
often contain identical values.

### Why keep the `.bpm` manifest file off the pipeline?

The `.bpm` (Bead Pool Manifest) is Illumina's proprietary binary format used
by GenomeStudio and Illumina-specific software to decode raw IDAT signals.
Meffil ships with its own internal manifest for EPIC and EPICv2 (accessed via
`meffil.get.features("epicv2")`) that includes probe names, chromosomal
coordinates (GRCh37/hg19), and allele information. The pipeline does not
require the `.bpm` file at any step.

---

## Repository Structure

```
epicflow/
├── main_meffil.nf              # Main Nextflow pipeline
├── nextflow_meffil.config      # Resource profiles (local, PBS)
├── env_meffil_epicv2.yml       # Conda environment specification
├── README.md                   # This file
└── bin/                        # All scripts called by the pipeline
    ├── group_by_plate.r         # Step 1: group IDATs into plate manifests
    ├── create_samplesheet.r     # Step 2: build per-plate samplesheet CSVs
    ├── plate_qc_meffil.r        # Step 3: per-plate meffil QC
    ├── combined_normalize_meffil.r  # Step 4: combined normalization + ComBat
    ├── extract_snp_betas.r      # Step 5: extract SNP betas from QC objects
    ├── pairwise_concordance.py  # Step 5: within-study all-vs-all concordance
    ├── concordance_report.py    # Step 5: HTML report for within-study check
    ├── match_h3a_snps.r         # Step 6/pre: map EPIC probes to H3A .bim
    ├── cross_study_concordance.py  # Step 6/pre: EPIC × H3A identity check
    └── snp-names.txt            # Authoritative list of 65 rs-probe SNP IDs
```

> `snp-names.txt` lists one rsID per line. It is the authoritative probe list
> used to query the H3A `.bim`, independent of which probes survived EPIC QC.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Glo-Bimoko/EPICflow.git
cd EPICflow
```

### 2. Create the conda environment

```bash
conda env create -f env_meffil_epicv2.yml
conda activate meffil_epicv2
```

### 3. Install meffil

Meffil is not available on conda. Install it from GitHub inside R after
activating the environment:

```bash
conda activate meffil_epicv2
Rscript -e "devtools::install_github('perishky/meffil')"
```

Verify the installation:

```bash
Rscript -e "library(meffil); cat(meffil.list.featuresets(), sep='\n')"
```

One of the listed featuresets should be `epicv2`. If you see a different name
(e.g. `IlluminaHumanMethylationEPICv2`), note it, you will need to pass it
via `--h3a_featureset`.

### 4. Install PLINK 1.9 (required for cross-study concordance only)

```bash
# Download from https://www.cog-genomics.org/plink/
# or via conda:
conda install -c bioconda plink
```

If plink is not on `PATH`, pass the full path via `--plink /full/path/to/plink`.

### 5. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
# Move to a directory on PATH, e.g.:
mv nextflow ~/.local/bin/
```

---

## Input Requirements

### IDAT files

Raw IDAT files must be organized into plate-level subdirectories under a single
`idats/` folder. Each subdirectory corresponds to one plate (12 chips, 96
samples). The IDAT filenames must follow Illumina's standard convention:

```
idats/
├── Plate_01/
│   ├── 209042110017_R04C01_Grn.idat
│   ├── 209042110017_R04C01_Red.idat
│   ├── 209042110017_R05C01_Grn.idat
│   ├── 209042110017_R05C01_Red.idat
│   └── ... (96 samples × 2 files = 192 IDATs per plate)
├── Plate_02/
│   └── ...
└── ...
```

The filename format is `<BeadChipBarcode>_<SentrixPosition>_<Color>.idat`.
The BeadChip barcode and Sentrix position together uniquely identify a sample
well. The pipeline extracts the `Sample_ID` (e.g. `G007`) by reading the
samplesheet, see below.

### Samplesheet format

EPICflow uses the same samplesheet format as H3Aflow:

| Sample_ID | BeadChip_Barcode  | Sentrix_Position | Plate_Number | Well_Position | Collected_Gender |
|-----------|-------------------|------------------|--------------|---------------|------------------|
| G007      | 209042110017      | R04C01           | Plate_01     | D01           | Unknown          |
| G009      | 209042110017      | R05C01           | Plate_01     | E01           | Unknown          |
| G012      | 209042110017      | R06C01           | Plate_01     | F01           | Unknown          |

The `Sample_ID` column is the study participant identifier. It must match the
IDs used in the H3Aflow output (`.fam` IID column) and in the id_map CSV.
`Collected_Gender` can be `Unknown` if sex is not recorded; meffil will predict
sex from X/Y chromosome methylation regardless.

### Sample map CSV (`--sample_map`, recommended)

EPICflow can accept a sample map that bridges raw Illumina array identifiers
(BeadChip barcode + Sentrix position) to study-level participant IDs. When
provided, every samplesheet, QC report, beta-value matrix column, and
concordance output will carry participant IDs (e.g. `G001`) instead of raw
barcode strings (e.g. `208789590164_R01C01`), making all results
participant-centric and directly mergeable with clinical metadata.

**Required columns** (column names are matched case-insensitively):

| Column | Description |
|--------|-------------|
| `Sample ID` | Study participant identifier (e.g. `G001`) |
| `BeadChip Barcode` | Illumina chip barcode (numeric, e.g. `208789590164`) |
| `Sentrix Position` | Position on chip (e.g. `R01C01`) |

Optional columns recognised automatically:

| Column | Description |
|--------|-------------|
| `Collected Gender` / `Sex` / `Gender` | Biological sex; mapped to `M` / `F` / `NA` for meffil |
| `Plate Number` | Plate identifier; used to label the `Plate_ID` column |

The file may cover all plates at once — a single CSV for the entire study is
the intended use. The pipeline filters to the rows relevant to each plate
during samplesheet creation.

```
--sample_map ./samplesheets/all_6_epic.csv
```

If `--sample_map` is not provided, the pipeline falls back to using raw
barcode basenames as `Sample_Name` (original behaviour).

### H3A PLINK dataset (for cross-study concordance only)

The H3Africa genotyping dataset must be in PLINK binary format
(`.bed`, `.bim`, `.fam`) with GRCh37/hg19 coordinates in the `.bim` POS
column. This is the standard output of H3Aflow.

The `.bim` file uses positional variant IDs, not rsIDs, in column 2:

```
1   seq-h3a_37_1_61989_G_C    0   61989   C   G
1   h3a_37_1_108310_T_C       0   108310  C   T
1   kgp15717912                0   534247  T   C
```

This is expected and handled by `match_h3a_snps.r` through positional
(chr:pos) lookup against the meffil manifest.

Provide the prefix path (without extension) via `--h3a_bfile`:

```
--h3a_bfile /path/to/h3a_output/qc/h3a_samples
```

This assumes `h3a_samples.bed`, `h3a_samples.bim`, and `h3a_samples.fam`
all exist at that location.

### Identity map CSV (for classified cross-study concordance only)

A CSV file mapping H3Africa Sample_IDs to EPIC Sample_IDs for participants
enrolled in both studies:

```csv
h3a_id,epic_id
G007,G007
G009,G009
G012,G012
G017,G017
```

Both columns use `Sample_ID` values. In practice, because both pipelines use
the same participant ID scheme, the two columns are often identical. The file
must have a header row with exactly these column names (`h3a_id` and `epic_id`,
case-insensitive). H3A participants who were not enrolled for methylation should
not appear in this file — their presence in the H3A dataset but absence from
the map is what enables `outside_dataset` detection.

---

## Pipeline Parameters

All parameters can be set on the command line with `--param_name value` or in
a `nextflow.config` file.

### Core parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--idat_dir` | `./idats` | Directory containing plate subdirectories of IDAT files |
| `--out_dir` | `./results` | Root output directory |
| `--qc_thresh` | `0.95` | Per-sample QC call rate threshold. Samples with fewer than this fraction of probes detected (detection p < 0.01) are excluded from normalization |
| `--sample_map` | `""` | Path to a CSV mapping `BeadChip Barcode` + `Sentrix Position` → `Sample ID`. When provided, all outputs use participant IDs instead of raw barcode strings. See [Sample map CSV](#sample-map-csv---sample_map-recommended) |
| `--conda_env` | `env_meffil_epicv2.yml` | Path to the conda environment YAML |
| `--rscript` | `Rscript` | Path to the Rscript binary (override if using a specific R installation) |

### Within-study concordance parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--conc_threshold` | `0.99` | Concordance at or above this value flags two EPIC samples as likely duplicates or swaps |
| `--conc_min_snps` | `50` | Minimum number of jointly-called SNP positions required for a pairwise comparison to be considered reliable |

### Cross-study concordance parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--h3a_bfile` | `""` | PLINK bfile prefix for the H3Africa genotyping dataset. Leave empty to skip all cross-study steps |
| `--h3a_conc_threshold` | `0.80` | Concordance value used for reporting in the summary. Does **not** control the four-category classification (see thresholds below) |
| `--h3a_conc_min_snps` | `20` | Minimum jointly-called SNPs for a reliable EPIC × H3A comparison |
| `--h3a_featureset` | `"epicv2"` | Meffil featureset name for retrieving EPIC probe GRCh37 coordinates. Run `meffil.list.featuresets()` to verify the correct name on your installation |
| `--plink` | `"plink"` | Path or name of the PLINK 1.9 binary |
| `--snp_names_file` | `bin/snp-names.txt` | Path to the authoritative list of 65 rs-probe SNP IDs. All probes in this file are sought in the H3A `.bim` regardless of whether any were dropped during EPIC QC |
| `--id_map` | `""` | Path to the `h3a_id,epic_id` CSV bridge file. Required for assigned-vs-best classification. Leave empty to report best-match only |

---

## Pipeline Steps

```
IDAT files
    │
    ▼
[Step 1] GROUP_BY_PLATE
    Scans the idat_dir for plate subdirectories.
    Produces one manifest text file per plate listing
    the chip IDs it contains.
    │
    ▼
[Step 2] CREATE_SAMPLESHEETS  (parallel, one per plate)
    Reads each plate manifest and the IDAT filenames
    to build a CSV samplesheet in meffil format.
    │
    ▼
[Step 3] PLATE_QC  (parallel, one per plate)
    Runs meffil.qc() on each plate independently.
    Performs: detection p-value calling, sex prediction,
    bisulfite conversion QC, control probe analysis,
    chip-to-chip batch detection.
    Produces: QC objects (.rds), passed sample lists,
    QC metrics CSV, HTML report per plate.
    │
    ├──────────────────────────────────────────────────────┐
    │                                                      │
    ▼                                                      ▼
[Step 3B] MERGE_SAMPLESHEETS             [Step 5-PRE] IDENTITY_CHECK_PRE_QC
    Combines per-plate CSVs into one            (only if --h3a_bfile set)
    combined_samplesheet.csv for               Uses ALL QC objects including
    joint normalization.                       failures. Extracts SNP betas,
                                               maps to H3A dataset, runs
    ▼                                          EPIC × H3A concordance before
[Step 4] COMBINED_NORMALIZE                    any sample is excluded.
    Runs meffil functional normalization       Catches swaps in samples that
    across all plates together.               failed QC due to the swap itself.
    Applies ComBat automatically if
    cross-plate batch effect detected.
    Produces: BetaValues, MValues,
    PCA plots, normalization report.
    │
    ▼
[Step 5] GENOTYPE_CONCORDANCE
    Extracts 65 rs-probe SNP betas from
    QC-passed samples only. Converts to
    genotype calls (0/1/2). Runs all-vs-all
    pairwise concordance within the EPIC
    cohort. Flags pairs ≥ conc_threshold
    as likely duplicates.
    │
    ▼
[Step 6] H3A_CONCORDANCE
    (only if --h3a_bfile set)
    Maps 65 EPIC probes to H3A .bim by
    chr:pos. Extracts H3A variants with
    PLINK. Runs all EPIC × all H3A
    concordance. Classifies each EPIC
    sample using assigned-vs-best logic.
```

### Step details

**Step 1: GROUP_BY_PLATE**

Scans `--idat_dir` for subdirectories. Each subdirectory is treated as one
plate. Chips within a plate are identified from the IDAT filenames (the barcode
portion, e.g. `209042110017`). A manifest file listing chip IDs is written for
each plate to `results/plate_manifests/`.

**Step 2: CREATE_SAMPLESHEETS**

For each plate manifest, reads the IDAT filenames and constructs a samplesheet
CSV in the format meffil expects: `Sample_Name`, `Basename` (full path to IDAT
prefix), `Sentrix_ID`, `Sentrix_Position`. Runs in parallel for all plates.

**Step 3: PLATE_QC**

Calls `plate_qc_meffil.r` on each plate independently. This script:
- Runs `meffil.qc()` to produce one QC object per sample
- Applies the `--qc_thresh` call rate threshold
- Generates an HTML QC report with control probe heatmaps, sex check
  scatter plot, bisulfite conversion bars, and call rate distributions
- Writes `*_passed_samples.txt` listing samples that passed
- Saves QC objects to `*_qc_objects.rds` for use downstream

Plates run in parallel. A plate with 96 samples typically takes 30–60 minutes.

**Step 3B: MERGE_SAMPLESHEETS**

Concatenates all per-plate samplesheets into a single
`combined_samplesheet.csv` for use by the combined normalization step.

**Step 4: COMBINED_NORMALIZE**

Calls `combined_normalize_meffil.r`, which:
- Loads all QC objects for QC-passed samples
- Runs `meffil.normalize.samples()` across all plates together using
  functional normalization (control probe PCs)
- Tests for significant cross-plate batch effects using ANOVA on PC1
- Applies ComBat batch correction automatically if the plate effect
  p-value is < 0.05
- Writes final beta values, M values, sample info, and PCA plots

**Step 5: GENOTYPE_CONCORDANCE**

Calls `extract_snp_betas.r` on QC-passed samples to extract `$snp.betas`
from each QC object, convert continuous betas to genotype calls
(beta < 0.25 → 0, 0.25–0.75 → 1, > 0.75 → 2), and write a PLINK-style
`.traw` file. Then calls `pairwise_concordance.py` to compute all-vs-all
concordance across EPIC samples. Pairs at or above `--conc_threshold`
(default 0.99) are flagged as likely duplicates.

**Step 5-PRE: IDENTITY_CHECK_PRE_QC**

Identical to the H3A concordance logic (Step 6) but runs on all samples
before QC filtering by passing `--pre_qc` to `extract_snp_betas.r`. This
skips the passed-samples filter so even QC-failed samples contribute their
SNP fingerprint to the check. A swapped sample that failed QC precisely
because it was swapped will be identified here. This step runs in parallel
with Step 3B/4 because it depends only on `PLATE_QC` output.

**Step 6: H3A_CONCORDANCE**

First, `match_h3a_snps.r`:
- Loads the meffil manifest for `--h3a_featureset` to get GRCh37 chr:pos
  for each of the 65 EPIC rs-probes
- Builds a `CHR:POS` index of the H3A `.bim`
- Matches each EPIC probe to the H3A variant at the same position
- Writes an extract list and calls PLINK to produce an H3A `.traw`
- Relabels H3A positional IDs to EPIC rs-probe names in the `.traw`
- Cleans PLINK's `FID_IID` sample column encoding to plain IID

Then, `cross_study_concordance.py`:
- Intersects the SNPs present in both `.traw` files
- Computes concordance for all EPIC × H3A pairs
- For each EPIC sample, records `assigned_score` (vs the mapped H3A ID)
  and `best_score` (vs the best-matching H3A sample across the full cohort)
- Classifies each EPIC sample (see classification section below)

---

## Running the Pipeline

### Minimal run (no cross-study checking)

To be used when one only has EPIC methylation data and want QC, normalization,
and within-study identity checking:

```bash
conda activate meffil_epicv2

nextflow run main_meffil.nf \
    -profile local \
    --idat_dir  ./idats \
    --out_dir   ./results \
    --sample_map ./samplesheets/all_6_epic.csv
```

Omit `--sample_map` to use raw barcode basenames (original behaviour):

```bash
nextflow run main_meffil.nf \
    -profile local \
    --idat_dir ./idats \
    --out_dir  ./results
```

### Full run with cross-study H3Africa identity checking

Use this when you have the H3Aflow-processed H3Africa genotyping dataset and
have prepared the id_map CSV:

```bash
nextflow run main_meffil.nf \
    -profile local \
    --idat_dir   ./idats \
    --out_dir    ./results \
    --sample_map ./samplesheets/all_6_epic.csv \
    --h3a_bfile  ./h3a_output/qc/h3a_samples \
    --id_map     ./samplesheets/id_bridge.csv
```

### Running on CHPC Lengau (PBS)

```bash
nextflow run main_meffil.nf \
    -profile pbs \
    -c nextflow_meffil.config \
    --idat_dir   /scratch/users/yourname/idats \
    --out_dir    /scratch/users/yourname/results \
    --sample_map /scratch/users/yourname/samplesheets/all_6_epic.csv \
    --h3a_bfile  /scratch/users/yourname/h3a_output/qc/h3a_samples \
    --id_map     /scratch/users/yourname/samplesheets/id_bridge.csv
```

### Resuming an interrupted run

Nextflow caches completed processes. If a run fails, add `-resume` to restart
from where it stopped:

```bash
nextflow run main_meffil.nf -profile local -resume \
    --idat_dir ./idats \
    --out_dir  ./results
```

### Adjusting QC stringency

The default call rate threshold of 0.95 is appropriate for most EPIC datasets.
If many samples fail and you want to investigate rather than exclude:

```bash
nextflow run main_meffil.nf -profile local \
    --idat_dir   ./idats \
    --out_dir    ./results \
    --qc_thresh  0.85
```

Open the per-plate HTML reports to understand the failure pattern before
permanently lowering the threshold.

### Specifying a different meffil featureset name

If `meffil.list.featuresets()` shows a different name than `epicv2`:

```bash
nextflow run main_meffil.nf -profile local \
    --idat_dir       ./idats \
    --out_dir        ./results \
    --h3a_bfile      ./h3a_output/qc/h3a_samples \
    --h3a_featureset "epic"
```

---

## Output Files

```
results/
│
├── plate_manifests/
│   └── Plate_01.txt                   # chip IDs per plate
│
├── samplesheets/
│   ├── Plate_01_samplesheet.csv       # per-plate samplesheet
│   └── combined_samplesheet.csv       # all plates merged
│
├── qc_results/                        # ⭐ review these for lab QC
│   ├── Plate_01_QC_report.html        # interactive per-plate QC report
│   ├── Plate_01_qc_metrics.csv        # per-sample QC metrics table
│   ├── Plate_01_passed_samples.txt    # sample IDs that passed QC
│   ├── Plate_01_qc_objects.rds        # meffil QC objects (used downstream)
│   └── ... (one set per plate)
│
├── normalized_combined/               # ⭐ primary outputs for analysis
│   ├── BetaValues_all_plates_combined.csv   # final beta values (use this)
│   ├── MValues_all_plates_combined.csv      # final M values (use this)
│   ├── Sample_info_all.csv                  # sample metadata
│   ├── normalization_report.html            # normalization QC report
│   ├── PCA_before_combat.png                # PCA coloured by plate (before)
│   ├── PCA_after_combat.png                 # PCA coloured by plate (after)
│   └── by_plate/                            # per-plate normalized subsets
│
├── concordance/                       # within-study identity check
│   ├── methylation_snps_concordance_report.html  # ⭐ review for duplicates
│   ├── methylation_snps_flagged_duplicates.tsv   # pairs ≥ conc_threshold
│   ├── methylation_snps_concordance_matrix.csv   # N×N concordance matrix
│   ├── methylation_snps_concordance_pairs.tsv    # all pairwise results
│   ├── methylation_snps_concordance_summary.txt  # plain-text summary
│   ├── methylation_snps_snp_calls.traw           # genotype calls (.traw)
│   └── methylation_snps_snp_betas.csv            # raw SNP betas
│
├── cross_study_pre_qc/                # ⭐ review FIRST if --h3a_bfile set
│   ├── pre_qc_concordance_report.html           # all samples including failures
│   ├── pre_qc_per_sample_verdict.tsv            # per-sample classification
│   ├── pre_qc_classification_summary.tsv        # category counts
│   ├── pre_qc_concordance_matrix.csv
│   ├── pre_qc_concordance_pairs.tsv
│   └── h3a_snps_snp_match_report.txt            # SNP mapping summary
│
└── cross_study_concordance/           # post-QC cross-study check
    ├── cross_study_concordance_report.html      # QC-passed samples only
    ├── cross_study_per_sample_verdict.tsv
    ├── cross_study_classification_summary.tsv
    ├── cross_study_concordance_matrix.csv
    ├── cross_study_concordance_pairs.tsv
    └── h3a_snps_snp_match_report.txt
```

---

## Interpreting Results

### QC report (`qc_results/Plate_*_QC_report.html`)

Open in a browser and share with the wet lab. Key sections:

**Call rate distribution** -> histogram of per-sample call rates. The vertical
line shows `--qc_thresh`. Samples to the left of the line are excluded.
A cluster of failures separated from the main distribution suggests a specific
chip failure; review `bisulfite_by_chip.png` for that plate.

**Sex prediction** -> scatter plot of X vs Y chromosome methylation. Samples
should cluster into two groups (XX upper-left, XY lower-right). Any point in
the wrong cluster, or isolated between clusters, is a likely sample swap or
sex mismatch. Cross-reference with the `Collected_Gender` column.

**Control probe heatmap** -> rows are the 7 control probe categories (bisulfite
conversion, extension, staining, hybridisation, specificity, negative, target
removal). Columns are chips. A column with uniformly poor control scores
indicates a chip-level failure.

**Bisulfite conversion by chip** -> bisulfite conversion efficiency should be
high (> 80%) and consistent across chips. A chip with low conversion failed in
the laboratory; all samples from it should be considered unreliable.

### Normalization report (`normalized_combined/normalization_report.html`)

**PCA before ComBat** -> if samples cluster visibly by plate colour, cross-plate
batch effects are present. ComBat correction is applied automatically.

**PCA after ComBat** -> samples should no longer cluster by plate. If they still
do, the plates may be confounded with biological variation (e.g. all cases on
one plate, all controls on another). Review plate assignment before proceeding.

**ComBat application status** -> the report states whether ComBat was applied
and the p-value that triggered it (< 0.05 threshold).

### Within-study concordance (`concordance/methylation_snps_concordance_report.html`)

The histogram shows the distribution of all pairwise concordance values within
the EPIC cohort. Unrelated individuals typically cluster around 0.25–0.45
(expected sharing by chance with 65 SNPs). Pairs at or above 0.99 are flagged
as likely duplicates.

**Flagged pairs** -> inspect `methylation_snps_flagged_duplicates.tsv`. For each
flagged pair, check:
- Are the `Sample_ID`s the same? If yes, true duplicate, remove one.
- Are they different IDs with concordance = 1.00? Likely a sample swap
  or mislabelling. Both should be excluded or relabelled.
- Concordance between 0.95 and 0.99 may indicate first-degree relatives
  (parent-child or siblings). This is not an error but should be noted.

### Cross-study identity check (`cross_study_pre_qc/`, `cross_study_concordance/`)

The pre-QC report (`cross_study_pre_qc/`) is the more important of the two
because it includes QC-failed samples. Review the `per_sample_verdict.tsv`
first. The post-QC report (`cross_study_concordance/`) confirms that the
remaining passed samples are correctly identified.

---

## Cross-Study Identity Checking

### How it works

`match_h3a_snps.r` retrieves GRCh37 chromosomal coordinates for the 65 EPIC
rs-probe SNPs from the meffil manifest, then scans the H3Africa `.bim` for
variants at the same `CHR:POS`. Because the H3Africa array uses positional
variant IDs (`seq-h3a_37_1_61989_G_C`) rather than rsIDs, coordinate matching
is the primary strategy. PLINK extracts the matched variants and writes a
`.traw` with relabelled SNP IDs (H3A positional IDs replaced by EPIC rs-probe
names). `cross_study_concordance.py` then computes allele-call concordance for
every EPIC sample × every H3A sample pair.

### The id_map file

The id_map CSV (`h3a_id,epic_id`) is what enables the pipeline to ask "does
this EPIC sample match the genotype we expect it to match?" rather than just
"who is this EPIC sample's closest H3A match?". Preparing the id_map is a
one-time step:

```bash
# Example: generate id_map from your participant list
# Assuming both pipelines use the same Sample_ID values:
echo "h3a_id,epic_id" > samplesheets/id_bridge.csv
# Add one line per participant enrolled in both studies:
# G007,G007
# G009,G009
# ... etc.
```

### Classification categories

For each EPIC sample, the pipeline records:
- `assigned_score` -> concordance with the H3A sample specified in the id_map
- `best_score` -> concordance with the best-matching H3A sample in the whole cohort

These two values together drive the classification:

| Category | Condition | Interpretation |
|----------|-----------|----------------|
| `correct_assignment` | best match = assigned H3A ID; score ≥ 0.90 | Sample correctly identified ✓ |
| `swap_high_confidence` | best match ≠ assigned H3A ID; score ≥ 0.90 | Sample likely belongs to a different participant; correct ID can probably be recovered ⚠ |
| `swap_moderate_confidence` | best match ≠ assigned H3A ID; score 0.70–0.90 | Possible swap; investigate before correcting ⚠ |
| `outside_dataset` | best match is an H3A ID **not in the id_map** | EPIC sample's DNA matches a participant who was genotyped but never enrolled for methylation; possible contamination or enrolment error 🔍 |
| `low_confidence` | best score < 0.70 | No reliable match found; DNA possibly absent from H3A cohort, or sample severely degraded ✗ |

Note: the classification thresholds (high ≥ 0.90, moderate ≥ 0.70) are fixed
in `cross_study_concordance.py` and are not configurable by parameter. The
`--h3a_conc_threshold` parameter controls only the reporting cutoff in the
summary statistics, not the classification logic.

### What to do with each category

**`correct_assignment`** - no action required.

**`swap_high_confidence`** - check the `best_h3a` column in
`per_sample_verdict.tsv`. The correct H3A ID is the `best_h3a` value. Check
whether any other EPIC sample claims that H3A ID as its assigned match (a
reciprocal swap). Update your sample tracking accordingly and relabel or
exclude as appropriate.

**`swap_moderate_confidence`** - do not relabel automatically. Pull the raw
concordance values from `concordance_pairs.tsv` and inspect the distribution
for this sample. If the best score is clearly separated from the second-best
(e.g. 0.85 vs 0.40), the swap is likely real. If the top several H3A matches
are all close together (e.g. 0.83, 0.81, 0.79), the signal is ambiguous,
consider excluding the sample.

**`outside_dataset`** - the `best_h3a` column contains the H3A ID of the
participant whose DNA is in the tube. Check your enrolment records. This can
indicate: (a) a tube labelling error where the wrong participant's sample was
collected, (b) cross-contamination during DNA extraction, or (c) a participant
who was genotyped but whose methylation sample was excluded from the EPIC
run for an unrelated reason (in which case this is a false positive,
check the id_map to confirm they are absent).

**`low_confidence`** - most commonly this means the participant was not
present in the H3A dataset at all (recruited for methylation only, or the
H3A data failed QC for that person). It can also indicate severe DNA
degradation. Cross-reference with the within-study concordance report, if
the same sample has low concordance with all other EPIC samples too, the
DNA itself is likely degraded.

### SNP matching report (`h3a_snps_snp_match_report.txt`)

Review this file after every run that includes `--h3a_bfile`. It reports:
- Total EPIC probes sought (up to 65)
- Number matched by position
- Number matched by direct name (rare for H3Africa data)
- Number unmatched (probes where the H3Africa array has no variant at
  that chromosomal position)

If fewer than ~40 probes match, concordance results will be unreliable.
The most common cause is a genome build mismatch: meffil uses GRCh37/hg19
coordinates, and the H3A `.bim` must also use GRCh37/hg19 in the POS column.
If H3Aflow lifted over to GRCh38, coordinate matching will fail. Check by
comparing a known probe position: `rs10033147` is at chr1:161476960 on hg19.

---

## Troubleshooting

### Meffil cannot find EPICv2 featureset

```
Error: Could not load meffil SNP annotations for any featureset.
```

Run `meffil.list.featuresets()` in R to see what is installed. Pass the
correct name:

```bash
--h3a_featureset "epic"
```

### No SNPs match the H3A .bim

```
No EPIC probes could be matched to the H3A dataset.
```

Almost always a genome build mismatch. Verify:

```r
# Check meffil coordinates for a known probe
snp <- meffil.get.features("epicv2")
snp[snp$name == "rs10033147", c("name","chromosome","position")]
# Should show: chr1, 161476960 (hg19)
```

Then check the `.bim`:
```bash
grep "161476960" /path/to/h3a_samples.bim
# Should find a match if the .bim is hg19
```

### Many samples fail QC

If more than 10% of samples fail the `--qc_thresh` filter, first check the
QC report before changing the threshold. Look for:

- All failures on one chip → that chip may have failed during hybridisation
- Failures clustered by plate → plate-level laboratory problem
- Failures with low bisulfite conversion → insufficient bisulfite treatment

Only lower `--qc_thresh` after understanding the cause.

### PLINK exits with a non-zero code

```
PLINK exited with code 127.
```

PLINK is not on `PATH`. Either install it (`conda install -c bioconda plink`)
or pass the full path:

```bash
--plink /full/path/to/plink
```

### Memory errors during COMBINED_NORMALIZE

The default PBS allocation is 64 GB. For cohorts > 500 samples, increase it
in `nextflow_meffil.config`:

```groovy
withName: 'COMBINED_NORMALIZE' {
    memory = '128 GB'
    time   = '12h'
}
```

### Plate effects remain after ComBat

If the post-normalization PCA still shows clustering by plate, the plates may
be confounded with biological variables. Check whether cases and controls are
evenly distributed across plates. If confounding exists, ComBat must be run
with a biological covariate:

```r
# In combined_normalize_meffil.r, find the ComBat call and modify:
mod <- model.matrix(~Disease_Status, data = pheno_data)
combat_beta <- ComBat(dat = beta_matrix, batch = plate_vector, mod = mod)
```

---

## Citations

If you use this pipeline, please cite:

**Meffil**
Min, J.L., Hemani, G., Davey Smith, G., Relton, C., & Suderman, M. (2018).
Meffil: efficient normalization and analysis of very large DNA methylation
datasets. *Bioinformatics*, 34(23), 3983–3989.
https://doi.org/10.1093/bioinformatics/bty594

**Minfi** (underlying dependency)
Aryee, M.J., et al. (2014). Minfi: a flexible and comprehensive Bioconductor
package for the analysis of Infinium DNA methylation microarrays.
*Bioinformatics*, 30(10), 1363–1369.
https://doi.org/10.1093/bioinformatics/btu049

**ComBat** (batch correction)
Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in
microarray expression data using empirical Bayes methods.
*Biostatistics*, 8(1), 118–127.
https://doi.org/10.1093/biostatistics/kxj037


**Nextflow**
Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational
workflows. *Nature Biotechnology*, 35, 316–319.
https://doi.org/10.1038/nbt.3820

**H3Aflow** (companion pipeline for H3Africa SNP array processing)
https://github.com/Glo-Bimoko/H3Aflow

---

## Future Direction: MetaFlow

EPICflow and H3Aflow are designed to eventually be unified into **MetaFlow**,
a single end-to-end pipeline that processes H3Africa SNP array IDATs and
EPIC methylation array IDATs together from raw files through GWAS, EWAS, and
multi-omics integration. MetaFlow will subsume the cross-study concordance
steps currently handled by both pipelines separately, and support joint analysis of the full H3Africa
consortium dataset across 30+ countries.
