# Meffil-Based Methylation QC Pipeline (Plate-Based)

## Overview

This pipeline processes Illumina EPIC array methylation data using **meffil** with a **plate-based approach** (12 chips = 96 samples per plate). This design is optimal for detecting and correcting batch effects.

## Key Features

✅ **Plate-based QC**: Processes 12 chips (96 samples) together  
✅ **Comprehensive QC**: Sex prediction, control probes, bisulfite conversion  
✅ **Within-plate batch correction**: Meffil functional normalization  
✅ **Cross-plate batch correction**: Automatic ComBat if needed  
✅ **Lab process validation**: Detailed QC reports for lab review  
✅ **Automatic batch detection**: Statistical testing at each step  

## Files Needed

- `env_meffil_epicv2.yml` - Conda environment
- `main_meffil.nf` - Main Nextflow pipeline
- `nextflow_meffil.config` - Configuration
- `group_by_plate.r` - Groups chips into plates
- `plate_qc_meffil.r` - QC per plate
- `plate_normalize_meffil.r` - Normalize per plate
- `merge_plates.r` - Merge all plates



## Directory Structure

```
project/
├── main_meffil.nf
├── nextflow_meffil.config
├── env_meffil_epicv2.yml
├── group_by_plate.r
├── plate_qc_meffil.r
├── plate_normalize_meffil.r
├── merge_plates.r
├── idats/
│   ├── 200000001234_R01C01_Grn.idat
│   ├── 200000001234_R01C01_Red.idat
│   └── ... (all IDAT files)
└── results/ (created by pipeline)
```

## Installation

### 1. Create Conda Environment

```bash
conda env create -f env_meffil_epicv2.yml
conda activate meffil_epicv2
```

### 2. Install Meffil

Meffil is not in conda, so install via R:

```r
# In R console
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("perishky/meffil")
```

### 3. Verify Installation

```r
library(meffil)
library(minfi)
```

## Usage

### Basic Usage (Local)

```bash
nextflow run main_meffil.nf \
    -profile local \
    --idat_dir ./idats \
    --out_dir ./results \
    --qc_thresh 0.95 \
    --chips_per_plate 12
```

### PBS Cluster Usage

```bash
nextflow run main_meffil.nf \
    -profile pbs \
    -c nextflow_meffil.config \
    --idat_dir /path/to/idats \
    --out_dir /path/to/results \
    --qc_thresh 0.95 \
    --chips_per_plate 12
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--idat_dir` | `./idats` | Directory containing IDAT files |
| `--out_dir` | `./results` | Output directory |
| `--qc_thresh` | `0.95` | QC threshold (call rate) |
| `--chips_per_plate` | `12` | Number of chips per plate (8 samples/chip = 96 samples/plate) |
| `--conc_threshold` | `0.99` | Within-study duplicate flag threshold |
| `--conc_min_snps` | `50` | Min jointly-called SNPs for a reliable within-study comparison |
| `--h3a_bfile` | `""` | PLINK bfile prefix for H3Africa dataset (skips cross-study step if empty) |
| `--h3a_conc_threshold` | `0.80` | Cross-study identity match threshold (lower than within-study) |
| `--h3a_conc_min_snps` | `20` | Min SNPs for a reliable cross-study comparison |
| `--snp_names_file` | `bin/snp-names.txt` | Authoritative list of rs-probe SNP IDs (one rsID per line). When set, all listed SNPs are sought in the H3A `.bim` regardless of EPIC QC dropout. |
| `--h3a_featureset` | `epicv2` | meffil featureset name for positional SNP fallback |
| `--plink` | `plink` | Path to PLINK 1.9 binary |

---

## Cross-Study Concordance Against H3Africa Genotyping Data

The EPIC/methylation participants are a subset of a larger H3Africa cohort that
was also genotyped with an H3Africa SNP array.  This step confirms sample
identity across the two platforms — catching swaps, contamination, or
mis-labelling between the methylation and genotyping arms of your study.

### Recommended workflow

Run H3Aflow on the H3Africa IDATs first to produce a quality-controlled PLINK
dataset, then point this pipeline at the output bfile:

```bash
# 1. Process H3Africa IDATs through H3Aflow
#    → produces h3a_output/qc/h3a_samples.bed/.bim/.fam (or similar)
nextflow run Glo-Bimoko/H3Aflow \
    --input_dir ./h3a_idats \
    --outdir    ./h3a_output \
    -profile    local

# 2. Run this pipeline with the H3Aflow bfile as input
nextflow run main_meffil.nf \
    -profile local \
    --idat_dir      ./epic_idats \
    --out_dir       ./results \
    --h3a_bfile     ./h3a_output/qc/h3a_samples \
    --snp_names_file bin/snp-names.txt
```

### How the SNP filtering works

`snp-names.txt` (in `bin/`) holds the canonical 65 rsID probes present on
the EPIC/EPICv2 array.  `match_h3a_snps.r` uses this list as its extraction
target, meaning:

- **All 65 SNPs are sought in the H3A `.bim`**, even if some were flagged
  and dropped during EPIC methylation QC (low call rate, etc.).
- Probes matched by direct rsID in the `.bim` are extracted as-is.
- Probes not found by rsID fall back to chromosomal position matching via
  the meffil array manifest.
- Any probes still unmatched are reported but do not cause the step to fail.
- PLINK extracts only the matched variants; concordance is then computed on
  the intersection present in both `.traw` files.

Because H3Africa cohorts typically have more participants than the EPIC subset,
the cross-study concordance output reports **all EPIC × all H3A pairs**, not
just 1:1.  This lets you identify:
- ✅ Confirmed matches (EPIC sample found in H3A at expected concordance)
- ⚠ Missing matches (EPIC sample has no H3A counterpart → possible swap or
  enrolment mismatch)
- 🚨 Unexpected matches (EPIC sample matches a different H3A participant →
  contamination or labelling error)

## Pipeline Steps

### Step 1: Group by Plate
- Groups chips into plates (12 chips = 1 plate)
- Creates manifest files for each plate
- **Output**: `results/plate_manifests/plate_*.txt`

### Step 2: Plate QC
- Runs meffil QC on each plate independently
- Detects within-plate batch effects (chip-to-chip)
- Comprehensive control probe analysis
- Sex prediction for sample swap detection
- **Output**: 
  - `results/qc_results/plate_*_QC_detailed.csv`
  - `results/qc_results/plate_*_QC_report.html`
  - `results/qc_results/plate_*_passed_samples.txt`
  - QC plots (call rates, sex check, controls, etc.)

### Step 3: Plate Normalization
- Meffil functional normalization per plate
- Corrects within-plate batch effects
- Uses control probes to model technical variation
- **Output**:
  - `results/normalized_plates/plate_*_beta.csv`
  - `results/normalized_plates/plate_*_M.csv`
  - `results/normalized_plates/plate_*_norm_report.html`
  - `results/normalized_plates/plate_*_PCA_normalized.png`

### Step 4: Merge Plates
- Combines all normalized plates
- Tests for cross-plate batch effects
- Applies ComBat if needed (p < 0.05)
- **Output**:
  - `results/final/BetaValues_FINAL.csv` ⭐ **USE THIS**
  - `results/final/MValues_FINAL.csv` ⭐ **USE THIS**
  - `results/final/Sample_info_all.csv`
  - `results/final/BetaValues_merged_plates.csv` (intermediate)
  - `results/final/PCA_all_plates_*.png`

## Output Files

### Final Files (for downstream analysis)
```
results/final/
├── BetaValues_FINAL.csv          ⭐ Main output (beta values)
├── MValues_FINAL.csv              ⭐ Main output (M values)
├── Sample_info_all.csv            📋 Sample metadata
└── PCA_all_plates_*.png           📊 QC plots
```

### QC Reports (for lab review)
```
results/qc_results/
├── plate_1_QC_report.html         📄 Interactive QC report
├── plate_1_QC_detailed.csv        📊 Detailed QC metrics
├── plate_1_callrate.png           📈 Call rate plots
├── plate_1_sex_check.png          🔍 Sample swap detection
├── plate_1_control_heatmap.png    🔥 Control probe performance
└── plate_1_bisulfite_by_chip.png  🧪 Lab process validation
```

## Understanding the Results

### QC Metrics

**Detection P-value**: Compares probe signal to negative controls
- **p < 0.01**: Probe detected (good signal)
- **Call rate**: % probes with p < 0.01
- **Threshold**: Default 0.95 (95% probes detected)

**Sex Prediction**: 
- Checks X/Y chromosome methylation
- Flags sample swaps if predicted ≠ expected
- XY difference plot shows separation

**Control Probes**:
- **Bisulfite conversion**: Should be high (>1000)
- **Extension**: Red/Green should be similar
- **Staining**: Checks labeling worked
- Low values = lab process failure

### Batch Effects

**Within-plate (chip-to-chip)**:
- Corrected by meffil functional normalization
- Check: `plate_*_PCA_normalized.png`
- Should NOT cluster by chip after normalization

**Cross-plate**:
- Tested automatically (p-value reported)
- ComBat applied if p < 0.05
- Check: `PCA_all_plates_after_combat.png`
- Should NOT cluster by plate after correction

### When to Worry

🚨 **High failure rate** (>10% samples fail):
- Check bisulfite conversion controls
- Review lab protocols
- Consider re-running failed samples

🚨 **Sex mismatches**:
- Sample swap occurred
- Check sample tracking
- Exclude or re-label affected samples

🚨 **Significant chip effects persist after normalization**:
- One chip may have failed
- Consider excluding entire chip
- Review processing for that chip

🚨 **Strong plate effects after ComBat**:
- Plates may be confounded with biology
- Check plate assignments
- May need different normalization strategy

## Advantages Over Per-Chip Approach

| Feature | Old (Per-Chip) | New (Plate-Based) |
|---------|----------------|-------------------|
| **Batch unit** | 8 samples | 96 samples |
| **Within-plate detection** | ❌ No | ✅ Yes |
| **Chip-to-chip variation** | Not visible | Visible & corrected |
| **Lab QC validation** | Basic | Comprehensive |
| **Sample swap detection** | ❌ No | ✅ Yes |
| **Control probe QC** | Limited | Detailed |
| **Normalization** | All samples together | Per plate, then merge |
| **Batch correction** | Manual | Automatic |

## Meffil vs Minfi Summary

| Feature | Minfi | Meffil |
|---------|-------|--------|
| **QC comprehensiveness** | Basic | ⭐ Extensive |
| **Sample swap detection** | ❌ | ✅ Sex prediction |
| **Control probe analysis** | Limited | ⭐ Detailed (7 types) |
| **Lab process validation** | Basic | ⭐ Excellent |
| **HTML reports** | Manual | ✅ Automatic |
| **Batch handling** | Manual ComBat | ⭐ Built-in functional norm |
| **Speed** | Moderate | ⭐ Faster |
| **Memory efficiency** | Good | ⭐ Better |

## Troubleshooting

### Problem: Meffil installation fails

**Solution**:
```r
# Try installing dependencies first
BiocManager::install(c("minfi", "preprocessCore", "GenomicRanges"))
devtools::install_github("perishky/meffil")
```

### Problem: "chip not recognized" error

**Solution**: Meffil may not recognize EPICv2. In scripts, try:
```r
qc_objects <- meffil.qc(..., chip = "epic")  # Try epic first
# If fails, try: chip = "epicv2"
```

### Problem: Memory errors during normalization

**Solution**: Increase memory allocation in `nextflow_meffil.config`:
```groovy
withName: 'PLATE_NORMALIZE' {
    memory = '128 GB'  // Increase from 64 GB
}
```

### Problem: No samples pass QC

**Solution**: Lower threshold temporarily to diagnose:
```bash
--qc_thresh 0.85  # Lower from 0.95
```
Check QC reports to identify specific failures.

### Problem: Plate effects remain after ComBat

**Solution**:
1. Check if plates confounded with biology
2. Include biological covariates in ComBat model
3. Edit `merge_plates.r`, line with `model.matrix`:
   ```r
   # Add phenotype data
   mod <- model.matrix(~Disease_Status, data = pheno_data)
   ```

## Expected Runtime

| Step | 96 samples (1 plate) | 192 samples (2 plates) |
|------|---------------------|------------------------|
| Group by plate | <1 min | <1 min |
| Plate QC | 30-60 min | 30-60 min each |
| Plate normalization | 1-2 hours | 1-2 hours each |
| Merge plates | 30 min | 1 hour |
| **Total** | **~2-4 hours** | **~4-8 hours** |

*Times assume 4-8 cores per process*

## Citation

If using this pipeline, please cite:

**Meffil**:
- Min, J.L., Hemani, G., Davey Smith, G., Relton, C., & Suderman, M. (2018). Meffil: efficient normalization and analysis of very large DNA methylation datasets. *Bioinformatics*, 34(23), 3983-3989.

**Minfi** (underlying dependency):
- Aryee, M.J., et al. (2014). Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. *Bioinformatics*, 30(10), 1363-1369.

**ComBat**:
- Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118-127.

## Support

For issues specific to:
- **Meffil**: https://github.com/perishky/meffil/issues
- **Nextflow**: https://www.nextflow.io/docs/latest/
- **This pipeline**: Check logs in `results/reports/`

## Quick Start Example

```bash
# 1. Activate environment
conda activate meffil_epicv2

# 2. Test with small dataset
nextflow run main_meffil.nf \
    -profile local \
    --idat_dir ./test_idats \
    --out_dir ./test_results \
    --chips_per_plate 12

# 3. Check results
ls test_results/final/BetaValues_FINAL.csv

# 4. Review QC
firefox test_results/qc_results/plate_1_QC_report.html
```

## Notes

- **Empty wells**: Not used as controls; built-in negative controls are better
- **Reference samples**: Consider adding same sample to each plate for cross-plate calibration
- **Chip format**: Assumes format `ChipID_Position` (e.g., `200000001234_R01C01`)
- **Parallel processing**: Pipeline parallelizes across plates automatically

## Key strength of Pipeline

1. ✅ **Processes 12 chips together** (not individually)
2. ✅ **Uses meffil** (not just minfi)
3. ✅ **Automatic batch detection** (statistical testing)
4. ✅ **Comprehensive QC** (sex, controls, swaps)
5. ✅ **HTML reports** (interactive, shareable)
6. ✅ **Two-level batch correction** (within-plate + cross-plate)
7. ✅ **Better for lab QC validation**