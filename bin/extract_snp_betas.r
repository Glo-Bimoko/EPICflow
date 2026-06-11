#!/usr/bin/env Rscript
# extract_snp_betas.r
# Extract rs-probe (SNP) betas from meffil QC objects and convert to genotype
# calls (0/1/2) for pairwise concordance analysis.
#
# SNP probes on the EPIC/EPICv2 array are NOT subject to bisulfite conversion
# in the same way as CpG probes. They measure allele-specific signal directly
# and cluster tightly at ~0.0 (AA), ~0.5 (AB), ~1.0 (BB), making them ideal
# for identity checking.
#
# Output: a PLINK-style .traw file (transposed genotype matrix) that is
# consumed directly by pairwise_concordance.py
#
# Usage:
#   Rscript extract_snp_betas.r <qc_objects_dir> <passed_samples_file> <out_prefix>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste0(
    "Usage: Rscript extract_snp_betas.r <qc_objects_dir> <passed_samples_dir> <out_prefix> [--pre_qc]\n\n",
    "  <qc_objects_dir>    : directory of *_qc_objects.rds files\n",
    "  <passed_samples_dir>: directory of *_passed_samples.txt files\n",
    "  <out_prefix>        : output file prefix\n",
    "  --pre_qc            : export ALL samples regardless of QC pass/fail status\n",
    "                        (for pre-QC identity checking before sample exclusion).\n",
    "                        When set, the passed_samples filter is skipped.\n"
  ))
}

qc_objects_dir     <- args[1]
passed_samples_dir <- args[2]
out_prefix         <- args[3]

# Optional --pre_qc flag: use all samples regardless of QC outcome
pre_qc_mode <- "--pre_qc" %in% args

cat("\n", rep("=", 60), "\n", sep = "")
cat("SNP BETA EXTRACTION FOR GENOTYPE CONCORDANCE\n")
cat(rep("=", 60), "\n\n")
if (pre_qc_mode) {
  cat("MODE: PRE-QC  — all samples included regardless of QC outcome\n")
  cat("      (identity checking before sample exclusion)\n\n")
} else {
  cat("MODE: POST-QC — only QC-passed samples included\n\n")
}

suppressPackageStartupMessages({
  library(meffil)
})

# ============================================================
# STEP 1: LOAD QC OBJECTS
# ============================================================
cat("--- STEP 1: Loading QC objects ---\n")

qc_files <- list.files(qc_objects_dir,
                       pattern = "_qc_objects\\.rds$",
                       full.names = TRUE)

if (length(qc_files) == 0) {
  stop("No QC object files found in: ", qc_objects_dir)
}

cat("Found", length(qc_files), "plate QC object file(s)\n")

all_qc <- list()
for (f in qc_files) {
  plate_id <- sub("_qc_objects\\.rds$", "", basename(f))
  cat("  Loading", plate_id, "...\n")
  plate_qc <- readRDS(f)
  for (s in names(plate_qc)) {
    all_qc[[s]] <- plate_qc[[s]]
  }
}
cat("Total QC objects loaded:", length(all_qc), "\n")

# ============================================================
# STEP 2: FILTER TO PASSED SAMPLES (skipped in pre_qc_mode)
# ============================================================
if (pre_qc_mode) {
  cat("\n--- STEP 2: Pre-QC mode — retaining all", length(all_qc), "samples ---\n")
  cat("  (passed_samples filter skipped; identity check runs on ALL samples)\n")
} else {
  cat("\n--- STEP 2: Filtering to QC-passed samples ---\n")

  passed_files <- list.files(passed_samples_dir,
                             pattern = "_passed_samples\\.txt$",
                             full.names = TRUE)

  passed_samples <- character(0)
  for (f in passed_files) {
    samps <- readLines(f)
    samps <- samps[nchar(trimws(samps)) > 0]
    passed_samples <- c(passed_samples, samps)
  }
  passed_samples <- unique(passed_samples)
  cat("Passed samples across all plates:", length(passed_samples), "\n")

  available <- intersect(passed_samples, names(all_qc))
  cat("QC objects available for passed samples:", length(available), "\n")

  if (length(available) == 0) {
    stop("No passed samples have QC objects available.")
  }

  all_qc <- all_qc[available]
}

# ============================================================
# STEP 3: EXTRACT SNP BETAS
# ============================================================
cat("\n--- STEP 3: Extracting SNP betas (rs-probe intensities) ---\n")
cat("Note: rs probes measure allelic signal AFTER bisulfite conversion.\n")
cat("      Because these probes target non-CpG SNP positions, bisulfite\n")
cat("      treatment does not alter the underlying base (the SNP allele);\n")
cat("      it only converts unmethylated cytosines nearby. Signal clusters\n")
cat("      remain at 0.0/0.5/1.0 irrespective of bisulfite treatment.\n\n")

# meffil stores SNP betas inside each QC object as $snp.betas
# These are the same rs-prefixed probes used internally for concordance QC
snp_betas_list <- lapply(all_qc, function(qc_obj) {
  if (!is.null(qc_obj$snp.betas)) {
    return(qc_obj$snp.betas)
  }
  return(NULL)
})

# Drop NULLs
snp_betas_list <- Filter(Negate(is.null), snp_betas_list)
cat("Samples with SNP betas:", length(snp_betas_list), "\n")

if (length(snp_betas_list) == 0) {
  stop("No SNP betas found in QC objects. ",
       "Ensure meffil.qc() was run with chip that includes rs probes.")
}

# Build matrix: probes x samples
all_probes <- Reduce(union, lapply(snp_betas_list, names))
rs_probes  <- all_probes[grepl("^rs", all_probes)]
cat("rs-prefixed SNP probes found:", length(rs_probes), "\n")

if (length(rs_probes) == 0) {
  # Fallback: use whatever probes are present (all are SNP probes here)
  rs_probes <- all_probes
  cat("Using all", length(rs_probes), "SNP probes (no rs prefix found)\n")
}

snp_matrix <- matrix(NA_real_,
                     nrow = length(rs_probes),
                     ncol = length(snp_betas_list),
                     dimnames = list(rs_probes, names(snp_betas_list)))

for (samp in names(snp_betas_list)) {
  betas <- snp_betas_list[[samp]]
  shared <- intersect(rs_probes, names(betas))
  snp_matrix[shared, samp] <- betas[shared]
}

cat("SNP beta matrix dimensions:", nrow(snp_matrix), "probes x",
    ncol(snp_matrix), "samples\n")

# Save raw betas for reference
beta_out <- paste0(out_prefix, "_snp_betas.csv")
write.csv(snp_matrix, beta_out, row.names = TRUE)
cat("✓ Raw SNP betas saved:", beta_out, "\n")

# ============================================================
# STEP 4: CONVERT BETAS TO GENOTYPE CALLS (0/1/2)
# ============================================================
cat("\n--- STEP 4: Converting betas to genotype calls ---\n")
cat("Thresholds: AA = beta < 0.25 → 0\n")
cat("            AB = 0.25 ≤ beta ≤ 0.75 → 1\n")
cat("            BB = beta > 0.75 → 2\n")
cat("            Missing / ambiguous → NA\n\n")

beta_to_genotype <- function(beta_vec) {
  gt <- rep(NA_integer_, length(beta_vec))
  gt[!is.na(beta_vec) & beta_vec <  0.25] <- 0L
  gt[!is.na(beta_vec) & beta_vec >= 0.25 & beta_vec <= 0.75] <- 1L
  gt[!is.na(beta_vec) & beta_vec >  0.75] <- 2L
  gt
}

gt_matrix <- apply(snp_matrix, 2, beta_to_genotype)
rownames(gt_matrix) <- rownames(snp_matrix)

# Report call rate per sample
call_rates <- apply(gt_matrix, 2, function(x) sum(!is.na(x)) / length(x))
cat("SNP call rate summary across samples:\n")
print(summary(call_rates))

low_call <- names(call_rates)[call_rates < 0.80]
if (length(low_call) > 0) {
  cat("\n⚠ Samples with SNP call rate < 80%:\n")
  for (s in low_call) cat("  -", s, sprintf("(%.1f%%)\n", call_rates[s]*100))
}

# Genotype frequency summary (to flag degenerate probes)
probe_missing <- apply(gt_matrix, 1, function(x) mean(is.na(x)))
bad_probes <- names(probe_missing)[probe_missing > 0.20]
if (length(bad_probes) > 0) {
  cat("\n⚠", length(bad_probes), "probes with >20% missing calls — will be excluded:\n")
  cat("  ", paste(bad_probes, collapse = ", "), "\n")
  gt_matrix <- gt_matrix[!rownames(gt_matrix) %in% bad_probes, , drop = FALSE]
}

cat("\nFinal genotype matrix:", nrow(gt_matrix), "probes x",
    ncol(gt_matrix), "samples\n")

# ============================================================
# STEP 5: WRITE PLINK-STYLE .traw FILE
# ============================================================
# Format: CHR  SNP  CM  POS  COUNTED  ALT  <sample1>  <sample2> ...
# CM and POS set to 0/1 as placeholders (not used downstream)
# COUNTED allele = B (coded as 2), ALT = A (coded as 0)
cat("\n--- STEP 5: Writing .traw genotype file ---\n")

traw_out <- paste0(out_prefix, "_snp_calls.traw")

# Build header
sample_ids <- colnames(gt_matrix)
header <- paste(c("CHR", "SNP", "CM", "POS", "COUNTED", "ALT", sample_ids),
                collapse = "\t")

# Build rows
traw_rows <- apply(gt_matrix, 1, function(row) {
  snp_id <- names(row)[1]  # will be set per probe below
  calls <- ifelse(is.na(row), "NA", as.character(row))
  paste(calls, collapse = "\t")
})

lines <- c(header)
for (probe in rownames(gt_matrix)) {
  calls <- gt_matrix[probe, ]
  call_str <- ifelse(is.na(calls), "NA", as.character(calls))
  line <- paste(c("0", probe, "0", "0", "B", "A", call_str), collapse = "\t")
  lines <- c(lines, line)
}

writeLines(lines, traw_out)
cat("✓ .traw genotype file written:", traw_out, "\n")
cat("  Rows (SNPs):", nrow(gt_matrix), "\n")
cat("  Columns (samples):", ncol(gt_matrix), "\n")

# ============================================================
# STEP 6: WRITE SAMPLE ID FILE (for python script)
# ============================================================
sample_id_out <- paste0(out_prefix, "_sample_ids.txt")
writeLines(sample_ids, sample_id_out)
cat("✓ Sample ID list written:", sample_id_out, "\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("SNP EXTRACTION COMPLETED\n")
cat(rep("=", 60), "\n")
cat("  Samples processed  :", ncol(gt_matrix), "\n")
cat("  SNP probes used    :", nrow(gt_matrix), "\n")
cat("  Output .traw       :", traw_out, "\n")
cat("  Raw betas CSV      :", beta_out, "\n")
cat(rep("=", 60), "\n")
