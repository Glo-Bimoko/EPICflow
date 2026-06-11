#!/usr/bin/env Rscript
# plate_qc_meffil.R
# QC for entire plate following meffil best practices
# Based on: https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plate_qc_meffil.R <samplesheet.csv> <qc_out_dir> <qc_threshold>")
}

samplesheet_file <- args[1]
qc_out_dir <- args[2]
qc_threshold <- as.numeric(args[3])

dir.create(qc_out_dir, showWarnings = FALSE, recursive = TRUE)

# Extract plate ID from samplesheet filename
plate_id <- sub("_samplesheet\\.csv$", "", basename(samplesheet_file))

cat("\n", rep("=", 60), "\n", sep = "")
cat("PLATE QC (Meffil Best Practices):", plate_id, "\n")
cat(rep("=", 60), "\n\n")

# Setup threading
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 
                      Sys.getenv("PBS_NUM_PPN", 4)))
cat("Using", n_cores, "cores\n")
options(mc.cores = n_cores)

# Install meffil if needed
suppressPackageStartupMessages({
  if (!requireNamespace("meffil", quietly = TRUE)) {
    cat("Installing meffil from GitHub...\n")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    devtools::install_github("perishky/meffil", quiet = TRUE)
  }
  library(meffil)
  library(ggplot2)
})

# ============================================================
# STEP 0: READ SAMPLESHEET AND SETUP
# ============================================================
cat("--- STEP 0: Reading samplesheet ---\n")

# Read samplesheet
samplesheet <- read.csv(samplesheet_file, stringsAsFactors = FALSE)
cat("Found", nrow(samplesheet), "samples in samplesheet\n")

# Check required columns
required_cols <- c("Sample_Name", "Basename")
missing_cols <- setdiff(required_cols, colnames(samplesheet))
if (length(missing_cols) > 0) {
  stop("Samplesheet missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# STEP 1: CREATE QC OBJECTS (Meffil recommended approach)
# ============================================================
cat("\n--- STEP 1: Creating QC objects ---\n")

# Set QC parameters as recommended by meffil
qc_params <- meffil.qc.parameters(
  detectionp.samples.threshold = 1 - qc_threshold,
  detectionp.cpgs.threshold = 0.05,
  beadnum.samples.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,  # More permissive as recommended
  snp.concordance.threshold = 0.9,  # More permissive
  sample.genotype.concordance.threshold = 0.8
)

# Create QC objects
cat("Creating QC objects for", nrow(samplesheet), "samples...\n")
qc_objects <- meffil.qc(samplesheet, verbose = TRUE)

cat("✓ QC objects created for", length(qc_objects), "samples\n")

# ============================================================
# STEP 2: GENERATE QC SUMMARY AND IDENTIFY OUTLIERS
# ============================================================
cat("\n--- STEP 2: Generating QC summary and identifying outliers ---\n")

qc_summary <- meffil.qc.summary(qc_objects, parameters = qc_params, verbose = TRUE)

# Get outlier information from the QC summary
outlier_samples <- qc_summary$bad.samples
all_samples <- names(qc_objects)

cat("Outlier detection results:\n")
cat("  Total samples:", length(all_samples), "\n")
cat("  Outliers detected:", nrow(outlier_samples), "\n")

if (nrow(outlier_samples) > 0) {
  cat("\nOutlier samples:\n")
  print(outlier_samples)
}

# ============================================================
# STEP 3: SAVE QC RESULTS
# ============================================================
cat("\n--- STEP 3: Saving QC results ---\n")

# Save QC objects for later use (as recommended by meffil)
qc_objects_file <- file.path(qc_out_dir, paste0(plate_id, "_qc_objects.rds"))
saveRDS(qc_objects, qc_objects_file)
cat("✓ QC objects saved to:", basename(qc_objects_file), "\n")

# Save QC summary
qc_summary_file <- file.path(qc_out_dir, paste0(plate_id, "_qc_summary.rds"))
saveRDS(qc_summary, qc_summary_file)
cat("✓ QC summary saved to:", basename(qc_summary_file), "\n")

# Save outlier information
outliers_file <- file.path(qc_out_dir, paste0(plate_id, "_outliers.rds"))
saveRDS(outlier_samples, outliers_file)
cat("✓ Outliers saved to:", basename(outliers_file), "\n")

# Create and save passed samples list
passed_samples <- setdiff(all_samples, outlier_samples$sample.name)
writeLines(passed_samples, file.path(qc_out_dir, paste0(plate_id, "_passed_samples.txt")))
cat("✓ Passed samples list saved\n")

# ============================================================
# STEP 4: GENERATE QC REPORT
# ============================================================
cat("\n--- STEP 4: Generating QC report ---\n")

tryCatch({
  meffil.qc.report(
    qc_summary,
    output.file = file.path(qc_out_dir, paste0(plate_id, "_qc_report.html")),
    author = "Meffil Pipeline",
    study = paste("Plate QC:", plate_id)
  )
  cat("✓ HTML QC report generated\n")
}, error = function(e) {
  cat("Warning: Could not generate HTML report:", conditionMessage(e), "\n")
})

# ============================================================
# STEP 5: CREATE SIMPLE QC METRICS TABLE
# ============================================================
cat("\n--- STEP 5: Creating QC metrics table ---\n")

# Extract basic metrics - handle vectors properly
qc_metrics <- data.frame(
  Sample = names(qc_objects),
  Plate = plate_id,
  stringsAsFactors = FALSE
)

# Add call rate (handle potential vectors)
qc_metrics$CallRate <- sapply(names(qc_objects), function(sample_name) {
  x <- qc_objects[[sample_name]]
  if (is.null(x$bad.probes.detectionp)) return(NA)
  # If it's a vector, take the mean; if scalar, use as-is
  if (length(x$bad.probes.detectionp) > 1) {
    return(mean(1 - x$bad.probes.detectionp, na.rm = TRUE))
  } else {
    return(1 - x$bad.probes.detectionp)
  }
})

# Add predicted sex
qc_metrics$PredictedSex <- sapply(names(qc_objects), function(sample_name) {
  x <- qc_objects[[sample_name]]
  if (is.null(x$predicted.sex)) return(NA)
  return(as.character(x$predicted.sex))
})

# Add outlier status
qc_metrics$IsOutlier <- qc_metrics$Sample %in% outlier_samples$sample.name
qc_metrics$QC_Status <- ifelse(qc_metrics$IsOutlier, "FAIL", "PASS")

# Add failure reasons
qc_metrics$FailureReason <- sapply(qc_metrics$Sample, function(s) {
  if (s %in% outlier_samples$sample.name) {
    idx <- which(outlier_samples$sample.name == s)
    if (length(idx) > 0) {
      paste(outlier_samples$issue[idx], collapse = "; ")
    } else {
      "Unknown"
    }
  } else {
    "PASS"
  }
})

write.csv(qc_metrics, file.path(qc_out_dir, paste0(plate_id, "_qc_metrics.csv")), row.names = FALSE)
cat("✓ QC metrics table saved\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("PLATE QC COMPLETED:", plate_id, "\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nResults:\n")
cat("  Total samples:", nrow(qc_metrics), "\n")
cat("  PASS:", sum(qc_metrics$QC_Status == "PASS"), "\n")
cat("  FAIL:", sum(qc_metrics$QC_Status == "FAIL"), "\n")
cat("  Output directory:", qc_out_dir, "\n")

cat("\nOutput files:\n")
cat("  -", paste0(plate_id, "_qc_objects.rds\n"))
cat("  -", paste0(plate_id, "_qc_summary.rds\n"))
cat("  -", paste0(plate_id, "_outliers.rds\n"))
cat("  -", paste0(plate_id, "_passed_samples.txt\n"))
cat("  -", paste0(plate_id, "_qc_metrics.csv\n"))
cat("  -", paste0(plate_id, "_qc_report.html\n"))

cat("\n", rep("=", 60), "\n", sep = "")
