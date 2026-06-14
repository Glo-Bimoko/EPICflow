#!/usr/bin/env Rscript
# bin/plate_qc_meffil.r
# QC for an entire plate following meffil best practices.
#
# References:
#   https://github.com/perishky/meffil/wiki/Sample-QC
#   https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets
#
# Usage:
#   Rscript plate_qc_meffil.r <samplesheet.csv> <qc_out_dir> <qc_threshold> [max_gb]
#
# Arguments:
#   samplesheet.csv  : meffil-format samplesheet (Sample_Name, Sex, Basename, ...)
#   qc_out_dir       : directory for output files
#   qc_threshold     : call-rate threshold for the pipeline's PASS/FAIL decision
#                      (e.g. 0.95 means a sample must have ≥95 % detected probes).
#                      This is SEPARATE from meffil's internal QC parameters.
#   max_gb           : (optional) total RAM available in GB (default 8).
#                      Used only to compute max.bytes per mclapply fork inside
#                      meffil.qc(); does NOT affect which samples pass QC.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plate_qc_meffil.r <samplesheet.csv> <qc_out_dir> <qc_threshold> [max_gb]")
}

samplesheet_file <- args[1]
qc_out_dir       <- args[2]
qc_threshold     <- as.numeric(args[3])
max_gb           <- if (length(args) >= 4) as.numeric(args[4]) else 8

dir.create(qc_out_dir, showWarnings = FALSE, recursive = TRUE)

plate_id <- sub("_samplesheet\\.csv$", "", basename(samplesheet_file))

cat("\n", rep("=", 60), "\n", sep = "")
cat("PLATE QC (Meffil Best Practices):", plate_id, "\n")
cat(rep("=", 60), "\n\n")

# ── Threading ──────────────────────────────────────────────────────────────
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                      Sys.getenv("PBS_NUM_PPN", 4)))
cat("Cores           :", n_cores, "\n")
cat("Memory ceiling  :", max_gb, "GB\n")
options(mc.cores = n_cores)

# FIX 1: use as.numeric() (64-bit double) not as.integer() (32-bit signed).
# as.integer() overflows to NA for values > 2,147,483,647 (~2.1 GB), which
# happens whenever (max_gb * 0.75 * 1024^3) / n_cores exceeds that ceiling —
# e.g. max_gb=28, n_cores=4 yields ~5.6 GB → NA → meffil.qc() crashes with
# "Error in if (n.fun < 1): missing value where TRUE/FALSE needed"
max_bytes <- as.numeric((max_gb * 0.75 * 1024^3) / n_cores)

if (is.na(max_bytes) || max_bytes < 1e6) {
  warning(sprintf(
    "max_bytes calculation produced an invalid value (%s). Falling back to 1 GB per fork.",
    max_bytes
  ))
  max_bytes <- 1e9
}

cat("max.bytes/fork  :", format(max_bytes, big.mark = ",", scientific = FALSE), "bytes\n\n")

# ── Load meffil ────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(meffil)
  library(ggplot2)
})

# ============================================================
# STEP 0: READ SAMPLESHEET
# ============================================================
cat("--- STEP 0: Reading samplesheet ---\n")

samplesheet <- read.csv(samplesheet_file, stringsAsFactors = FALSE)
cat("Found", nrow(samplesheet), "samples in samplesheet\n")

required_cols <- c("Sample_Name", "Basename")
missing_cols  <- setdiff(required_cols, colnames(samplesheet))
if (length(missing_cols) > 0)
  stop("Samplesheet missing required columns: ", paste(missing_cols, collapse = ", "))

# ============================================================
# STEP 1: CREATE QC OBJECTS
# ============================================================
cat("\n--- STEP 1: Creating QC objects ---\n")
cat("Creating QC objects for", nrow(samplesheet), "samples...\n")

qc_objects <- meffil.qc(samplesheet, verbose = TRUE, max.bytes = max_bytes)

cat("✓ QC objects created for", length(qc_objects), "samples\n")

# Verify objects are named — meffil names the list by Sample_Name
if (is.null(names(qc_objects)) || any(names(qc_objects) == "")) {
  # Fall back: name by the sample.name field stored inside each object
  names(qc_objects) <- sapply(qc_objects, function(x) x$sample.name)
  cat("NOTE: qc_objects were unnamed; named from internal sample.name field\n")
}

# ============================================================
# STEP 2: QC SUMMARY AND OUTLIER DETECTION
# ============================================================
cat("\n--- STEP 2: Generating QC summary and identifying outliers ---\n")

qc_params <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1,
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

qc_summary <- meffil.qc.summary(qc_objects, parameters = qc_params, verbose = TRUE)

all_bad     <- qc_summary$bad.samples
all_samples <- names(qc_objects)

cat("Outlier detection results:\n")
cat("  Total samples     :", length(all_samples), "\n")
cat("  Flagged (all)     :", nrow(all_bad), "\n")

if (nrow(all_bad) > 0) {
  cat("\nAll flagged samples (before issue filter):\n")
  print(all_bad)
  cat("\nIssue counts:\n")
  print(table(all_bad$issue))
}

# ============================================================
# STEP 3: DETERMINE PASSED SAMPLES
# ============================================================
cat("\n--- STEP 3: Determining passed samples ---\n")

EXCLUDE_ISSUES <- c(
  "Control probe (dye.bias)",
  "Methylated vs Unmethylated",
  "X-Y ratio outlier",
  "Low bead numbers",
  "Detection p-value",
  "Sex mismatch",
  "Genotype mismatch",
  "Control probe (bisulfite1)",
  "Control probe (bisulfite2)"
)

if (nrow(all_bad) > 0) {
  outlier_filtered <- all_bad[all_bad$issue %in% EXCLUDE_ISSUES, ]
} else {
  outlier_filtered <- all_bad
}

cat("Samples excluded by issue filter:", nrow(outlier_filtered), "\n")
if (nrow(outlier_filtered) > 0) print(outlier_filtered)

samples_to_exclude <- unique(outlier_filtered$sample.name)

# ── Call-rate computation ──────────────────────────────────────────────────
# FIX 2: bad.probes.detectionp is a NAMED NUMERIC VECTOR of p-values for the
# probes that FAILED detection (i.e. p > threshold) — it is NOT a logical
# vector over all probes.  The previous code did mean(dp) which computed the
# mean p-value of failing probes (≈ 0.5), yielding call_rate ≈ 0.5 for every
# sample regardless of actual quality, and failing the 0.95 threshold for all.
#
# Correct formula (mirrors meffil.plot.detectionp.samples in qc-report.r):
#   call_rate = 1 - (number of failed probes / total probes for the featureset)
# Y-chromosome probes are excluded for female samples, matching meffil's own
# plot function which does setdiff(probes, y.probes) for females.
# ──────────────────────────────────────────────────────────────────────────

featureset <- tryCatch(qc_objects[[1]]$featureset, error = function(e) "epic")
if (is.null(featureset)) featureset <- "epic"

all_probe_names <- meffil.get.features(featureset)$name
y_probe_names   <- meffil.get.y.sites(featureset)
n_total_probes  <- length(all_probe_names)
n_non_y_probes  <- length(setdiff(all_probe_names, y_probe_names))

cat("\nFeatureset        :", featureset, "\n")
cat("Total probes      :", n_total_probes, "\n")
cat("Non-Y probes      :", n_non_y_probes, "(used for female call-rate denominator)\n\n")

call_rate <- sapply(names(qc_objects), function(s) {
  obj          <- qc_objects[[s]]
  bad_dp       <- obj$bad.probes.detectionp   # named numeric vector of failing p-values
  predicted_sex <- tryCatch(as.character(obj$predicted.sex), error = function(e) "NA")

  if (is.null(bad_dp)) {
    # No bad probes at all — perfect call rate
    return(1.0)
  }

  # Intersect with the featureset probes so the denominator is consistent
  bad_in_featureset <- intersect(names(bad_dp), all_probe_names)

  # For females, exclude Y-probe failures from both numerator and denominator
  if (!is.na(predicted_sex) && predicted_sex == "F") {
    bad_in_featureset <- setdiff(bad_in_featureset, y_probe_names)
    denom <- n_non_y_probes
  } else {
    denom <- n_total_probes
  }

  1 - length(bad_in_featureset) / denom
})

cat("Call-rate summary (n =", length(call_rate), "samples):\n")
cat("  Min  :", round(min(call_rate,  na.rm = TRUE), 4), "\n")
cat("  Mean :", round(mean(call_rate, na.rm = TRUE), 4), "\n")
cat("  Max  :", round(max(call_rate,  na.rm = TRUE), 4), "\n")

low_call_rate_samples <- names(call_rate)[!is.na(call_rate) & call_rate < qc_threshold]

if (length(low_call_rate_samples) > 0) {
  cat("\nSamples below call-rate threshold (", qc_threshold, "):\n", sep = "")
  for (s in low_call_rate_samples)
    cat(sprintf("  %s  call_rate=%.4f\n", s, call_rate[s]))
  samples_to_exclude <- union(samples_to_exclude, low_call_rate_samples)
}

passed_samples <- setdiff(all_samples, samples_to_exclude)
cat("\nFinal QC summary:\n")
cat("  Total samples   :", length(all_samples), "\n")
cat("  Excluded        :", length(samples_to_exclude), "\n")
cat("  Passed          :", length(passed_samples), "\n")

# ── Hard stop if nothing passed ────────────────────────────────────────────
# Exit with status 1 so Nextflow marks this task FAILED and does NOT cache
# the result.  Without this, a silent bug (e.g. wrong call-rate formula)
# produces an empty passed_samples.txt with exit 0, which Nextflow caches
# and reuses on every subsequent -resume run, making the bug invisible.
if (length(passed_samples) == 0) {
  stop(sprintf(paste(
    "No samples passed QC for plate %s.",
    "  call-rate threshold : %s",
    "  call-rate range     : %.4f – %.4f",
    "  meffil flagged      : %d samples",
    "  issue-filter removed: %d samples",
    "Check the QC report and metrics CSV for details.",
    sep = "\n  "
  ), plate_id, qc_threshold,
     min(call_rate, na.rm = TRUE), max(call_rate, na.rm = TRUE),
     nrow(all_bad), nrow(outlier_filtered)))
}

# ============================================================
# STEP 4: SAVE QC OBJECTS AND SUMMARIES
# ============================================================
cat("\n--- STEP 4: Saving QC results ---\n")

qc_objects_file <- file.path(qc_out_dir, paste0(plate_id, "_qc_objects.rds"))
saveRDS(qc_objects, qc_objects_file)
cat("✓ QC objects saved to:", basename(qc_objects_file), "\n")

saveRDS(qc_summary, file.path(qc_out_dir, paste0(plate_id, "_qc_summary.rds")))
cat("✓ QC summary saved\n")

saveRDS(all_bad, file.path(qc_out_dir, paste0(plate_id, "_outliers.rds")))
cat("✓ Flagged samples saved\n")

writeLines(passed_samples,
           file.path(qc_out_dir, paste0(plate_id, "_passed_samples.txt")))
cat("✓ Passed samples list saved (", length(passed_samples), "samples)\n")

# ============================================================
# STEP 5: HTML QC REPORT
# ============================================================
cat("\n--- STEP 5: Generating QC report ---\n")

tryCatch({
  meffil.qc.report(
    qc_summary,
    output.file = file.path(qc_out_dir, paste0(plate_id, "_qc_report.html")),
    author = "Meffil Pipeline",
    study  = paste("Plate QC:", plate_id)
  )
  cat("✓ HTML QC report generated\n")
}, error = function(e) {
  cat("Warning: Could not generate HTML report:", conditionMessage(e), "\n")
})

# ============================================================
# STEP 6: QC METRICS TABLE
# ============================================================
cat("\n--- STEP 6: Creating QC metrics table ---\n")

qc_metrics <- data.frame(
  Sample           = all_samples,
  Plate            = plate_id,
  CallRate         = call_rate[all_samples],
  PredictedSex     = sapply(all_samples, function(s) {
    ps <- tryCatch(as.character(qc_objects[[s]]$predicted.sex), error = function(e) NA_character_)
    if (is.null(ps)) NA_character_ else ps
  }),
  Meffil_Flagged   = all_samples %in% all_bad$sample.name,
  Pipeline_Exclude = all_samples %in% samples_to_exclude,
  QC_Status        = ifelse(all_samples %in% samples_to_exclude, "FAIL", "PASS"),
  stringsAsFactors = FALSE
)

qc_metrics$MeffIssues <- sapply(all_samples, function(s) {
  if (!s %in% all_bad$sample.name) return("")
  paste(all_bad$issue[all_bad$sample.name == s], collapse = "; ")
})

write.csv(qc_metrics,
          file.path(qc_out_dir, paste0(plate_id, "_qc_metrics.csv")),
          row.names = FALSE)
cat("✓ QC metrics table saved\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("PLATE QC COMPLETED:", plate_id, "\n")
cat(rep("=", 60), "\n")
cat("  Total samples   :", nrow(qc_metrics), "\n")
cat("  PASS            :", sum(qc_metrics$QC_Status == "PASS"), "\n")
cat("  FAIL            :", sum(qc_metrics$QC_Status == "FAIL"), "\n")
cat("  Call-rate thresh:", qc_threshold, "\n")
cat("  Output dir      :", qc_out_dir, "\n\n")
cat("Output files:\n")
for (f in list.files(qc_out_dir, pattern = paste0("^", plate_id)))
  cat("  -", f, "\n")
cat("\n")
