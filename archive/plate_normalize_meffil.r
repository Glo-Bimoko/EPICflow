#!/usr/bin/env Rscript
# plate_normalize_meffil.R
# Normalization following meffil best practices
# With proper preprocessCore threading fix

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plate_normalize_meffil.R <samplesheet.csv> <qc_objects.rds> <passed_samples.txt> <out_dir>")
}

samplesheet_file <- args[1]
qc_objects_file <- args[2]
passed_samples_file <- args[3]
out_dir <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plate_id <- sub("_samplesheet\\.csv$", "", basename(samplesheet_file))

cat("\n", rep("=", 60), "\n", sep = "")
cat("PLATE NORMALIZATION (Meffil Best Practices):", plate_id, "\n")
cat(rep("=", 60), "\n\n")

# CRITICAL: Set single-threaded mode BEFORE loading any packages
# This is the key to avoiding preprocessCore pthread errors
# See: https://github.com/perishky/meffil/wiki/Common-problems
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(NUMEXPR_NUM_THREADS = 1)
Sys.setenv(OMP_THREAD_LIMIT = 1)

cat("Threading environment set to single-threaded mode\n")
options(mc.cores = 1)

suppressPackageStartupMessages({
  library(meffil)
  library(ggplot2)
})

# ============================================================
# STEP 1: LOAD DATA
# ============================================================
cat("--- STEP 1: Loading data ---\n")

# Read passed samples
passed_samples <- readLines(passed_samples_file)
cat("Passed QC samples:", length(passed_samples), "\n")

if (length(passed_samples) == 0) {
  stop("No samples passed QC for normalization")
}

# Load QC objects
qc_objects <- readRDS(qc_objects_file)
cat("Loaded QC objects for", length(qc_objects), "samples\n")

# Filter QC objects to only include passed samples
qc_objects_passed <- qc_objects[names(qc_objects) %in% passed_samples]
cat("Samples for normalization:", length(qc_objects_passed), "\n")

# ============================================================
# STEP 2: FUNCTIONAL NORMALIZATION (Meffil recommended)
# ============================================================
cat("\n--- STEP 2: Functional normalization ---\n")

# Determine number of PCs to use (meffil recommendation: 2-10)
n_pcs <- min(10, max(2, floor(length(qc_objects_passed) / 10)))
cat("Using", n_pcs, "principal components\n")

# Perform functional normalization
cat("Computing normalization parameters...\n")

tryCatch({
  norm_params <- meffil.normalize.quantiles(
    qc.objects = qc_objects_passed,
    number.pcs = n_pcs,
    verbose = TRUE
  )
  
  cat("Normalization parameters created successfully\n")
  cat("  Class:", class(norm_params), "\n")
  cat("  Length:", length(norm_params), "\n")
  
  cat("Applying normalization...\n")
  beta_matrix <- meffil.normalize.samples(
    norm_params,
    just.beta = TRUE,  # Get beta values first
    verbose = TRUE
  )
  
  cat("Beta normalization complete\n")
  cat("  Class:", class(beta_matrix), "\n")
  cat("  Dimensions:", paste(dim(beta_matrix), collapse=" x "), "\n")
  
  # Get M values separately
  cat("Computing M values...\n")
  M_matrix <- meffil.normalize.samples(
    norm_params,
    just.beta = FALSE,  # This actually returns M values when combined properly
    verbose = FALSE
  )
  
  # Handle different return types from meffil
  if (is.list(M_matrix) && "M" %in% names(M_matrix)) {
    M_matrix <- M_matrix$M
  } else if (is.list(M_matrix) && "beta" %in% names(M_matrix)) {
    # Calculate M values from beta if needed
    M_matrix <- log2(M_matrix$beta / (1 - M_matrix$beta))
  } else if (!is.matrix(M_matrix)) {
    # Calculate M from beta
    M_matrix <- log2(beta_matrix / (1 - beta_matrix))
  }
  
  # Create output structure
  norm_matrix <- list(
    beta = beta_matrix,
    M = M_matrix
  )
  
  # Validate
  if (is.null(norm_matrix$beta) || ncol(norm_matrix$beta) == 0) {
    cat("ERROR: Beta matrix has no samples!\n")
    stop("Normalization failed: no samples in output matrix!")
  }
  
  cat("✓ Normalization completed\n")
  cat("  Samples:", ncol(norm_matrix$beta), "\n")
  cat("  Probes:", nrow(norm_matrix$beta), "\n")
  
}, error = function(e) {
  cat("\n✗ ERROR during normalization:", conditionMessage(e), "\n\n")
  
  if (grepl("pthread", conditionMessage(e), ignore.case = TRUE)) {
    cat("This is the pthread_create error!\n")
    cat("\nSOLUTION: You need to reinstall preprocessCore without threading:\n")
    cat("  1. In R console, run:\n")
    cat("     remove.packages('preprocessCore')\n")
    cat("     BiocManager::install('preprocessCore', configure.args='--disable-threading')\n")
    cat("  2. Or run the fix script: bash fix_meffil_threading.sh\n")
    cat("\nSee: https://github.com/perishky/meffil/wiki/Common-problems\n\n")
  }
  
  stop("Normalization failed. See error message above.")
})

# ============================================================
# STEP 3: SAVE NORMALIZED DATA
# ============================================================
cat("\n--- STEP 3: Saving normalized data ---\n")

# Save normalized values
write.csv(norm_matrix$beta, 
          file.path(out_dir, paste0(plate_id, "_beta.csv")),
          row.names = TRUE)
cat("✓ Beta values saved\n")

write.csv(norm_matrix$M,
          file.path(out_dir, paste0(plate_id, "_M.csv")),
          row.names = TRUE)
cat("✓ M values saved\n")

# Save normalization parameters for later use
saveRDS(norm_params, file.path(out_dir, paste0(plate_id, "_norm_params.rds")))
cat("✓ Normalization parameters saved\n")

# Save sample info
if (!is.null(norm_matrix$beta) && ncol(norm_matrix$beta) > 0) {
  sample_info <- data.frame(
    Sample = colnames(norm_matrix$beta),
    Plate = plate_id,
    stringsAsFactors = FALSE
  )
  write.csv(sample_info,
            file.path(out_dir, paste0(plate_id, "_sample_info.csv")),
            row.names = FALSE)
  cat("✓ Sample information saved\n")
} else {
  cat("⚠ Warning: Could not save sample info (no samples in matrix)\n")
}

# ============================================================
# STEP 4: GENERATE NORMALIZATION REPORT
# ============================================================
cat("\n--- STEP 4: Generating normalization report ---\n")

tryCatch({
  norm_summary <- meffil.normalization.summary(norm_params, verbose = FALSE)
  meffil.normalization.report(
    norm_summary,
    output.file = file.path(out_dir, paste0(plate_id, "_norm_report.html")),
    author = "Meffil Pipeline",
    study = paste("Plate", plate_id)
  )
  cat("✓ Normalization report generated\n")
}, error = function(e) {
  cat("Warning: Could not generate normalization report:", conditionMessage(e), "\n")
})

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("NORMALIZATION COMPLETED:", plate_id, "\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nSummary:\n")
cat("  Samples normalized:", ncol(norm_matrix$beta), "\n")
cat("  Probes:", nrow(norm_matrix$beta), "\n")
cat("  Method: Meffil functional normalization\n")
cat("  PCs used:", n_pcs, "\n")

cat("\nOutput files:\n")
cat("  -", paste0(plate_id, "_beta.csv\n"))
cat("  -", paste0(plate_id, "_M.csv\n"))
cat("  -", paste0(plate_id, "_norm_params.rds\n"))
cat("  -", paste0(plate_id, "_sample_info.csv\n"))
cat("  -", paste0(plate_id, "_norm_report.html\n"))

cat("\n", rep("=", 60), "\n", sep = "")
