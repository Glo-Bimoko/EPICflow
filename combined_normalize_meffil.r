#!/usr/bin/env Rscript
# combined_normalize_meffil.R
# Normalize all plates together (combined normalization)
# This provides a single unified normalization across all samples

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript combined_normalize_meffil.R <qc_objects_dir> <passed_samples_dir> <combined_samplesheet> <out_dir> <qc_threshold>")
}

qc_objects_dir <- args[1]
passed_samples_dir <- args[2]
combined_samplesheet_file <- args[3]
out_dir <- args[4]
qc_threshold <- as.numeric(args[5])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("COMBINED NORMALIZATION (All Plates Together)\n")
cat(rep("=", 60), "\n\n")
cat("QC threshold:", qc_threshold, "\n")

# CRITICAL: Set single-threaded mode BEFORE loading any packages
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
# STEP 1: LOAD ALL QC OBJECTS AND PASSED SAMPLES
# ============================================================
cat("--- STEP 1: Loading data ---\n")

# Load combined samplesheet
cat("Loading combined samplesheet:", basename(combined_samplesheet_file), "\n")
combined_samplesheet <- read.csv(combined_samplesheet_file, stringsAsFactors = FALSE)
cat("Combined samplesheet has", nrow(combined_samplesheet), "samples\n")

# Find all QC object files
qc_object_files <- list.files(qc_objects_dir, 
                               pattern = "_qc_objects\\.rds$",
                               full.names = TRUE,
                               recursive = FALSE)

if (length(qc_object_files) == 0) {
  stop("No QC object files found in ", qc_objects_dir)
}

cat("Found", length(qc_object_files), "plate(s)\n")

# Load all QC objects and combine
all_qc_objects <- list()
plate_info <- data.frame(
  Sample = character(),
  Plate = character(),
  stringsAsFactors = FALSE
)

for (qc_file in qc_object_files) {
  plate_id <- sub("_qc_objects\\.rds$", "", basename(qc_file))
  cat("  Loading", plate_id, "...\n")
  
  plate_qc <- readRDS(qc_file)
  
  # Add to combined list
  for (sample_name in names(plate_qc)) {
    all_qc_objects[[sample_name]] <- plate_qc[[sample_name]]
    
    plate_info <- rbind(plate_info, data.frame(
      Sample = sample_name,
      Plate = plate_id,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\nTotal samples loaded:", length(all_qc_objects), "\n")

# ============================================================
# STEP 2: LOAD PASSED SAMPLES FROM ALL PLATES
# ============================================================
cat("\n--- STEP 2: Loading passed samples ---\n")

passed_sample_files <- list.files(passed_samples_dir,
                                  pattern = "_passed_samples\\.txt$",
                                  full.names = TRUE,
                                  recursive = FALSE)

all_passed_samples <- character()

for (passed_file in passed_sample_files) {
  plate_id <- sub("_passed_samples\\.txt$", "", basename(passed_file))
  samples <- readLines(passed_file)
  all_passed_samples <- c(all_passed_samples, samples)
  cat("  ", plate_id, ":", length(samples), "passed samples\n")
}

cat("\nTotal passed samples:", length(all_passed_samples), "\n")

# Filter QC objects to only include passed samples
qc_objects_passed <- all_qc_objects[names(all_qc_objects) %in% all_passed_samples]
cat("QC objects for normalization:", length(qc_objects_passed), "\n")

if (length(qc_objects_passed) == 0) {
  stop("No samples passed QC for normalization")
}

# ============================================================
# STEP 3: COMBINED FUNCTIONAL NORMALIZATION
# ============================================================
cat("\n--- STEP 3: Combined functional normalization ---\n")
cat("Normalizing all", length(qc_objects_passed), "samples together\n")

# Determine number of PCs
n_pcs <- min(10, max(2, floor(length(qc_objects_passed) / 10)))
cat("Using", n_pcs, "principal components\n")

# Perform functional normalization on all samples together
cat("Computing normalization parameters...\n")

tryCatch({
  norm_params <- meffil.normalize.quantiles(
    qc.objects = qc_objects_passed,
    number.pcs = n_pcs,
    verbose = TRUE
  )
  
  cat("Normalization parameters created successfully\n")
  
  cat("Applying normalization...\n")
  beta_matrix <- meffil.normalize.samples(
    norm_params,
    just.beta = TRUE,
    verbose = TRUE
  )
  
  cat("Beta normalization complete\n")
  cat("  Dimensions:", paste(dim(beta_matrix), collapse=" x "), "\n")
  
  # Get M values
  cat("Computing M values...\n")
  M_matrix <- meffil.normalize.samples(
    norm_params,
    just.beta = FALSE,
    verbose = FALSE
  )
  
  # Handle different return types
  if (is.list(M_matrix) && "M" %in% names(M_matrix)) {
    M_matrix <- M_matrix$M
  } else if (is.list(M_matrix) && "beta" %in% names(M_matrix)) {
    M_matrix <- log2(M_matrix$beta / (1 - M_matrix$beta))
  } else if (!is.matrix(M_matrix)) {
    M_matrix <- log2(beta_matrix / (1 - beta_matrix))
  }
  
  norm_matrix <- list(
    beta = beta_matrix,
    M = M_matrix
  )
  
  if (is.null(norm_matrix$beta) || ncol(norm_matrix$beta) == 0) {
    stop("Normalization failed: no samples in output matrix!")
  }
  
  cat("✓ Combined normalization completed\n")
  cat("  Samples:", ncol(norm_matrix$beta), "\n")
  cat("  Probes:", nrow(norm_matrix$beta), "\n")
  
}, error = function(e) {
  cat("\n✗ ERROR during normalization:", conditionMessage(e), "\n\n")
  
  if (grepl("pthread", conditionMessage(e), ignore.case = TRUE)) {
    cat("This is the pthread_create error!\n")
    cat("\nSOLUTION: Reinstall preprocessCore without threading:\n")
    cat("  remove.packages('preprocessCore')\n")
    cat("  BiocManager::install('preprocessCore', configure.args='--disable-threading')\n\n")
  }
  
  stop("Normalization failed. See error message above.")
})

# ============================================================
# STEP 4: SAVE COMBINED NORMALIZED DATA
# ============================================================
cat("\n--- STEP 4: Saving combined normalized data ---\n")

# Save combined data
write.csv(norm_matrix$beta, 
          file.path(out_dir, "BetaValues_all_plates_combined.csv"),
          row.names = TRUE)
cat("✓ Combined beta values saved\n")

write.csv(norm_matrix$M,
          file.path(out_dir, "MValues_all_plates_combined.csv"),
          row.names = TRUE)
cat("✓ Combined M values saved\n")

# Save normalization parameters
saveRDS(norm_params, file.path(out_dir, "combined_norm_params.rds"))
cat("✓ Normalization parameters saved\n")

# Add plate information to sample info
sample_info_combined <- plate_info[plate_info$Sample %in% colnames(norm_matrix$beta), ]
sample_info_combined <- sample_info_combined[match(colnames(norm_matrix$beta), 
                                                    sample_info_combined$Sample), ]

write.csv(sample_info_combined,
          file.path(out_dir, "sample_info_combined.csv"),
          row.names = FALSE)
cat("✓ Sample information saved\n")

# ============================================================
# STEP 5: SPLIT BY PLATE (optional for plate-level analysis)
# ============================================================
cat("\n--- STEP 5: Creating plate-specific subsets ---\n")

plates <- unique(sample_info_combined$Plate)
cat("Creating", length(plates), "plate-specific files\n")

plate_dir <- file.path(out_dir, "by_plate")
dir.create(plate_dir, showWarnings = FALSE, recursive = TRUE)

for (plate in plates) {
  plate_samples <- sample_info_combined$Sample[sample_info_combined$Plate == plate]
  
  beta_plate <- norm_matrix$beta[, plate_samples, drop = FALSE]
  M_plate <- norm_matrix$M[, plate_samples, drop = FALSE]
  
  write.csv(beta_plate,
            file.path(plate_dir, paste0(plate, "_beta.csv")),
            row.names = TRUE)
  
  write.csv(M_plate,
            file.path(plate_dir, paste0(plate, "_M.csv")),
            row.names = TRUE)
  
  cat("  ", plate, ":", ncol(beta_plate), "samples\n")
}

cat("✓ Plate-specific files saved to:", plate_dir, "\n")

# ============================================================
# STEP 6: GENERATE COMBINED QC REPORT (like individual plates)
# ============================================================
cat("\n--- STEP 6: Generating combined QC report ---\n")
cat("This will be the same style as individual plate reports\n")
cat("but treating all plates as one combined dataset\n\n")

tryCatch({
  cat("Creating combined QC summary from all QC objects...\n")
  
  # Set QC parameters (same as used in plate_qc_meffil.r)
  qc_threshold <- 0.95  # This should match your pipeline parameter
  qc_params <- meffil.qc.parameters(
    detectionp.samples.threshold = 1 - qc_threshold,
    detectionp.cpgs.threshold = 0.05,
    beadnum.samples.threshold = 0.1,
    beadnum.cpgs.threshold = 0.1,
    sex.outlier.sd = 5,
    snp.concordance.threshold = 0.9,
    sample.genotype.concordance.threshold = 0.8
  )
  
  cat("Using only passed samples for combined QC report...\n")
  cat("  QC objects available:", length(qc_objects_passed), "\n")
  
  # Generate QC summary for all passed samples together
  cat("Generating QC summary (this may take a few minutes)...\n")
  combined_qc_summary <- meffil.qc.summary(
    qc.objects = qc_objects_passed,
    parameters = qc_params,
    verbose = TRUE
  )
  
  cat("✓ Combined QC summary created\n")
  
  # Generate the combined QC report (same as individual plates)
  cat("Generating combined QC HTML report...\n")
  
  qc_report_file <- file.path(out_dir, "combined_qc_report.html")
  
  meffil.qc.report(
    qc.summary = combined_qc_summary,
    output.file = qc_report_file,
    author = "Meffil Combined Pipeline",
    study = paste0("Combined QC Report - All ", 
                   length(unique(sample_info_combined$Plate)), 
                   " Plates Together (",
                   length(qc_objects_passed), " samples)")
  )
  
  if (file.exists(qc_report_file)) {
    cat("✓ Combined QC report generated successfully!\n")
    cat("  Location:", qc_report_file, "\n")
    cat("  This report shows all samples as if they were one plate\n")
  } else {
    cat("⚠ QC report function completed but file not found\n")
  }
  
}, error = function(e) {
  cat("\n⚠ Combined QC report generation failed\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("\nDebug information:\n")
  cat("  Number of QC objects:", length(qc_objects_passed), "\n")
  cat("  QC object names (first 5):", paste(head(names(qc_objects_passed), 5), collapse=", "), "\n")
})

# ============================================================
# STEP 7: GENERATE NORMALIZATION REPORT
# ============================================================
cat("\n--- STEP 7: Generating normalization report ---\n")

tryCatch({
  cat("Creating normalization summary...\n")
  
  norm_summary <- meffil.normalization.summary(
    norm.objects = norm_params,
    pcs = n_pcs
  )
  
  cat("✓ Normalization summary created\n")
  cat("Generating normalization HTML report...\n")
  
  norm_report_file <- file.path(out_dir, "combined_normalization_report.html")
  
  meffil.normalization.report(
    normalization.summary = norm_summary,
    output.file = norm_report_file,
    author = "Meffil Pipeline",
    study = paste0("Combined Normalization - ", 
                   length(unique(sample_info_combined$Plate)), " plates, ",
                   nrow(sample_info_combined), " samples")
  )
  
  if (file.exists(norm_report_file)) {
    cat("✓ Normalization report generated\n")
  }
  
}, error = function(e) {
  cat("⚠ Normalization report failed:", conditionMessage(e), "\n")
})

# Create a comprehensive HTML summary (always, as backup or main report)
cat("\nCreating comprehensive HTML summary...\n")
tryCatch({
  # Calculate some summary statistics
  mean_call_rate <- mean(apply(norm_matrix$beta, 2, function(x) sum(!is.na(x))/length(x)))
  
  html_content <- paste0(
    "<!DOCTYPE html>\n<html>\n<head>\n",
    "<title>Combined Normalization Report</title>\n",
    "<meta charset='UTF-8'>\n",
    "<style>\n",
    "body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background-color: #f5f5f5; }\n",
    ".container { max-width: 1200px; margin: 0 auto; padding: 20px; }\n",
    ".header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }\n",
    "h1 { margin: 0; font-size: 2.5em; }\n",
    ".subtitle { opacity: 0.9; margin-top: 10px; }\n",
    ".card { background: white; padding: 25px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); margin-bottom: 20px; }\n",
    "h2 { color: #667eea; border-bottom: 3px solid #667eea; padding-bottom: 10px; }\n",
    ".stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }\n",
    ".stat-box { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; }\n",
    ".stat-value { font-size: 2.5em; font-weight: bold; margin: 10px 0; }\n",
    ".stat-label { opacity: 0.9; font-size: 0.9em; }\n",
    "table { width: 100%; border-collapse: collapse; margin: 20px 0; }\n",
    "th { background-color: #667eea; color: white; padding: 15px; text-align: left; }\n",
    "td { padding: 12px; border-bottom: 1px solid #e0e0e0; }\n",
    "tr:hover { background-color: #f8f9fa; }\n",
    ".visualization { text-align: center; margin: 30px 0; }\n",
    ".visualization img { max-width: 100%; border: 1px solid #ddd; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }\n",
    ".file-list { list-style: none; padding: 0; }\n",
    ".file-list li { padding: 10px; margin: 5px 0; background: #f8f9fa; border-left: 4px solid #667eea; }\n",
    ".success { color: #27ae60; font-weight: bold; }\n",
    ".info-box { background: #e3f2fd; border-left: 4px solid #2196f3; padding: 15px; margin: 20px 0; border-radius: 4px; }\n",
    "</style>\n</head>\n<body>\n",
    "<div class='container'>\n",
    "<div class='header'>\n",
    "<h1>🧬 Combined Normalization Report</h1>\n",
    "<div class='subtitle'>Meffil Functional Normalization - All Plates Combined</div>\n",
    "</div>\n"
  )
  
  # Summary statistics
  html_content <- paste0(html_content,
    "<div class='card'>\n<h2>📊 Summary Statistics</h2>\n",
    "<div class='stat-grid'>\n",
    "<div class='stat-box'>\n<div class='stat-value'>", ncol(norm_matrix$beta), "</div>\n<div class='stat-label'>Total Samples</div>\n</div>\n",
    "<div class='stat-box'>\n<div class='stat-value'>", nrow(norm_matrix$beta), "</div>\n<div class='stat-label'>CpG Probes</div>\n</div>\n",
    "<div class='stat-box'>\n<div class='stat-value'>", length(unique(sample_info_combined$Plate)), "</div>\n<div class='stat-label'>Plates</div>\n</div>\n",
    "<div class='stat-box'>\n<div class='stat-value'>", n_pcs, "</div>\n<div class='stat-label'>PCs Used</div>\n</div>\n",
    "</div>\n",
    "<div class='info-box'>\n",
    "<strong>Normalization Method:</strong> Meffil functional normalization with ", n_pcs, " principal components<br>\n",
    "<strong>Mean Probe Detection:</strong> ", round(mean_call_rate * 100, 2), "%\n",
    "</div>\n</div>\n"
  )
  
  # Samples by plate
  html_content <- paste0(html_content,
    "<div class='card'>\n<h2>🔬 Samples by Plate</h2>\n<table>\n",
    "<tr><th>Plate ID</th><th>Number of Samples</th><th>Percentage</th></tr>\n"
  )
  
  plate_counts <- table(sample_info_combined$Plate)
  total_samples <- sum(plate_counts)
  for (plate in names(sort(plate_counts, decreasing = TRUE))) {
    pct <- round(100 * plate_counts[plate] / total_samples, 1)
    html_content <- paste0(html_content,
      "<tr><td><strong>", plate, "</strong></td><td>", plate_counts[plate], "</td><td>", pct, "%</td></tr>\n"
    )
  }
  html_content <- paste0(html_content, "</table>\n</div>\n")
  
  # PCA visualization
  if (file.exists(file.path(out_dir, "PCA_combined_normalization.png"))) {
    html_content <- paste0(html_content,
      "<div class='card'>\n<h2>📈 Principal Component Analysis</h2>\n",
      "<p>Visualization showing sample clustering and plate effects after combined normalization:</p>\n",
      "<div class='visualization'>\n",
      "<img src='PCA_combined_normalization.png' alt='PCA Plot'>\n",
      "</div>\n</div>\n"
    )
  }
  
  # Output files
  html_content <- paste0(html_content,
    "<div class='card'>\n<h2>📁 Output Files</h2>\n",
    "<ul class='file-list'>\n",
    "<li><strong>BetaValues_all_plates_combined.csv</strong> - Beta values (0-1 scale) for all samples</li>\n",
    "<li><strong>MValues_all_plates_combined.csv</strong> - M values (log2 ratio) for all samples</li>\n",
    "<li><strong>sample_info_combined.csv</strong> - Sample metadata with plate assignments</li>\n",
    "<li><strong>combined_norm_params.rds</strong> - Normalization parameters for reproducibility</li>\n",
    "<li><strong>by_plate/</strong> - Directory containing plate-specific subsets:\n",
    "<ul style='margin-left: 20px; margin-top: 10px;'>\n"
  )
  
  for (plate in names(plate_counts)) {
    html_content <- paste0(html_content,
      "<li>", plate, "_beta.csv (", plate_counts[plate], " samples)</li>\n",
      "<li>", plate, "_M.csv (", plate_counts[plate], " samples)</li>\n"
    )
  }
  
  html_content <- paste0(html_content,
    "</ul></li>\n</ul>\n</div>\n",
    "<div class='card'>\n<h2>✅ Pipeline Status</h2>\n",
    "<p class='success'>✓ Combined normalization completed successfully</p>\n",
    "<p>All plates have been normalized together using a unified reference, ensuring comparability across the entire dataset.</p>\n",
    "</div>\n",
    "</div>\n</body>\n</html>"
  )
  
  writeLines(html_content, file.path(out_dir, "normalization_summary.html"))
  cat("✓ Comprehensive HTML summary created: normalization_summary.html\n")
  
}, error = function(e) {
  cat("⚠ Could not create HTML summary:", conditionMessage(e), "\n")
})

# ============================================================
# STEP 7: PCA TO VISUALIZE PLATE EFFECTS
# ============================================================
cat("\n--- STEP 7: Creating PCA visualization ---\n")

beta_complete <- norm_matrix$beta[complete.cases(norm_matrix$beta), ]
probe_vars <- apply(beta_complete, 1, var)
top_probes <- names(sort(probe_vars, decreasing = TRUE)[1:min(10000, length(probe_vars))])
beta_for_pca <- beta_complete[top_probes, ]

pca_result <- prcomp(t(beta_for_pca), scale. = TRUE)
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(beta_for_pca),
  stringsAsFactors = FALSE
)

# Add plate info
pca_df <- merge(pca_df, sample_info_combined, by = "Sample")

# Test plate effect
if (length(unique(pca_df$Plate)) > 1) {
  plate_pval <- anova(lm(PC1 ~ Plate, data = pca_df))$'Pr(>F)'[1]
  cat("Plate effect test (PC1 ~ Plate): p =", format.pval(plate_pval, digits = 3), "\n")
  
  # Plot PCA
  png(file.path(out_dir, "PCA_combined_normalization.png"),
      width = 1600, height = 1000, res = 150)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate, shape = Plate)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(title = "PCA: Combined Normalization (All Plates Together)",
         subtitle = paste("Plate effect p-value =", format.pval(plate_pval, digits = 3)),
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)"))
  print(p)
  dev.off()
  
  cat("✓ PCA plot saved\n")
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("COMBINED NORMALIZATION COMPLETED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nSummary:\n")
cat("  Total samples normalized:", ncol(norm_matrix$beta), "\n")
cat("  Total probes:", nrow(norm_matrix$beta), "\n")
cat("  Plates included:", length(unique(sample_info_combined$Plate)), "\n")
cat("  Method: Meffil functional normalization (combined)\n")
cat("  PCs used:", n_pcs, "\n")

cat("\nOutput files:\n")
cat("  - BetaValues_all_plates_combined.csv (ALL samples)\n")
cat("  - MValues_all_plates_combined.csv (ALL samples)\n")
cat("  - sample_info_combined.csv\n")
cat("  - combined_norm_params.rds\n")
cat("  - combined_normalization_report.html\n")
cat("  - by_plate/ (individual plate subsets)\n")

if (length(unique(sample_info_combined$Plate)) > 1) {
  cat("\nPlate effects in combined normalization:\n")
  cat("  p-value:", format.pval(plate_pval, digits = 3), "\n")
  if (plate_pval < 0.05) {
    cat("  Note: Significant plate effects remain after combined normalization\n")
    cat("        Consider using ComBat for additional batch correction\n")
  } else {
    cat("  ✓ No significant plate effects\n")
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
