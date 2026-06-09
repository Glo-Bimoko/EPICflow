#!/usr/bin/env Rscript
# multi_plate_qc_report.R
# Generate comprehensive QC report across all plates
# Identifies true outliers when considering all samples together

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript multi_plate_qc_report.r <qc_results_dir> <out_dir> <qc_threshold>")
}

qc_results_dir <- args[1]
out_dir <- args[2]
qc_threshold <- as.numeric(args[3])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("MULTI-PLATE QC REPORT GENERATION\n")
cat(rep("=", 60), "\n\n")

suppressPackageStartupMessages({
  library(meffil)
  library(ggplot2)
  library(gridExtra)
})

# ============================================================
# STEP 1: LOAD ALL PLATE QC OBJECTS
# ============================================================
cat("--- STEP 1: Loading QC objects from all plates ---\n")

# Find all QC object files in the input directory
qc_object_files <- list.files(qc_results_dir, 
                               pattern = "_qc_objects\\.rds$",
                               full.names = TRUE,
                               recursive = FALSE)

if (length(qc_object_files) == 0) {
  stop("No QC object files found in ", qc_results_dir)
}

cat("Found", length(qc_object_files), "plate(s) with QC objects\n")

# Load and combine all QC objects
all_qc_objects <- list()
plate_mapping <- data.frame(
  Sample = character(),
  Plate = character(),
  stringsAsFactors = FALSE
)

for (qc_file in qc_object_files) {
  plate_id <- sub("_qc_objects\\.rds$", "", basename(qc_file))
  cat("  Loading", plate_id, "...\n")
  
  plate_qc <- readRDS(qc_file)
  
  # Add to combined list with unique names
  for (sample_name in names(plate_qc)) {
    # Create unique key: plate_sample
    unique_key <- paste0(plate_id, "_", sample_name)
    all_qc_objects[[unique_key]] <- plate_qc[[sample_name]]
    
    # Track plate membership
    plate_mapping <- rbind(plate_mapping, data.frame(
      Sample = unique_key,
      Original_Sample = sample_name,
      Plate = plate_id,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\nTotal samples across all plates:", length(all_qc_objects), "\n")
cat("Plates:", length(unique(plate_mapping$Plate)), "\n")

# ============================================================
# STEP 2: GENERATE MULTI-PLATE QC SUMMARY
# ============================================================
cat("\n--- STEP 2: Generating multi-plate QC summary ---\n")

# Set QC parameters (same as individual plates)
qc_params <- meffil.qc.parameters(
  detectionp.samples.threshold = 1 - qc_threshold,
  detectionp.cpgs.threshold = 0.05,
  beadnum.samples.threshold = 0.1,
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,
  snp.concordance.threshold = 0.9,
  sample.genotype.concordance.threshold = 0.8
)

cat("Generating QC summary for all", length(all_qc_objects), "samples...\n")
multi_qc_summary <- meffil.qc.summary(all_qc_objects, 
                                       parameters = qc_params, 
                                       verbose = TRUE)

# Get outlier information
multi_outliers <- multi_qc_summary$bad.samples
all_samples <- names(all_qc_objects)

cat("\nMulti-plate outlier detection results:\n")
cat("  Total samples:", length(all_samples), "\n")
cat("  Outliers detected:", nrow(multi_outliers), "\n")
cat("  Pass rate:", round(100 * (1 - nrow(multi_outliers)/length(all_samples)), 2), "%\n")

if (nrow(multi_outliers) > 0) {
  cat("\nOutlier breakdown by plate:\n")
  outlier_by_plate <- merge(multi_outliers, plate_mapping, 
                            by.x = "sample.name", by.y = "Sample")
  outlier_summary <- table(outlier_by_plate$Plate)
  print(outlier_summary)
  
  cat("\nOutlier breakdown by issue:\n")
  print(table(multi_outliers$issue))
}

# ============================================================
# STEP 3: EXTRACT QC METRICS FOR ALL SAMPLES
# ============================================================
cat("\n--- STEP 3: Extracting QC metrics ---\n")

# Initialize metrics dataframe
qc_metrics_all <- data.frame(
  Sample = character(),
  Original_Sample = character(),
  Plate = character(),
  CallRate = numeric(),
  PredictedSex = character(),
  MedianMethylated = numeric(),
  MedianUnmethylated = numeric(),
  IsOutlier = logical(),
  FailureReason = character(),
  stringsAsFactors = FALSE
)

for (sample_name in names(all_qc_objects)) {
  qc_obj <- all_qc_objects[[sample_name]]
  
  # Get plate info
  plate_info <- plate_mapping[plate_mapping$Sample == sample_name, ]
  
  # Extract call rate
  call_rate <- NA
  if (!is.null(qc_obj$bad.probes.detectionp)) {
    if (length(qc_obj$bad.probes.detectionp) > 1) {
      call_rate <- mean(1 - qc_obj$bad.probes.detectionp, na.rm = TRUE)
    } else {
      call_rate <- 1 - qc_obj$bad.probes.detectionp
    }
  }
  
  # Extract predicted sex
  pred_sex <- if (!is.null(qc_obj$predicted.sex)) {
    as.character(qc_obj$predicted.sex)
  } else {
    NA
  }
  
  # Extract signal intensities
  median_meth <- if (!is.null(qc_obj$median.m.signal)) {
    qc_obj$median.m.signal
  } else {
    NA
  }
  
  median_unmeth <- if (!is.null(qc_obj$median.u.signal)) {
    qc_obj$median.u.signal
  } else {
    NA
  }
  
  # Check outlier status
  is_outlier <- sample_name %in% multi_outliers$sample.name
  
  failure_reason <- if (is_outlier) {
    idx <- which(multi_outliers$sample.name == sample_name)
    paste(multi_outliers$issue[idx], collapse = "; ")
  } else {
    "PASS"
  }
  
  # Add to dataframe
  qc_metrics_all <- rbind(qc_metrics_all, data.frame(
    Sample = sample_name,
    Original_Sample = plate_info$Original_Sample,
    Plate = plate_info$Plate,
    CallRate = call_rate,
    PredictedSex = pred_sex,
    MedianMethylated = median_meth,
    MedianUnmethylated = median_unmeth,
    IsOutlier = is_outlier,
    FailureReason = failure_reason,
    stringsAsFactors = FALSE
  ))
}

qc_metrics_all$QC_Status <- ifelse(qc_metrics_all$IsOutlier, "FAIL", "PASS")

cat("Extracted metrics for", nrow(qc_metrics_all), "samples\n")

# ============================================================
# STEP 4: CREATE COMPARISON VISUALIZATIONS
# ============================================================
cat("\n--- STEP 4: Creating comparison visualizations ---\n")

# Plot 1: Call rate distribution by plate
p1 <- ggplot(qc_metrics_all, aes(x = Plate, y = CallRate, fill = QC_Status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = qc_threshold, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("PASS" = "#4CAF50", "FAIL" = "#F44336")) +
  labs(title = "Call Rate Distribution Across All Plates",
       subtitle = paste("Red line: QC threshold (", qc_threshold, ")", sep = ""),
       y = "Call Rate",
       x = "Plate")

# Plot 2: Sex prediction distribution
sex_counts <- table(qc_metrics_all$Plate, qc_metrics_all$PredictedSex)
sex_df <- as.data.frame(sex_counts)
colnames(sex_df) <- c("Plate", "Sex", "Count")

p2 <- ggplot(sex_df, aes(x = Plate, y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predicted Sex Distribution by Plate",
       y = "Number of Samples",
       x = "Plate")

# Plot 3: Signal intensity comparison
p3 <- ggplot(qc_metrics_all, aes(x = MedianMethylated, y = MedianUnmethylated, 
                                  color = Plate, shape = QC_Status)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal(base_size = 12) +
  scale_shape_manual(values = c("PASS" = 16, "FAIL" = 4)) +
  labs(title = "Signal Intensities: Methylated vs Unmethylated",
       subtitle = "Shapes: circle = PASS, X = FAIL",
       x = "Median Methylated Signal",
       y = "Median Unmethylated Signal")

# Plot 4: QC status by plate
status_counts <- table(qc_metrics_all$Plate, qc_metrics_all$QC_Status)
status_df <- as.data.frame(status_counts)
colnames(status_df) <- c("Plate", "Status", "Count")

# Calculate percentages
status_df <- do.call(rbind, lapply(split(status_df, status_df$Plate), function(df) {
  df$Percentage <- 100 * df$Count / sum(df$Count)
  df
}))

p4 <- ggplot(status_df, aes(x = Plate, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("PASS" = "#4CAF50", "FAIL" = "#F44336")) +
  labs(title = "QC Pass/Fail Rates by Plate",
       y = "Number of Samples",
       x = "Plate")

# Save plots
cat("Saving visualizations...\n")

png(file.path(out_dir, "multi_plate_call_rates.png"),
    width = 1600, height = 1000, res = 150)
print(p1)
dev.off()

png(file.path(out_dir, "multi_plate_sex_distribution.png"),
    width = 1600, height = 1000, res = 150)
print(p2)
dev.off()

png(file.path(out_dir, "multi_plate_signal_intensities.png"),
    width = 1600, height = 1000, res = 150)
print(p3)
dev.off()

png(file.path(out_dir, "multi_plate_qc_status.png"),
    width = 1600, height = 1000, res = 150)
print(p4)
dev.off()

# Combined summary plot
png(file.path(out_dir, "multi_plate_summary.png"),
    width = 2400, height = 1600, res = 150)
grid.arrange(p1, p4, p2, p3, ncol = 2)
dev.off()

cat("✓ Visualizations saved\n")

# ============================================================
# STEP 5: COMPARE SINGLE-PLATE VS MULTI-PLATE OUTLIERS
# ============================================================
cat("\n--- STEP 5: Comparing single-plate vs multi-plate outliers ---\n")

# Load individual plate outliers
single_plate_outliers <- character()

for (qc_file in qc_object_files) {
  plate_id <- sub("_qc_objects\\.rds$", "", basename(qc_file))
  outliers_file <- file.path(dirname(qc_file), paste0(plate_id, "_outliers.rds"))
  
  if (file.exists(outliers_file)) {
    plate_outliers <- readRDS(outliers_file)
    if (nrow(plate_outliers) > 0) {
      # Add plate prefix to match our naming
      prefixed_names <- paste0(plate_id, "_", plate_outliers$sample.name)
      single_plate_outliers <- c(single_plate_outliers, prefixed_names)
    }
  }
}

multi_plate_outliers <- multi_outliers$sample.name

# Find differences
only_single_plate <- setdiff(single_plate_outliers, multi_plate_outliers)
only_multi_plate <- setdiff(multi_plate_outliers, single_plate_outliers)
both <- intersect(single_plate_outliers, multi_plate_outliers)

cat("\nOutlier comparison:\n")
cat("  Flagged in single-plate QC only:", length(only_single_plate), "\n")
cat("  Flagged in multi-plate QC only:", length(only_multi_plate), "\n")
cat("  Flagged in both:", length(both), "\n")

# Create comparison dataframe
comparison_df <- data.frame(
  Sample = qc_metrics_all$Sample,
  Original_Sample = qc_metrics_all$Original_Sample,
  Plate = qc_metrics_all$Plate,
  SinglePlate_Outlier = qc_metrics_all$Sample %in% single_plate_outliers,
  MultiPlate_Outlier = qc_metrics_all$Sample %in% multi_plate_outliers,
  CallRate = qc_metrics_all$CallRate,
  stringsAsFactors = FALSE
)

comparison_df$Outlier_Category <- "PASS"
comparison_df$Outlier_Category[comparison_df$SinglePlate_Outlier & 
                                !comparison_df$MultiPlate_Outlier] <- "Single-plate only"
comparison_df$Outlier_Category[!comparison_df$SinglePlate_Outlier & 
                                comparison_df$MultiPlate_Outlier] <- "Multi-plate only"
comparison_df$Outlier_Category[comparison_df$SinglePlate_Outlier & 
                                comparison_df$MultiPlate_Outlier] <- "Both"

write.csv(comparison_df, 
          file.path(out_dir, "outlier_comparison.csv"),
          row.names = FALSE)

cat("✓ Outlier comparison saved\n")

# ============================================================
# STEP 6: SAVE RESULTS
# ============================================================
cat("\n--- STEP 6: Saving results ---\n")

# Save multi-plate QC summary
saveRDS(multi_qc_summary, 
        file.path(out_dir, "multi_plate_qc_summary.rds"))
cat("✓ Multi-plate QC summary saved\n")

# Save outliers
saveRDS(multi_outliers,
        file.path(out_dir, "multi_plate_outliers.rds"))
cat("✓ Multi-plate outliers saved\n")

# Save all metrics
write.csv(qc_metrics_all,
          file.path(out_dir, "multi_plate_qc_metrics.csv"),
          row.names = FALSE)
cat("✓ Multi-plate QC metrics saved\n")

# Save list of samples passing multi-plate QC
passed_samples_multi <- qc_metrics_all$Sample[!qc_metrics_all$IsOutlier]
writeLines(passed_samples_multi,
           file.path(out_dir, "multi_plate_passed_samples.txt"))
cat("✓ Multi-plate passed samples list saved\n")

# ============================================================
# STEP 7: GENERATE HTML REPORT
# ============================================================
cat("\n--- STEP 7: Generating HTML QC report ---\n")

tryCatch({
  meffil.qc.report(
    multi_qc_summary,
    output.file = file.path(out_dir, "multi_plate_qc_report.html"),
    author = "Meffil Multi-Plate Pipeline",
    study = paste("Multi-Plate QC:", length(unique(plate_mapping$Plate)), "plates")
  )
  cat("✓ HTML QC report generated\n")
}, error = function(e) {
  cat("Warning: Could not generate HTML report:", conditionMessage(e), "\n")
})

# ============================================================
# STEP 8: CREATE TEXT SUMMARY REPORT
# ============================================================
cat("\n--- STEP 8: Creating text summary report ---\n")

summary_report <- c(
  paste(rep("=", 70), collapse = ""),
  "MULTI-PLATE QC REPORT",
  paste(rep("=", 70), collapse = ""),
  "",
  "OVERVIEW",
  paste0("  Total plates analyzed: ", length(unique(plate_mapping$Plate))),
  paste0("  Total samples: ", length(all_samples)),
  paste0("  QC threshold: ", qc_threshold),
  "",
  "RESULTS",
  paste0("  Samples passing QC: ", sum(!qc_metrics_all$IsOutlier), 
         " (", round(100 * sum(!qc_metrics_all$IsOutlier) / nrow(qc_metrics_all), 2), "%)"),
  paste0("  Samples failing QC: ", sum(qc_metrics_all$IsOutlier),
         " (", round(100 * sum(qc_metrics_all$IsOutlier) / nrow(qc_metrics_all), 2), "%)"),
  "",
  "PLATE BREAKDOWN"
)

# Add per-plate statistics
for (plate in unique(qc_metrics_all$Plate)) {
  plate_data <- qc_metrics_all[qc_metrics_all$Plate == plate, ]
  n_pass <- sum(!plate_data$IsOutlier)
  n_fail <- sum(plate_data$IsOutlier)
  pass_pct <- round(100 * n_pass / nrow(plate_data), 1)
  
  summary_report <- c(summary_report,
    paste0("  ", plate, ":"),
    paste0("    Total: ", nrow(plate_data)),
    paste0("    Pass: ", n_pass, " (", pass_pct, "%)"),
    paste0("    Fail: ", n_fail, " (", 100 - pass_pct, "%)")
  )
}

summary_report <- c(summary_report,
  "",
  "OUTLIER REASONS (Top 5)"
)

if (nrow(multi_outliers) > 0) {
  issue_table <- sort(table(multi_outliers$issue), decreasing = TRUE)
  top_issues <- head(issue_table, 5)
  for (i in seq_along(top_issues)) {
    summary_report <- c(summary_report,
      paste0("  ", i, ". ", names(top_issues)[i], ": ", top_issues[i], " samples")
    )
  }
} else {
  summary_report <- c(summary_report, "  No outliers detected")
}

summary_report <- c(summary_report,
  "",
  "SINGLE-PLATE VS MULTI-PLATE COMPARISON",
  paste0("  Outliers only in single-plate QC: ", length(only_single_plate)),
  paste0("  Outliers only in multi-plate QC: ", length(only_multi_plate)),
  paste0("  Outliers in both: ", length(both)),
  "",
  "INTERPRETATION",
  "  - Samples flagged only in single-plate QC may be outliers relative",
  "    to their plate but acceptable in the full dataset context",
  "  - Samples flagged only in multi-plate QC are true outliers that",
  "    may have been masked by plate-level variation",
  "  - Use multi-plate results for final sample exclusion decisions",
  "",
  paste(rep("=", 70), collapse = ""),
  "END OF REPORT",
  paste(rep("=", 70), collapse = "")
)

writeLines(summary_report, file.path(out_dir, "multi_plate_summary_report.txt"))
cat("✓ Text summary report saved\n")

# ============================================================
# FINAL SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("MULTI-PLATE QC REPORT COMPLETED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nSummary:\n")
cat("  Total plates:", length(unique(plate_mapping$Plate)), "\n")
cat("  Total samples:", length(all_samples), "\n")
cat("  Samples passing multi-plate QC:", sum(!qc_metrics_all$IsOutlier), "\n")
cat("  Samples failing multi-plate QC:", sum(qc_metrics_all$IsOutlier), "\n")
cat("  Overall pass rate:", 
    round(100 * sum(!qc_metrics_all$IsOutlier) / length(all_samples), 2), "%\n")

cat("\nOutput files in:", out_dir, "\n")
cat("  - multi_plate_qc_summary.rds\n")
cat("  - multi_plate_outliers.rds\n")
cat("  - multi_plate_qc_metrics.csv\n")
cat("  - multi_plate_passed_samples.txt\n")
cat("  - multi_plate_qc_report.html\n")
cat("  - multi_plate_summary_report.txt\n")
cat("  - outlier_comparison.csv\n")
cat("  - Visualization PNGs (5 files)\n")

cat("\n", rep("=", 60), "\n", sep = "")
