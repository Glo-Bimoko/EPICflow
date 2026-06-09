#!/usr/bin/env Rscript
# merge_plates.R
# Merge normalized data from all plates into final dataset
# Optionally apply cross-plate batch correction

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript merge_plates.R <normalized_dir> <out_dir>")
}

normalized_dir <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("MERGING NORMALIZED PLATES\n")
cat(rep("=", 60), "\n\n")

suppressPackageStartupMessages({
  library(sva)  # For ComBat if needed
  library(ggplot2)
})

# ============================================================
# STEP 1: FIND AND READ PLATE DATA
# ============================================================
cat("--- STEP 1: Loading plate data ---\n")

# Find all beta files
beta_files <- list.files(normalized_dir, pattern = "_beta\\.csv$", 
                        full.names = TRUE, recursive = FALSE)
M_files <- list.files(normalized_dir, pattern = "_M\\.csv$",
                     full.names = TRUE, recursive = FALSE)
info_files <- list.files(normalized_dir, pattern = "_sample_info\\.csv$",
                        full.names = TRUE, recursive = FALSE)

cat("Found", length(beta_files), "plate(s) with normalized data\n")

if (length(beta_files) == 0) {
  stop("No normalized plate data found in ", normalized_dir)
}

# Read all plates
beta_list <- list()
M_list <- list()
info_list <- list()

for (i in seq_along(beta_files)) {
  plate_name <- sub("_beta\\.csv$", "", basename(beta_files[i]))
  cat("  Loading", plate_name, "...\n")
  
  beta_list[[plate_name]] <- read.csv(beta_files[i], row.names = 1, check.names = FALSE)
  M_list[[plate_name]] <- read.csv(M_files[i], row.names = 1, check.names = FALSE)
  info_list[[plate_name]] <- read.csv(info_files[i], stringsAsFactors = FALSE)
}

cat("✓ All plates loaded\n")

# ============================================================
# STEP 2: MERGE PLATES
# ============================================================
cat("\n--- STEP 2: Merging plates ---\n")

# Get common probes across all plates
common_probes <- Reduce(intersect, lapply(beta_list, rownames))
cat("Common probes across all plates:", length(common_probes), "\n")

# Merge beta values
cat("Merging beta values...\n")
beta_merged <- do.call(cbind, lapply(beta_list, function(x) x[common_probes, ]))

# Merge M values
cat("Merging M values...\n")
M_merged <- do.call(cbind, lapply(M_list, function(x) x[common_probes, ]))

# Merge sample info
sample_info_merged <- do.call(rbind, info_list)

cat("✓ Merged data:\n")
cat("  Total samples:", ncol(beta_merged), "\n")
cat("  Total probes:", nrow(beta_merged), "\n")
cat("  Plates:", length(unique(sample_info_merged$Plate)), "\n")

# Save merged raw data (plate-normalized, not cross-plate corrected)
write.csv(beta_merged, 
          file.path(out_dir, "BetaValues_merged_plates.csv"),
          row.names = TRUE)
write.csv(M_merged,
          file.path(out_dir, "MValues_merged_plates.csv"),
          row.names = TRUE)
write.csv(sample_info_merged,
          file.path(out_dir, "Sample_info_all.csv"),
          row.names = FALSE)

# ============================================================
# STEP 3: ASSESS CROSS-PLATE BATCH EFFECTS
# ============================================================
cat("\n--- STEP 3: Assessing cross-plate batch effects ---\n")

# Only relevant if multiple plates
if (length(unique(sample_info_merged$Plate)) > 1) {
  
  cat("Analyzing variation between plates...\n")
  
  # PCA
  beta_complete <- beta_merged[complete.cases(beta_merged), ]
  probe_vars <- apply(beta_complete, 1, var)
  top_probes <- names(sort(probe_vars, decreasing = TRUE)[1:min(10000, length(probe_vars))])
  beta_for_pca <- beta_complete[top_probes, ]
  
  pca_result <- prcomp(t(beta_for_pca), scale. = TRUE)
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Plate = sample_info_merged$Plate,
    Sample = colnames(beta_for_pca),
    stringsAsFactors = FALSE
  )
  
  # Test plate effect
  plate_pval <- anova(lm(PC1 ~ Plate, data = pca_df))$'Pr(>F)'[1]
  cat("Plate effect test (PC1 ~ Plate): p =", format.pval(plate_pval, digits = 3), "\n")
  
  # Plot PCA colored by plate
  png(file.path(out_dir, "PCA_all_plates_before_combat.png"),
      width = 1400, height = 800, res = 150)
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate, shape = Plate)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(title = "PCA: All Plates (Before Cross-Plate Correction)",
         subtitle = paste("Plate effect p-value =", format.pval(plate_pval, digits = 3)),
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)"))
  print(p1)
  dev.off()
  
  # ============================================================
  # STEP 4: CROSS-PLATE BATCH CORRECTION (if needed)
  # ============================================================
  
  if (plate_pval < 0.05) {
    cat("\n--- STEP 4: Applying cross-plate batch correction ---\n")
    cat("⚠  Significant plate effects detected (p < 0.05)\n")
    cat("Applying ComBat for cross-plate correction...\n\n")
    
    # Create model matrix (intercept only - no biological covariates)
    # If you have phenotype data, add it here: mod <- model.matrix(~Disease, data=pheno)
    mod <- model.matrix(~1, data = data.frame(row.names = colnames(beta_merged)))
    
    # Apply ComBat
    cat("ComBat on beta values...\n")
    beta_combat <- ComBat(
      dat = beta_merged,
      batch = sample_info_merged$Plate,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
    
    cat("ComBat on M values...\n")
    M_combat <- ComBat(
      dat = M_merged,
      batch = sample_info_merged$Plate,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
    
    cat("✓ ComBat completed\n")
    
    # Post-ComBat PCA
    beta_complete_combat <- beta_combat[complete.cases(beta_combat), ]
    beta_for_pca_combat <- beta_complete_combat[top_probes, ]
    
    pca_combat <- prcomp(t(beta_for_pca_combat), scale. = TRUE)
    var_explained_combat <- summary(pca_combat)$importance[2, 1:2] * 100
    
    pca_df_combat <- data.frame(
      PC1 = pca_combat$x[, 1],
      PC2 = pca_combat$x[, 2],
      Plate = sample_info_merged$Plate,
      Sample = colnames(beta_for_pca_combat),
      stringsAsFactors = FALSE
    )
    
    plate_pval_after <- anova(lm(PC1 ~ Plate, data = pca_df_combat))$'Pr(>F)'[1]
    
    cat("\nBatch correction effectiveness:\n")
    cat("  Before ComBat: p =", format.pval(plate_pval, digits = 3), "\n")
    cat("  After ComBat:  p =", format.pval(plate_pval_after, digits = 3), "\n")
    
    # Plot post-ComBat PCA
    png(file.path(out_dir, "PCA_all_plates_after_combat.png"),
        width = 1400, height = 800, res = 150)
    p2 <- ggplot(pca_df_combat, aes(x = PC1, y = PC2, color = Plate, shape = Plate)) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal(base_size = 12) +
      labs(title = "PCA: All Plates (After Cross-Plate Correction)",
           subtitle = paste("Plate effect p-value =", format.pval(plate_pval_after, digits = 3)),
           x = paste0("PC1 (", round(var_explained_combat[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained_combat[2], 1), "%)"))
    print(p2)
    dev.off()
    
    # Save ComBat-corrected data as FINAL
    write.csv(beta_combat,
              file.path(out_dir, "BetaValues_FINAL.csv"),
              row.names = TRUE)
    write.csv(M_combat,
              file.path(out_dir, "MValues_FINAL.csv"),
              row.names = TRUE)
    
    combat_applied <- TRUE
    
  } else {
    cat("\n--- STEP 4: Cross-plate correction ---\n")
    cat("✓ No significant plate effects (p ≥ 0.05)\n")
    cat("Plate-normalized data is sufficient - ComBat not needed\n")
    
    # Use merged data as final
    write.csv(beta_merged,
              file.path(out_dir, "BetaValues_FINAL.csv"),
              row.names = TRUE)
    write.csv(M_merged,
              file.path(out_dir, "MValues_FINAL.csv"),
              row.names = TRUE)
    
    combat_applied <- FALSE
  }
  
} else {
  cat("Only one plate - no cross-plate batch effects possible\n")
  
  # Single plate - just copy to FINAL
  write.csv(beta_merged,
            file.path(out_dir, "BetaValues_FINAL.csv"),
            row.names = TRUE)
  write.csv(M_merged,
            file.path(out_dir, "MValues_FINAL.csv"),
            row.names = TRUE)
  
  combat_applied <- FALSE
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("FINAL DATASET CREATED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nSummary:\n")
cat("  Total samples:", ncol(beta_merged), "\n")
cat("  Total probes:", nrow(beta_merged), "\n")
cat("  Plates merged:", length(unique(sample_info_merged$Plate)), "\n")

if (length(unique(sample_info_merged$Plate)) > 1) {
  cat("\nBatch correction:\n")
  cat("  Within-plate: Meffil functional normalization\n")
  if (combat_applied) {
    cat("  Cross-plate: ComBat (applied)\n")
  } else {
    cat("  Cross-plate: Not needed\n")
  }
}

cat("\nOutput files in:", out_dir, "\n")
cat("  - Sample_info_all.csv (all sample metadata)\n")
cat("  - BetaValues_merged_plates.csv (plate-normalized)\n")
cat("  - MValues_merged_plates.csv (plate-normalized)\n")
cat("  - BetaValues_FINAL.csv (FINAL - use this for analysis)\n")
cat("  - MValues_FINAL.csv (FINAL - use this for analysis)\n")
if (length(unique(sample_info_merged$Plate)) > 1) {
  cat("  - PCA plots (before/after cross-plate correction)\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Pipeline completed successfully!\n")
cat("FINAL data ready for downstream analysis.\n")
cat(rep("=", 60), "\n\n")
