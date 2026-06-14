#!/usr/bin/env Rscript
# bin/combined_normalize_meffil.r
# Combined functional normalization across all plates using meffil.
#
# References:
#   https://github.com/perishky/meffil/wiki/Sample-QC
#   https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets
#   Min et al. (2018) Bioinformatics 34:3983-3989
#
# Usage:
#   Rscript combined_normalize_meffil.r \
#       <qc_objects_dir> <passed_samples_dir> <qc_summaries_dir> \
#       <combined_samplesheet> <out_dir> <qc_threshold> [max_gb] [n_pcs]
#
# Arguments:
#   qc_objects_dir        : directory containing per-plate *_qc_objects.rds
#   passed_samples_dir    : directory containing per-plate *_passed_samples.txt
#   qc_summaries_dir      : directory containing per-plate *_qc_summary.rds
#                           (used to collect bad CpGs for cpglist.remove)
#   combined_samplesheet  : merged samplesheet CSV from MERGE_SAMPLESHEETS
#   out_dir               : output directory
#   qc_threshold          : call-rate threshold (passed through for reporting)
#   max_gb                : (optional) RAM ceiling in GB, default 32.
#                           Controls max.bytes in meffil.normalize.samples().
#   n_pcs                 : (optional) number of control-probe PCs for FN.
#                           If 0 or omitted, determined automatically via
#                           meffil.plot.pc.fit() cross-validated scree plot.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript combined_normalize_meffil.r",
    "<qc_objects_dir> <passed_samples_dir> <qc_summaries_dir>",
    "<combined_samplesheet> <out_dir> <qc_threshold> [max_gb] [n_pcs]"
  ))
}

qc_objects_dir          <- args[1]
passed_samples_dir      <- args[2]
qc_summaries_dir        <- args[3]
combined_samplesheet_file <- args[4]
out_dir                 <- args[5]
qc_threshold            <- as.numeric(args[6])
max_gb                  <- if (length(args) >= 7) as.numeric(args[7]) else 32
n_pcs_arg               <- if (length(args) >= 8) as.integer(args[8]) else 0L

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("COMBINED NORMALIZATION (All Plates Together)\n")
cat(rep("=", 60), "\n\n")
cat("QC threshold  :", qc_threshold, "\n")
cat("Memory ceiling:", max_gb, "GB\n")
cat("n_pcs arg     :", if (n_pcs_arg == 0L) "auto" else n_pcs_arg, "\n\n")

# ── Threading ──────────────────────────────────────────────────────────────
# meffil.normalize.samples() parallelises safely — each sample is normalised
# independently (paper Section 2.2).  Use all available cores.
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                      Sys.getenv("PBS_NUM_PPN", 4)))
cat("Cores         :", n_cores, "\n")
options(mc.cores = n_cores)

# max.bytes per fork — same formula as plate_qc.
#
# FIX: use as.numeric() (64-bit double) instead of as.integer() (32-bit signed).
# as.integer() overflows to NA for values > 2,147,483,647 (~2.1 GB), which
# happens whenever (max_gb * 0.75 * 1024^3) / n_cores exceeds that ceiling —
# e.g. max_gb=28, n_cores=4 yields ~5.6 GB, silently becoming NA.
# meffil then crashes with:
#   "Error in if (n.fun < 1): missing value where TRUE/FALSE needed"
max_bytes <- as.numeric((max_gb * 0.75 * 1024^3) / n_cores)

# Defensive guard: if max_bytes is somehow NA or implausibly small,
# fall back to a safe 1 GB per fork rather than crashing deep inside meffil.
if (is.na(max_bytes) || max_bytes < 1e6) {
  warning(sprintf(
    "max_bytes calculation produced an invalid value (%s). Falling back to 1 GB per fork. Check that max_gb ('%s') is a valid positive number.",
    max_bytes, args[7]
  ))
  max_bytes <- 1e9
}

cat("max.bytes/fork:", format(max_bytes, big.mark = ",", scientific = FALSE), "bytes\n\n")

suppressPackageStartupMessages({
  library(meffil)
  library(ggplot2)
})

# ============================================================
# STEP 1: LOAD QC OBJECTS AND FILTER TO PASSED SAMPLES
# ============================================================
cat("--- STEP 1: Loading QC objects ---\n")

qc_object_files <- list.files(qc_objects_dir,
                               pattern = "_qc_objects\\.rds$",
                               full.names = TRUE)
if (length(qc_object_files) == 0)
  stop("No *_qc_objects.rds files found in: ", qc_objects_dir)

cat("Plates found:", length(qc_object_files), "\n")

all_qc_objects <- list()
plate_of_sample <- character(0)   # named vector: sample -> plate_id

for (f in qc_object_files) {
  plate_id  <- sub("_qc_objects\\.rds$", "", basename(f))
  plate_qc  <- readRDS(f)
  cat("  ", plate_id, ":", length(plate_qc), "QC objects\n")
  for (s in names(plate_qc)) {
    all_qc_objects[[s]]  <- plate_qc[[s]]
    plate_of_sample[[s]] <- plate_id
  }
  rm(plate_qc); gc()
}
cat("Total QC objects loaded:", length(all_qc_objects), "\n")

# ── Passed samples ─────────────────────────────────────────────────────────
cat("\n--- STEP 1b: Loading passed-sample lists ---\n")

passed_files <- list.files(passed_samples_dir,
                            pattern = "_passed_samples\\.txt$",
                            full.names = TRUE)

if (length(passed_files) == 0)
  stop("No *_passed_samples.txt files found in: ", passed_samples_dir)

all_passed <- character(0)
for (f in passed_files) {
  plate_id <- sub("_passed_samples\\.txt$", "", basename(f))

  # FIX: readLines() on an empty file (written by plate_qc when meffil.qc()
  # crashed in a prior run) returns character(0), not an error. Catch this
  # explicitly and warn so the user knows which plate produced no passing
  # samples — it almost always means the plate_qc step needs to be re-run
  # (e.g. after fixing the max_bytes integer-overflow bug in plate_qc_meffil.r).
  samples <- tryCatch(
    readLines(f, warn = FALSE),
    error = function(e) {
      cat("  WARNING: could not read", basename(f), "—", conditionMessage(e), "\n")
      character(0)
    }
  )

  # Strip any blank lines that could arise from platform line-ending differences
  samples <- samples[nzchar(trimws(samples))]

  cat("  ", plate_id, ":", length(samples), "passed\n")

  if (length(samples) == 0) {
    cat("  WARNING: no passed samples for", plate_id,
        "— was plate_qc_meffil.r run successfully for this plate?\n")
  }

  all_passed <- c(all_passed, samples)
}
cat("Total passed samples:", length(all_passed), "\n")

qc_objects_passed <- all_qc_objects[names(all_qc_objects) %in% all_passed]
cat("QC objects for normalization:", length(qc_objects_passed), "\n")

if (length(qc_objects_passed) == 0) {
  stop(paste0(
    "No samples passed QC — cannot normalise.\n",
    "  Possible causes:\n",
    "  1. plate_qc_meffil.r crashed before writing *_passed_samples.txt\n",
    "     (check per-plate logs; look for the 'NAs introduced by coercion'\n",
    "      warning which indicates the max_bytes integer-overflow bug).\n",
    "  2. All samples genuinely failed the call-rate threshold (", qc_threshold, ").\n",
    "     Inspect the per-plate *_qc_metrics.csv files in qc_results/.\n",
    "  3. Sample names in *_passed_samples.txt do not match names(qc_objects).\n",
    "     This can happen if the samplesheet Sample_Name column was changed\n",
    "     between plate_qc and combined_normalize runs."
  ))
}

# ============================================================
# STEP 2: COLLECT BAD CpGS FROM PER-PLATE QC SUMMARIES
# ============================================================
# The meffil wiki and paper recommend passing cpglist.remove to
# meffil.normalize.samples() to exclude probes flagged as poor quality
# during QC (high fraction of undetected calls or low bead numbers).
# We aggregate bad.cpgs$name across all per-plate QC summaries.
# ============================================================
cat("\n--- STEP 2: Collecting bad CpGs from QC summaries ---\n")

summary_files <- list.files(qc_summaries_dir,
                             pattern = "_qc_summary\\.rds$",
                             full.names = TRUE)

bad_cpgs <- character(0)

if (length(summary_files) == 0) {
  cat("WARNING: No *_qc_summary.rds files found in:", qc_summaries_dir, "\n")
  cat("         Proceeding without cpglist.remove — suboptimal but not fatal.\n")
} else {
  for (f in summary_files) {
    plate_id <- sub("_qc_summary\\.rds$", "", basename(f))
    qs       <- readRDS(f)
    if (!is.null(qs$bad.cpgs) && nrow(qs$bad.cpgs) > 0) {
      n_bad <- nrow(qs$bad.cpgs)
      bad_cpgs <- union(bad_cpgs, qs$bad.cpgs$name)
      cat("  ", plate_id, ":", n_bad, "bad CpGs\n")
    } else {
      cat("  ", plate_id, ": no bad CpGs\n")
    }
    rm(qs); gc()
  }
  cat("Total unique bad CpGs to remove:", length(bad_cpgs), "\n")
}

# ============================================================
# STEP 3: DETERMINE NUMBER OF CONTROL-PROBE PCs
# ============================================================
# The paper (Section 3.2.2) shows that PC choice matters greatly for
# downstream EWAS sensitivity/specificity.  The recommended approach is
# meffil.plot.pc.fit() which uses 10-fold cross-validation to find the
# number that minimises residual probe-quantile variance.
# The paper found 10 PCs optimal for their dataset; we use that as the
# default if automatic selection fails or produces an implausible result.
# ============================================================
cat("\n--- STEP 3: Determining number of normalization PCs ---\n")

PC_DEFAULT <- 10L

if (n_pcs_arg > 0L) {
  n_pcs <- n_pcs_arg
  cat("Using user-supplied n_pcs:", n_pcs, "\n")
} else {
  cat("Running cross-validated PC selection (meffil.plot.pc.fit)...\n")
  n_pcs <- tryCatch({
    pc_fit <- meffil.plot.pc.fit(qc_objects_passed)

    # ── Robustly identify the n.pcs and MSR columns ────────────────────────
    # meffil returns a data frame but the exact column names differ across
    # versions.  Rather than hardcoding, we:
    #   1. Find the column whose name contains "pc" (case-insensitive) —
    #      that is the number-of-PCs axis.
    #   2. Average over cross-validation groups so each PC count has one
    #      mean residual, then pick the minimum.
    # This is robust to "n.pcs"/"number.pcs" and "mean.sq"/"mean.residuals"
    # naming differences across meffil versions.
    df <- pc_fit$data
    cat("pc_fit$data columns:", paste(colnames(df), collapse = ", "), "\n")

    # Column containing the number of PCs (integer-like, named with "pc")
    pc_col <- colnames(df)[grep("^n\\.pcs$|^number|^pcs$|pc", colnames(df),
                                ignore.case = TRUE)][1]
    if (is.na(pc_col)) pc_col <- colnames(df)[sapply(df, is.numeric)][1]

    # Column containing the residual / MSR (numeric, NOT the pc column,
    # NOT a group/fold column)
    numeric_cols <- colnames(df)[sapply(df, is.numeric)]
    group_cols   <- colnames(df)[grep("group|fold|cross|iter", colnames(df),
                                      ignore.case = TRUE)]
    msr_col <- setdiff(numeric_cols, c(pc_col, group_cols))[1]

    cat("Using PC column:", pc_col, "  MSR column:", msr_col, "\n")

    # Average MSR per PC count (handles per-fold rows)
    mean_msr <- tapply(df[[msr_col]], df[[pc_col]], mean, na.rm = TRUE)
    best_pc  <- as.integer(names(mean_msr)[which.min(mean_msr)])

    cat("Cross-validated scree minimum at:", best_pc, "PCs\n")

    # Save the scree plot
    png(file.path(out_dir, "pc_fit_scree.png"), width = 800, height = 500, res = 120)
    print(pc_fit$plot)
    dev.off()
    cat("✓ PC fit scree plot saved: pc_fit_scree.png\n")

    best_pc
  }, error = function(e) {
    cat("WARNING: meffil.plot.pc.fit failed (", conditionMessage(e), ")\n")
    cat("Falling back to default:", PC_DEFAULT, "PCs\n")
    PC_DEFAULT
  })
  # Sanity bounds: meffil supports max 42 control-probe PCs
  if (is.null(n_pcs) || length(n_pcs) == 0 || is.na(n_pcs) ||
      n_pcs < 1L || n_pcs > 42L) {
    cat("PC selection returned", n_pcs, "— using default:", PC_DEFAULT, "\n")
    n_pcs <- PC_DEFAULT
  }
}
cat("Final n_pcs for normalization:", n_pcs, "\n")


# ============================================================
# STEP 4: NORMALISE QUANTILES (control-probe PCA step)
# ============================================================
cat("\n--- STEP 4: Normalising quantiles ---\n")
cat("Samples:", length(qc_objects_passed), "\n")

# ── Pre-flight: check control matrix variance ──────────────────────────────
# If all control probes are zero-variance, meffil.normalize.quantiles() will
# crash inside meffil.control.matrix(). Diagnose early with a clear message.
cat("Pre-flight: checking control matrix variance...\n")
ctrl_check <- tryCatch({
  cm <- meffil.control.matrix(qc_objects_passed, normalize = FALSE)
  col_vars <- apply(cm, 2, var, na.rm = TRUE)
  n_zero   <- sum(col_vars == 0, na.rm = TRUE)
  cat("  Control matrix:", nrow(cm), "samples x", ncol(cm), "probes\n")
  cat("  Zero-variance probes:", n_zero, "/", ncol(cm), "\n")
  if (n_zero == ncol(cm)) {
    cat("  CRITICAL: All control matrix probes have zero variance!\n")
    cat("  This almost always means the QC objects were built with a\n")
    cat("  mismatched featureset (e.g. 'epic' instead of 'epicv2').\n")
    cat("  Check: sapply(qc_objects_passed[1:3], function(x) x$featureset)\n")
    # Print featuresets seen
    fs <- unique(sapply(qc_objects_passed, function(x)
                        tryCatch(x$featureset, error = function(e) "UNKNOWN")))
    cat("  Featuresets in QC objects:", paste(fs, collapse = ", "), "\n")
  }
  list(ok = n_zero < ncol(cm), n_zero = n_zero, total = ncol(cm))
}, error = function(e) {
  cat("  WARNING: control matrix pre-check failed:", conditionMessage(e), "\n")
  list(ok = TRUE, n_zero = NA, total = NA)   # proceed and let meffil report it
})

norm_objects <- meffil.normalize.quantiles(
  qc.objects     = qc_objects_passed,
  number.pcs     = n_pcs,
  fixed.effects  = NULL,   # ← ADD THIS: prevents implicit batch variable injection
  random.effects = NULL,   # ← ADD THIS
  verbose        = TRUE
)
cat("✓ Quantile normalisation complete\n")

# ============================================================
# STEP 5: NORMALISE SAMPLES → GDS FILE
# ============================================================
# Writing to a GDS file (Genomic Data Structure) instead of keeping
# the full beta matrix in R memory.  For 935k probes × 288 samples,
# the in-memory matrix peaks at ~67 GB (paper Table 1); the GDS path
# peaks at ~3 GB.  Downstream operations use meffil.gds.methylation()
# or meffil.gds.apply() to access data one CpG at a time.
# ============================================================
cat("\n--- STEP 5: Normalising samples → GDS file ---\n")

gds_file <- file.path(out_dir, "beta_normalised.gds")
cat("GDS output path:", gds_file, "\n")
cat("remove.poor.signal: TRUE (sets NA where detection p-val too high or bead count too low)\n")
cat("cpglist.remove:", length(bad_cpgs), "bad CpGs\n\n")

meffil.normalize.samples(
  norm.objects       = norm_objects,
  just.beta          = TRUE,
  remove.poor.signal = TRUE,
  cpglist.remove     = if (length(bad_cpgs) > 0) bad_cpgs else NULL,
  gds.filename       = gds_file,
  max.bytes          = max_bytes,
  verbose            = TRUE
)
cat("✓ GDS file written:", gds_file, "\n")

# ============================================================
# STEP 6: PCA ON AUTOSOMAL SITES FOR NORMALISATION REPORT
# ============================================================
# The wiki shows this exact sequence for generating the normalisation
# report for large datasets:
#   autosomal.sites <- meffil.get.autosomal.sites("epicv2")
#   pcs             <- meffil.methylation.pcs(gds.filename, sites=autosomal.sites)
#   norm.summary    <- meffil.normalization.summary(norm.objects, pcs=pcs)
# ============================================================
cat("\n--- STEP 6: Computing methylation PCs for normalisation report ---\n")

# Determine featureset from first QC object
featureset <- tryCatch(
  qc_objects_passed[[1]]$featureset,
  error = function(e) "epic"
)
cat("Array featureset:", featureset, "\n")

autosomal_sites <- tryCatch(
  meffil.get.autosomal.sites(featureset),
  error = function(e) {
    cat("WARNING: could not retrieve autosomal sites for", featureset,
        "— falling back to 'epic'\n")
    meffil.get.autosomal.sites("epic")
  }
)
cat("Autosomal CpG sites:", length(autosomal_sites), "\n")

pcs_for_report <- tryCatch({
  cat("Computing methylation PCs from GDS (this reads the file row-wise)...\n")
  meffil.methylation.pcs(gds_file, sites = autosomal_sites, verbose = TRUE)
}, error = function(e) {
  cat("WARNING: meffil.methylation.pcs failed:", conditionMessage(e), "\n")
  NULL
})

# ── Normalisation report ───────────────────────────────────────────────────
cat("\n--- STEP 6b: Generating normalisation report ---\n")

tryCatch({
  norm_summary <- meffil.normalization.summary(
    norm.objects = norm_objects,
    pcs          = pcs_for_report
  )
  saveRDS(norm_summary, file.path(out_dir, "combined_norm_summary.rds"))

  norm_report <- file.path(out_dir, "combined_normalization_report.html")
  meffil.normalization.report(
    norm_summary,
    output.file = norm_report,
    author      = "Meffil Pipeline",
    study       = paste0("Combined normalisation — ",
                         length(unique(plate_of_sample[all_passed])),
                         " plates, ", length(all_passed), " samples")
  )
  cat("✓ Normalisation report:", basename(norm_report), "\n")
}, error = function(e) {
  cat("WARNING: normalisation report failed:", conditionMessage(e), "\n")
})

# ============================================================
# STEP 7: SAMPLE METADATA AND PLATE SUMMARY
# ============================================================
cat("\n--- STEP 7: Saving sample metadata ---\n")

# Build sample info table for passed samples
sample_info <- data.frame(
  Sample_Name = all_passed,
  Plate_ID    = plate_of_sample[all_passed],
  stringsAsFactors = FALSE
)

# Optionally join combined samplesheet columns
if (file.exists(combined_samplesheet_file)) {
  cs <- read.csv(combined_samplesheet_file, stringsAsFactors = FALSE)
  cs_cols <- setdiff(colnames(cs), "Sample_Name")
  sample_info <- merge(sample_info, cs[, c("Sample_Name", cs_cols), drop = FALSE],
                       by = "Sample_Name", all.x = TRUE)
}

write.csv(sample_info,
          file.path(out_dir, "sample_info_combined.csv"),
          row.names = FALSE)
cat("✓ sample_info_combined.csv written (", nrow(sample_info), "samples )\n")

# Plate summary counts
plate_counts <- table(sample_info$Plate_ID)
cat("\nSamples per plate after QC:\n")
for (pl in names(sort(plate_counts, decreasing = TRUE)))
  cat(sprintf("  %-25s %d\n", pl, plate_counts[pl]))

# ============================================================
# STEP 8: PCA VISUALISATION OF PLATE EFFECTS
# ============================================================
# Read a top-variance subset from the GDS file to keep this in memory.
cat("\n--- STEP 8: PCA visualisation of plate effects ---\n")

tryCatch({
  # Load all autosomal sites for variance calculation — use apply via GDS
  cat("Computing per-probe variance from GDS...\n")
  probe_vars <- meffil.gds.apply(
    gds_file, bysite = TRUE, type = "double",
    FUN = var, na.rm = TRUE
  )
  top_sites <- names(sort(probe_vars, decreasing = TRUE))[1:min(10000, length(probe_vars))]

  cat("Loading top", length(top_sites), "variable probes for PCA...\n")
  beta_subset <- meffil.gds.methylation(gds_file, sites = top_sites)

  # Remove rows with any NA (prcomp doesn't handle NAs)
  beta_subset <- beta_subset[complete.cases(beta_subset), ]
  cat("Complete-case probes for PCA:", nrow(beta_subset), "\n")

  pca_result  <- prcomp(t(beta_subset), scale. = TRUE)
  var_exp     <- summary(pca_result)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1         = pca_result$x[, 1],
    PC2         = pca_result$x[, 2],
    Sample_Name = rownames(pca_result$x),
    stringsAsFactors = FALSE
  )
  pca_df <- merge(pca_df, sample_info[, c("Sample_Name", "Plate_ID")],
                  by = "Sample_Name")

  plate_pval <- NA
  if (length(unique(pca_df$Plate_ID)) > 1) {
    plate_pval <- anova(lm(PC1 ~ Plate_ID, data = pca_df))$"Pr(>F)"[1]
    cat("Plate effect on PC1 (ANOVA p):", format.pval(plate_pval, digits = 3), "\n")
  }

  png(file.path(out_dir, "PCA_combined_normalization.png"),
      width = 1600, height = 1000, res = 150)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Plate_ID, shape = Plate_ID)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title    = "PCA: Combined Normalisation (All Plates)",
      subtitle = if (!is.na(plate_pval))
                   paste("Plate effect ANOVA p =", format.pval(plate_pval, digits = 3))
                 else "Single plate",
      x        = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y        = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )
  print(p)
  dev.off()
  cat("✓ PCA plot saved: PCA_combined_normalization.png\n")

}, error = function(e) {
  cat("WARNING: PCA visualisation failed:", conditionMessage(e), "\n")
})

# Save norm objects for reproducibility / downstream EWAS
saveRDS(norm_objects, file.path(out_dir, "combined_norm_objects.rds"))
cat("✓ Normalisation objects saved: combined_norm_objects.rds\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("COMBINED NORMALIZATION COMPLETED\n")
cat(rep("=", 60), "\n")
cat("  Samples normalised :", length(all_passed), "\n")
cat("  Plates             :", length(unique(plate_of_sample[all_passed])), "\n")
cat("  Control-probe PCs  :", n_pcs, "\n")
cat("  Bad CpGs removed   :", length(bad_cpgs), "\n")
cat("  Output directory   :", out_dir, "\n\n")
cat("Key outputs:\n")
cat("  beta_normalised.gds               — normalised beta values (GDS format)\n")
cat("  combined_normalization_report.html — meffil normalisation QC report\n")
cat("  combined_norm_objects.rds         — norm objects (for EWAS / reproducibility)\n")
cat("  sample_info_combined.csv          — sample metadata\n")
cat("  PCA_combined_normalization.png    — plate-effect visualisation\n")
cat("  pc_fit_scree.png                  — cross-validated PC selection plot\n")
cat("\nTo access beta values downstream:\n")
cat("  library(meffil)\n")
cat("  beta <- meffil.gds.methylation('beta_normalised.gds',\n")
cat("              sites=my_sites, samples=my_samples)\n")
cat("\n", rep("=", 60), "\n", sep = "")
