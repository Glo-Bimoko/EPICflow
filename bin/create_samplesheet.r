#!/usr/bin/env Rscript
# create_samplesheet.r
# Creates meffil samplesheet from all IDAT files in a plate folder.
# When a sample map CSV is supplied via --sample_map, the Sample_Name column
# is populated with study-level Sample IDs (e.g. G001) instead of raw
# BeadChip barcodes, making all downstream outputs participant-centric.
#
# Usage (without map):
#   Rscript create_samplesheet.r <plate_manifest.txt> <out_samplesheet.csv>
#
# Usage (with sample map):
#   Rscript create_samplesheet.r <plate_manifest.txt> <out_samplesheet.csv> \
#       --sample_map <all_samples.csv>
#
# Sample map CSV must contain (case-insensitive column names):
#   Sample ID | BeadChip Barcode | Sentrix Position
# The pipeline --sample_map parameter points to this file (e.g. all_6_epic.csv).

suppressPackageStartupMessages(library(methods))

# ── Argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_samplesheet.r <plate_manifest.txt> <out_samplesheet.csv> [--sample_map <map.csv>]")
}

plate_manifest  <- args[1]
out_samplesheet <- args[2]

# Optional --sample_map flag
sample_map_file <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--sample_map" && i < length(args)) {
    sample_map_file <- args[i + 1]
  }
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("CREATING MEFFIL SAMPLESHEET FOR PLATE\n")
cat(rep("=", 60), "\n\n")

plate_id <- sub("\\.txt$", "", basename(plate_manifest))
cat("Plate:", plate_id, "\n")

# ── Read IDAT file paths from manifest ───────────────────────────────────────
idat_files <- readLines(plate_manifest)
if (length(idat_files) == 0) {
  stop("Empty plate manifest: ", plate_manifest)
}
cat("Found", length(idat_files), "IDAT files in manifest\n")

# Get unique basenames (strip _Grn.idat / _Red.idat)
basenames <- unique(sub("_(Grn|Red)\\.idat$", "", idat_files, ignore.case = TRUE))
cat("Unique chip positions:", length(basenames), "\n")

# Extract barcode and position from each basename
# Expected filename pattern: <BeadChipBarcode>_<SentrixPosition>
bn_base     <- basename(basenames)             # e.g. "208789590164_R01C01"
barcode_vec <- sub("_([^_]+)$", "", bn_base)  # "208789590164"
pos_vec     <- sub(".*_",        "", bn_base)  # "R01C01"

cat("\nFirst few chip positions:\n")
for (i in seq_len(min(3, length(bn_base)))) {
  cat("  ", bn_base[i], "  barcode:", barcode_vec[i], "  position:", pos_vec[i], "\n")
}

# ── Load sample map and resolve Sample IDs ────────────────────────────────────
if (!is.null(sample_map_file)) {
  cat("\nLoading sample map:", sample_map_file, "\n")
  if (!file.exists(sample_map_file)) {
    stop("Sample map file not found: ", sample_map_file)
  }

  smap <- read.csv(sample_map_file, stringsAsFactors = FALSE, check.names = FALSE)

  # Normalise column names: lowercase + replace spaces/dots with underscores
  colnames(smap) <- tolower(gsub("[ .]", "_", colnames(smap)))

  required_cols <- c("sample_id", "beadchip_barcode", "sentrix_position")
  missing_cols  <- setdiff(required_cols, colnames(smap))
  if (length(missing_cols) > 0) {
    stop("Sample map is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         "\nFound columns: ", paste(colnames(smap), collapse = ", "))
  }

  # Normalise key columns for matching
  smap$beadchip_barcode    <- as.character(trimws(smap$beadchip_barcode))
  smap$sentrix_position    <- toupper(trimws(smap$sentrix_position))
  smap$sample_id           <- trimws(smap$sample_id)

  # Collect sex if available (Collected_Gender / Sex / Gender)
  sex_col <- intersect(c("collected_gender", "sex", "gender"), colnames(smap))[1]
  has_sex  <- !is.na(sex_col)

  # Build lookup key: barcode_position (same convention as IDAT basename)
  smap$lookup_key <- paste0(smap$beadchip_barcode, "_", smap$sentrix_position)

  cat("Sample map loaded:", nrow(smap), "rows\n")

  # Filter to rows matching this plate (optional guard — map may cover all plates)
  # We rely on the key match below; no explicit plate filter needed.

  # Match each IDAT position to a Sample ID
  idat_keys <- paste0(barcode_vec, "_", toupper(pos_vec))
  matched   <- smap[match(idat_keys, smap$lookup_key), ]

  n_matched   <- sum(!is.na(matched$sample_id))
  n_unmatched <- sum(is.na(matched$sample_id))

  cat("Matched:", n_matched, "/ Unmatched:", n_unmatched, "\n")

  if (n_unmatched > 0) {
    unmatched_keys <- idat_keys[is.na(matched$sample_id)]
    cat("WARNING – the following IDAT positions have no entry in the sample map:\n")
    for (k in unmatched_keys) cat("  ", k, "\n")
    cat("These positions will use the raw barcode as Sample_Name.\n")
  }

  # Resolve Sample_Name: use mapped ID where available, fall back to barcode
  sample_names <- ifelse(!is.na(matched$sample_id),
                         matched$sample_id,
                         bn_base)

  # Resolve Sex: use map value where available
  if (has_sex) {
    sex_values <- ifelse(!is.na(matched$sample_id),
                         as.character(matched[[sex_col]]),
                         "NA")
    # Standardise to M/F/NA — meffil expects M, F, or NA
    sex_values <- ifelse(tolower(sex_values) %in% c("male",   "m"), "M",
                  ifelse(tolower(sex_values) %in% c("female", "f"), "F", "NA"))
  } else {
    sex_values <- rep("NA", length(basenames))
  }

  # Plate info from the map (use map value for matched rows; fall back to plate_id)
  if ("plate_number" %in% colnames(smap)) {
    plate_values <- ifelse(!is.na(matched$sample_id),
                           matched$plate_number,
                           plate_id)
  } else {
    plate_values <- rep(plate_id, length(basenames))
  }

} else {
  # No sample map — use barcode basename as before
  cat("\nNo --sample_map provided; using raw IDAT basenames as Sample_Name.\n")
  sample_names <- bn_base
  sex_values   <- rep("NA", length(basenames))
  plate_values <- rep(plate_id, length(basenames))
}

# ── Build samplesheet ─────────────────────────────────────────────────────────
samplesheet <- data.frame(
  Sample_Name = sample_names,
  Sex         = sex_values,
  Basename    = basenames,
  Plate_ID    = plate_values,
  stringsAsFactors = FALSE
)

# Ensure unique Sample_Name values (shouldn't happen with a good map,
# but guards against duplicates in the raw barcode fallback)
if (any(duplicated(samplesheet$Sample_Name))) {
  cat("WARNING: Duplicate Sample_Name values detected; making them unique.\n")
  samplesheet$Sample_Name <- make.unique(samplesheet$Sample_Name, sep = "_dup")
}

cat("\nCreated samplesheet with", nrow(samplesheet), "samples\n")

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(dirname(out_samplesheet), showWarnings = FALSE, recursive = TRUE)
write.csv(samplesheet, out_samplesheet, row.names = FALSE, quote = TRUE)

cat("\nSamplesheet saved to:", out_samplesheet, "\n")
cat("\nSamplesheet preview:\n")
print(head(samplesheet[, c("Sample_Name", "Sex", "Plate_ID")], 6))

cat("\n", rep("=", 60), "\n", sep = "")
cat("SAMPLESHEET CREATION COMPLETED\n")
cat(rep("=", 60), "\n\n")

cat("Summary:\n")
cat("  Samples             :", nrow(samplesheet), "\n")
cat("  Mapped from file    :", if (!is.null(sample_map_file)) sum(samplesheet$Sample_Name != bn_base) else 0, "\n")
cat("  Fallback (barcode)  :", if (!is.null(sample_map_file)) sum(samplesheet$Sample_Name == bn_base) else nrow(samplesheet), "\n")
cat("  Sex values present  :", sum(samplesheet$Sex %in% c("M", "F")), "\n")
