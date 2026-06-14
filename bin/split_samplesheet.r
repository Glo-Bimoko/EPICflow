#!/usr/bin/env Rscript
# bin/split_samplesheet.r
#
# Reads the study-level samplesheet CSV, splits by Plate Number, and writes
# one meffil-ready CSV per plate.
#
# KEY DESIGN: directory-agnostic IDAT discovery
#   Plate folder names on disk almost never match the Plate Number values in
#   the samplesheet (e.g. "Plate_01" vs "Plate1" vs "1859ISC_Plate3").
#   Instead of requiring an exact name match, this script recursively scans
#   ALL subdirectories under <idat_dir> for *_Grn.idat files and builds a
#   global lookup table: Sentrix_Position -> (barcode, full_idat_path).
#   Each sample's IDAT is resolved by position match alone — the Plate Number
#   column is used only for logical grouping (which output CSV it goes into).
#
# ADDITIONAL: barcode recovery from IDAT filenames
#   When the CSV was saved from Excel, long numeric BeadChip barcodes may be
#   silently truncated to scientific notation (e.g. 208789590164 -> 2.08789E+11).
#   The full barcode is recovered directly from the IDAT filename on disk.
#
# COLUMN PASS-THROUGH
#   All columns from the study samplesheet are carried forward into each
#   per-plate CSV so that participant metadata (sample_id, collection site,
#   age, etc.) is available in combined_samplesheet.csv and ultimately in
#   sample_info_combined.csv.  The meffil-required columns are added or
#   renamed on top of the original columns:
#     Sample_Name  <- sample_id   (meffil uses this as the QC object key)
#     Sex          <- resolved from collected_gender / sex / gender
#     Basename     <- resolved from IDAT files on disk
#     Plate_ID     <- plate_number value
#
# Usage:
#   Rscript split_samplesheet.r <study_samplesheet.csv> <idat_dir> <out_dir>

suppressPackageStartupMessages(library(methods))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript split_samplesheet.r <study_samplesheet.csv> <idat_dir> <out_dir>")
}

ss_path  <- args[1]
idat_dir <- args[2]
out_dir  <- args[3]

cat("\n", rep("=", 60), "\n", sep = "")
cat("SPLITTING STUDY SAMPLESHEET BY PLATE\n")
cat(rep("=", 60), "\n\n")
cat("Study samplesheet:", ss_path, "\n")
cat("IDAT directory   :", idat_dir, "\n\n")

if (!file.exists(ss_path)) stop("Samplesheet not found: ", ss_path)
if (!dir.exists(idat_dir)) stop("IDAT directory not found: ", idat_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ss <- read.csv(ss_path, stringsAsFactors = FALSE, check.names = FALSE,
               colClasses = "character")

# Normalise column names: lowercase, spaces/dots -> underscores
colnames(ss) <- trimws(tolower(gsub("[. ]", "_", colnames(ss))))

required <- c("sample_id", "sentrix_position", "plate_number")
missing  <- setdiff(required, colnames(ss))
if (length(missing) > 0) {
  stop("Study samplesheet missing columns: ", paste(missing, collapse = ", "),
       "\nFound: ", paste(colnames(ss), collapse = ", "))
}

ss$sentrix_position <- toupper(trimws(ss$sentrix_position))
ss$sample_id        <- trimws(ss$sample_id)

# Resolve sex column
sex_col <- intersect(c("collected_gender", "sex", "gender"), colnames(ss))[1]
if (!is.na(sex_col)) {
  raw_sex <- tolower(trimws(ss[[sex_col]]))
  ss$Sex  <- ifelse(raw_sex %in% c("male",   "m"), "M",
             ifelse(raw_sex %in% c("female", "f"), "F", "NA"))
} else {
  ss$Sex <- "NA"
}

# Separate unassigned rows
unassigned <- ss[is.na(ss$plate_number) | trimws(ss$plate_number) == "", ]
assigned   <- ss[!is.na(ss$plate_number) & trimws(ss$plate_number) != "", ]

if (nrow(unassigned) > 0) {
  cat("WARNING:", nrow(unassigned),
      "rows have no Plate Number and will be excluded from analysis.\n")
  ua_file <- file.path(out_dir, "unassigned_samples.csv")
  write.csv(unassigned, ua_file, row.names = FALSE, quote = TRUE)
  cat("  Unassigned rows written to:", ua_file, "\n\n")
}

# ============================================================
# BUILD GLOBAL IDAT LOOKUP — scan ALL subdirs under idat_dir
# ============================================================
# Discover every *_Grn.idat file recursively.
# Build: sentrix_position -> list(barcode, full_basename_path)
# "basename" here = full path minus the _Grn.idat / _Red.idat suffix,
# i.e. what meffil expects as the Basename column.
# ============================================================
cat("Scanning IDAT directory recursively for *_Grn.idat files...\n")

all_grn <- list.files(idat_dir, pattern = "_Grn\\.idat$",
                       ignore.case = TRUE, full.names = TRUE, recursive = TRUE)

if (length(all_grn) == 0) {
  stop("No *_Grn.idat files found anywhere under: ", idat_dir)
}

cat("Total *_Grn.idat files found:", length(all_grn), "\n")

# Parse each path: extract barcode and position from filename
# Expected pattern: <barcode>_<SentrixPosition>_Grn.idat
# e.g. 208788870030_R01C01_Grn.idat
global_lookup <- list()  # key = "R01C01", value = list(barcode=..., path=...)
parse_warnings <- character(0)

for (fp in all_grn) {
  fn    <- basename(fp)
  parts <- strsplit(fn, "_")[[1]]
  # Need at least: <barcode> _ <RxxCxx> _ Grn.idat  → 3 parts minimum
  if (length(parts) < 3) {
    parse_warnings <- c(parse_warnings,
      sprintf("  Skipping unrecognised filename pattern: %s", fp))
    next
  }
  bc  <- parts[1]
  pos <- toupper(parts[2])  # e.g. R01C01
  # Basename = full path stripped of _Grn.idat
  bn  <- sub("_Grn\\.idat$", "", fp, ignore.case = TRUE)

  if (pos %in% names(global_lookup)) {
    # Duplicate position: means two plates have the same chip layout —
    # this is expected and fine; we resolve per-plate below.
    # Store as a list of hits.
    global_lookup[[pos]] <- c(global_lookup[[pos]],
                               list(list(barcode = bc, basename = bn, file = fp)))
  } else {
    global_lookup[[pos]] <- list(list(barcode = bc, basename = bn, file = fp))
  }
}

if (length(parse_warnings) > 0) {
  cat("Filename parse warnings:\n")
  for (w in parse_warnings) cat(w, "\n")
  cat("\n")
}

# Report discovered plate directories (top-level subdirs that contained IDATs)
cat("Plate directories found under idat_dir:\n")
for (pd in sort(unique(sapply(all_grn, function(f) {
  # Walk up from the idat file until we reach a direct child of idat_dir
  rel <- sub(paste0("^", normalizePath(idat_dir), .Platform$file.sep), "", normalizePath(f))
  strsplit(rel, .Platform$file.sep)[[1]][1]
})))) {
  cat("  -", pd, "\n")
}
cat("\n")

# ============================================================
# SPLIT BY PLATE
# ============================================================
plates <- unique(assigned$plate_number)
cat("Plates found in samplesheet:", length(plates), "\n\n")

for (pl in plates) {
  sub_ss <- assigned[assigned$plate_number == pl, ]
  cat("--- Plate:", pl, "---\n")
  cat("  Samples in samplesheet:", nrow(sub_ss), "\n")

  # For each sample, resolve its IDAT basename from the global lookup
  resolved_barcode  <- character(nrow(sub_ss))
  resolved_basename <- character(nrow(sub_ss))
  unresolved_pos    <- character(0)

  for (i in seq_len(nrow(sub_ss))) {
    pos   <- sub_ss$sentrix_position[i]
    hits  <- global_lookup[[pos]]

    if (is.null(hits) || length(hits) == 0) {
      unresolved_pos <- c(unresolved_pos, pos)
      resolved_barcode[i]  <- NA_character_
      resolved_basename[i] <- NA_character_
      next
    }

    if (length(hits) == 1) {
      # Unambiguous
      resolved_barcode[i]  <- hits[[1]]$barcode
      resolved_basename[i] <- hits[[1]]$basename
    } else {
      # Multiple plates share this position — use the beadchip_barcode column
      # from the samplesheet to disambiguate (if present and not corrupted).
      bc_col <- intersect(c("beadchip_barcode", "beadchip", "barcode", "chip_id"),
                          colnames(sub_ss))[1]
      csv_bc <- if (!is.na(bc_col)) trimws(sub_ss[[bc_col]][i]) else ""

      matched_hit <- NULL
      if (nchar(csv_bc) > 0 && !grepl("E[+-]", csv_bc, ignore.case = TRUE)) {
        # CSV barcode looks intact (not scientific notation) — use it
        matched_hit <- Filter(function(h) h$barcode == csv_bc, hits)
        matched_hit <- if (length(matched_hit) > 0) matched_hit[[1]] else NULL
      }

      if (is.null(matched_hit)) {
        # Barcode unavailable or ambiguous — fall back to first hit and warn
        matched_hit <- hits[[1]]
        cat(sprintf("  WARNING: position %s matches %d IDAT files; ",
                    pos, length(hits)))
        cat(sprintf("using %s (set BeadChip_Barcode column to disambiguate)\n",
                    matched_hit$barcode))
      }

      resolved_barcode[i]  <- matched_hit$barcode
      resolved_basename[i] <- matched_hit$basename
    }
  }

  if (length(unresolved_pos) > 0) {
    stop("Could not find IDAT files for ", length(unresolved_pos),
         " positions in plate ", pl, ":\n  ",
         paste(unresolved_pos, collapse = ", "),
         "\nMake sure these IDAT files exist somewhere under: ", idat_dir,
         "\nand follow the naming pattern <barcode>_<SentrixPosition>_Grn.idat")
  }

  # ── Build per-plate output CSV ──────────────────────────────────────────
  # Start from the full original samplesheet row so no metadata is lost.
  # Then set / overwrite the four meffil-required columns:
  #   Sample_Name  — the study participant ID (meffil uses this as the QC
  #                  object key; it must be unique within the plate)
  #   Sex          — normalised M/F/NA (already computed above in ss$Sex)
  #   Basename     — full path prefix resolved from IDAT files on disk
  #   Plate_ID     — plate identifier for grouping in downstream steps
  #
  # Columns that are internal artefacts of this script (normalised lowercase
  # originals, the raw barcode as read from CSV, real_barcode from disk) are
  # dropped to avoid confusion with the canonical output columns.
  out <- sub_ss

  # Attach resolved IDAT paths
  out$real_barcode <- resolved_barcode
  out$Basename     <- resolved_basename

  # Set meffil-required columns (canonical mixed-case names meffil expects)
  out$Sample_Name  <- out$sample_id
  out$Plate_ID     <- pl
  # Sex was already set as ss$Sex above and is present in out

  # Drop internal / redundant columns that would cause confusion downstream:
  #   sample_id        — kept as Sample_Name (same value, different case)
  #   plate_number     — kept as Plate_ID
  #   sentrix_position — array position, not needed after IDAT resolution
  #   real_barcode     — useful for debugging but not a downstream input
  # Keep everything else (participant metadata, collection site, age, etc.)
  drop_cols <- intersect(c("sample_id", "plate_number", "sentrix_position",
                            "real_barcode"),
                         colnames(out))
  out <- out[, setdiff(colnames(out), drop_cols), drop = FALSE]

  # Ensure meffil-required columns are first for readability
  priority_cols  <- intersect(c("Sample_Name", "Sex", "Basename", "Plate_ID"),
                               colnames(out))
  remaining_cols <- setdiff(colnames(out), priority_cols)
  out <- out[, c(priority_cols, remaining_cols), drop = FALSE]

  outfile <- file.path(out_dir, paste0(pl, "_samplesheet.csv"))
  write.csv(out, outfile, row.names = FALSE, quote = TRUE)

  cat("  Resolved barcodes   :", length(unique(resolved_barcode)), "unique chip(s)\n")
  cat("  Columns in output   :", ncol(out), "\n")
  cat("  Wrote", nrow(out), "samples ->", basename(outfile), "\n\n")
}

cat(rep("=", 60), "\n", sep = "")
cat("SPLIT COMPLETE\n")
cat(rep("=", 60), "\n\n")
