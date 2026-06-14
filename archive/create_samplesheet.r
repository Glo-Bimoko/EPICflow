#!/usr/bin/env Rscript
# create_samplesheet.r
# ====================
# Builds a meffil-format samplesheet for one plate by joining:
#   1. IDAT file paths from the plate manifest
#   2. The study samplesheet (Sample_ID, BeadChip_Barcode, Sentrix_Position,
#      Plate_Number, Well_Position, Collected_Gender)
#
# The study samplesheet is the authoritative source of Sample_ID.
# Without it, Sample_Name defaults to the Sentrix barcode string
# (e.g. 208789590164_R01C01) which loses participant identity throughout
# the entire pipeline including concordance outputs and the id_map lookup.
#
# Usage:
#   Rscript create_samplesheet.r \
#       <plate_manifest.txt> \
#       <study_samplesheet.csv> \
#       <out_samplesheet.csv>
#
# Arguments:
#   plate_manifest.txt   : list of IDAT file paths for this plate (from group_by_plate.r)
#   study_samplesheet.csv: study samplesheet with Sample_ID column.
#                          Pass "" or "NONE" to fall back to Sentrix basenames.
#   out_samplesheet.csv  : output path for the meffil samplesheet
#
# Study samplesheet required columns (case-insensitive):
#   Sample_ID            : participant identifier (e.g. G001, G003, ...)
#   BeadChip_Barcode     : chip barcode (e.g. 208789590164)
#   Sentrix_Position     : well position on chip (e.g. R01C01)
#
# Optional columns passed through to the meffil samplesheet:
#   Collected_Gender / Sex : used by meffil sex prediction check
#   Plate_Number, Well_Position : retained as metadata

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript create_samplesheet.r ",
    "<plate_manifest.txt> <study_samplesheet.csv> <out_samplesheet.csv>"
  )
}

plate_manifest_path  <- args[1]
study_ss_raw         <- args[2]
out_samplesheet_path <- args[3]

# Treat empty string or "NONE" as "no samplesheet provided"
study_ss_path <- if (nchar(trimws(study_ss_raw)) == 0 || toupper(trimws(study_ss_raw)) == "NONE") {
  NULL
} else {
  study_ss_raw
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("CREATING MEFFIL SAMPLESHEET FOR PLATE\n")
cat(rep("=", 60), "\n\n")

plate_id <- sub("\\.txt$", "", basename(plate_manifest_path))
cat("Plate            :", plate_id, "\n")
cat("Study samplesheet:", if (!is.null(study_ss_path)) study_ss_path else "(not provided — using Sentrix basenames)\n")

# ============================================================
# STEP 1: READ IDAT PATHS FROM PLATE MANIFEST
# ============================================================
cat("\n--- STEP 1: Reading plate manifest ---\n")
idat_files <- readLines(plate_manifest_path)
idat_files <- idat_files[nchar(trimws(idat_files)) > 0]
if (length(idat_files) == 0) stop("Empty plate manifest: ", plate_manifest_path)
cat("IDAT files in manifest:", length(idat_files), "\n")

# Derive unique Basenames (full path without _Grn.idat / _Red.idat)
basenames <- unique(sub("_(Grn|Red)\\.idat$", "", idat_files, ignore.case = TRUE))
cat("Unique samples   :", length(basenames), "\n")

# Parse BeadChip_Barcode and Sentrix_Position from filename
# Format: /path/to/<Barcode>_<Position>_Grn.idat
barcode_pos <- do.call(rbind, lapply(basenames, function(bn) {
  fname <- basename(bn)
  parts <- strsplit(fname, "_")[[1]]
  if (length(parts) < 2) {
    cat("  Warning: unexpected filename format:", fname, "\n")
    return(data.frame(Barcode = fname, Position = NA_character_,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    Barcode  = parts[1],
    Position = paste(parts[-1], collapse = "_"),  # handles R01C01 etc.
    stringsAsFactors = FALSE
  )
}))
barcode_pos$Basename <- basenames

cat("\nFirst few (Barcode | Position | Basename):\n")
for (i in seq_len(min(3, nrow(barcode_pos)))) {
  cat(sprintf("  %s | %s | %s\n",
              barcode_pos$Barcode[i], barcode_pos$Position[i],
              basename(barcode_pos$Basename[i])))
}

# ============================================================
# STEP 2: JOIN WITH STUDY SAMPLESHEET (if provided)
# ============================================================
cat("\n--- STEP 2: Joining with study samplesheet ---\n")

normalise_colname <- function(x) gsub("[^a-zA-Z0-9]", "_", tolower(trimws(x)))

if (!is.null(study_ss_path)) {
  if (!file.exists(study_ss_path)) stop("Study samplesheet not found: ", study_ss_path)

  study_ss <- read.csv(study_ss_path, stringsAsFactors = FALSE)
  colnames(study_ss) <- trimws(colnames(study_ss))
  cat("Study samplesheet rows :", nrow(study_ss), "\n")
  cat("Study samplesheet cols :", paste(colnames(study_ss), collapse = ", "), "\n")

  # Find the Sample_ID column (flexible naming)
  norm_cols   <- normalise_colname(colnames(study_ss))
  id_col      <- colnames(study_ss)[which(norm_cols == "sample_id")[1]]
  barcode_col <- colnames(study_ss)[which(norm_cols %in% c("beadchip_barcode", "beadchip", "barcode", "chip_id"))[1]]
  pos_col     <- colnames(study_ss)[which(norm_cols %in% c("sentrix_position", "sentrix_pos", "position", "well"))[1]]

  if (is.na(id_col))      stop("Cannot find Sample_ID column in study samplesheet. Found: ",
                                paste(colnames(study_ss), collapse = ", "))
  if (is.na(barcode_col)) stop("Cannot find BeadChip_Barcode column in study samplesheet.")
  if (is.na(pos_col))     stop("Cannot find Sentrix_Position column in study samplesheet.")

  cat("Matched columns -> Sample_ID:", id_col, "| Barcode:", barcode_col, "| Position:", pos_col, "\n")

  # Normalise barcode and position for joining (trim whitespace, uppercase position)
  study_ss$join_barcode  <- trimws(as.character(study_ss[[barcode_col]]))
  study_ss$join_position <- toupper(trimws(study_ss[[pos_col]]))
  barcode_pos$join_barcode  <- trimws(barcode_pos$Barcode)
  barcode_pos$join_position <- toupper(trimws(barcode_pos$Position))

  merged <- merge(
    barcode_pos,
    study_ss,
    by = c("join_barcode", "join_position"),
    all.x = TRUE
  )

  n_matched   <- sum(!is.na(merged[[id_col]]))
  n_unmatched <- sum(is.na(merged[[id_col]]))
  cat("Samples matched in study samplesheet:", n_matched, "\n")

  if (n_unmatched > 0) {
    cat("⚠  Unmatched IDAT files (no entry in study samplesheet):", n_unmatched, "\n")
    unmatched_rows <- merged[is.na(merged[[id_col]]), ]
    for (i in seq_len(nrow(unmatched_rows))) {
      cat(sprintf("   Barcode %s Position %s → Sample_Name set to Sentrix basename\n",
                  unmatched_rows$join_barcode[i], unmatched_rows$join_position[i]))
    }
  }

  # Set Sample_Name: use Sample_ID where matched, fall back to Sentrix basename
  merged$Sample_Name <- ifelse(
    !is.na(merged[[id_col]]) & nchar(trimws(merged[[id_col]])) > 0,
    trimws(merged[[id_col]]),
    paste0(merged$Barcode, "_", merged$Position)
  )

  # Resolve sex column (meffil expects "Sex" with M/F/"NA")
  sex_col <- colnames(study_ss)[which(normalise_colname(colnames(study_ss)) %in%
                                       c("collected_gender", "sex", "gender"))[1]]
  if (!is.na(sex_col)) {
    merged$Sex <- toupper(trimws(as.character(merged[[sex_col]])))
    # Normalise: Male->M, Female->F, Unknown/empty->"NA"
    merged$Sex <- ifelse(merged$Sex %in% c("M", "MALE"),   "M",
                  ifelse(merged$Sex %in% c("F", "FEMALE"), "F", "NA"))
  } else {
    merged$Sex <- "NA"
  }

  # Pull extra metadata columns if present
  plate_col    <- colnames(study_ss)[which(normalise_colname(colnames(study_ss)) == "plate_number")[1]]
  well_col     <- colnames(study_ss)[which(normalise_colname(colnames(study_ss)) == "well_position")[1]]

  merged$Plate_ID    <- plate_id
  merged$Plate_Number <- if (!is.na(plate_col))    trimws(as.character(merged[[plate_col]]))    else plate_id
  merged$Well_Position <- if (!is.na(well_col))    trimws(as.character(merged[[well_col]]))    else NA_character_

  samplesheet <- data.frame(
    Sample_Name    = merged$Sample_Name,
    Sex            = merged$Sex,
    Basename       = merged$Basename,
    Plate_ID       = merged$Plate_ID,
    Plate_Number   = merged$Plate_Number,
    Well_Position  = merged$Well_Position,
    Sentrix_ID     = merged$Barcode,
    Sentrix_Position = merged$Position,
    stringsAsFactors = FALSE
  )

} else {
  # No study samplesheet — fall back to Sentrix basename as Sample_Name
  cat("No study samplesheet provided. Sample_Name = Sentrix barcode_position.\n")
  cat("⚠  Downstream outputs (concordance, id_map) will use Sentrix IDs not participant IDs.\n")

  samplesheet <- data.frame(
    Sample_Name      = paste0(barcode_pos$Barcode, "_", barcode_pos$Position),
    Sex              = "NA",
    Basename         = barcode_pos$Basename,
    Plate_ID         = plate_id,
    Plate_Number     = plate_id,
    Well_Position    = NA_character_,
    Sentrix_ID       = barcode_pos$Barcode,
    Sentrix_Position = barcode_pos$Position,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# STEP 3: VALIDATE AND WRITE
# ============================================================
cat("\n--- STEP 3: Validating samplesheet ---\n")

if (any(duplicated(samplesheet$Sample_Name))) {
  dups <- samplesheet$Sample_Name[duplicated(samplesheet$Sample_Name)]
  cat("⚠  Duplicate Sample_Name values detected:", paste(unique(dups), collapse = ", "), "\n")
  cat("   Making them unique by appending suffix...\n")
  samplesheet$Sample_Name <- make.unique(samplesheet$Sample_Name, sep = "_dup")
}

cat("Final samplesheet:", nrow(samplesheet), "samples\n")
cat("\nPreview (first 5 rows):\n")
print(head(samplesheet[, c("Sample_Name", "Sex", "Sentrix_ID", "Sentrix_Position", "Plate_ID")], 5))

write.csv(samplesheet, out_samplesheet_path, row.names = FALSE, quote = TRUE)
cat("\n✓ Samplesheet saved to:", out_samplesheet_path, "\n")

cat("\n", rep("=", 60), "\n", sep = "")
cat("SAMPLESHEET CREATION COMPLETED\n")
cat(rep("=", 60), "\n")
cat("  Plate    :", plate_id, "\n")
cat("  Samples  :", nrow(samplesheet), "\n")
cat("  ID source:", if (!is.null(study_ss_path)) "study samplesheet (Sample_ID)" else "Sentrix basename (fallback)", "\n")
