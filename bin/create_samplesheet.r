#!/usr/bin/env Rscript
# create_samplesheet.r
# Creates meffil samplesheet from all IDAT files in a plate folder
# Usage: Rscript create_samplesheet.r <plate_manifest.txt> <out_samplesheet.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_samplesheet.r <plate_manifest.txt> <out_samplesheet.csv>")
}

plate_manifest <- args[1]
out_samplesheet <- args[2]

cat("\n", rep("=", 60), "\n", sep = "")
cat("CREATING MEFFIL SAMPLESHEET FOR PLATE\n")
cat(rep("=", 60), "\n\n")

plate_id <- sub("\\.txt$", "", basename(plate_manifest))
cat("Plate:", plate_id, "\n")

# Read IDAT file paths from manifest
idat_files <- readLines(plate_manifest)
if (length(idat_files) == 0) {
  stop("Empty plate manifest: ", plate_manifest)
}

cat("Found", length(idat_files), "IDAT files in manifest\n")

# Get unique basenames (remove _Grn.idat and _Red.idat)
basenames <- unique(sub("_(Grn|Red)\\.idat$", "", idat_files, ignore.case = TRUE))
cat("Unique samples:", length(basenames), "\n")

# Show some example files to understand the naming pattern
cat("\nFirst few sample basenames:\n")
for (i in 1:min(3, length(basenames))) {
  cat("  ", basename(basenames[i]), "\n")
}

# Create simple samplesheet manually
samplesheet <- data.frame(
  Sample_Name = basename(basenames),
  Sex = "NA",  # Important: Use string "NA" not R's NA value
  Basename = basenames,
  Plate_ID = plate_id,
  stringsAsFactors = FALSE
)

# Ensure unique Sample_Name values
if (any(duplicated(samplesheet$Sample_Name))) {
  cat("Making duplicate Sample_Name values unique...\n")
  samplesheet$Sample_Name <- make.unique(samplesheet$Sample_Name, sep = "_")
}

cat("Created samplesheet with", nrow(samplesheet), "samples\n")

# Save samplesheet
write.csv(samplesheet, out_samplesheet, row.names = FALSE, quote = TRUE)

cat("\nSamplesheet saved to:", out_samplesheet, "\n")
cat("\nSamplesheet preview:\n")
print(head(samplesheet[, c("Sample_Name", "Sex", "Plate_ID")], 3))

cat("\n", rep("=", 60), "\n", sep = "")
cat("SAMPLESHEET CREATION COMPLETED\n")
cat(rep("=", 60), "\n\n")

# Return info for verification
cat("Summary:\n")
cat("  Samples:", nrow(samplesheet), "\n")
cat("  Columns:", ncol(samplesheet), "\n")
cat("  Unique Sample_Name values:", length(unique(samplesheet$Sample_Name)), "\n")
cat("  Sex values (should all be 'NA'):\n")
print(table(samplesheet$Sex, useNA = "ifany"))
