#!/usr/bin/env Rscript
# group_by_plate.r
# Groups IDAT files by existing plate folders
# Usage: Rscript group_by_plate.r <idat_dir> <out_dir> <summary_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript group_by_plate.r <idat_dir> <out_dir> <summary_file>")
}

idat_dir <- args[1]
out_dir <- args[2]
summary_file <- args[3]

# Validate inputs
if (!dir.exists(idat_dir)) {
  stop("IDAT directory does not exist: ", idat_dir)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("GROUPING IDATs BY EXISTING PLATE FOLDERS\n")
cat(rep("=", 60), "\n\n")
cat("IDAT directory:", idat_dir, "\n")
cat("Output directory:", out_dir, "\n\n")

# Find all plate folders (folders containing IDAT files)
plate_folders <- list.dirs(idat_dir, full.names = TRUE, recursive = FALSE)
plate_folders <- plate_folders[grepl("plate", basename(plate_folders), ignore.case = TRUE)]

if (length(plate_folders) == 0) {
  # If no plate folders found, look for any folders containing IDAT files
  cat("No plate folders found. Searching for any folders with IDAT files...\n")
  all_folders <- list.dirs(idat_dir, full.names = TRUE, recursive = TRUE)
  idat_folders <- c()
  
  for (folder in all_folders) {
    idat_files <- list.files(folder, pattern = "\\.idat$", ignore.case = TRUE, full.names = TRUE)
    if (length(idat_files) > 0) {
      idat_folders <- c(idat_folders, folder)
    }
  }
  plate_folders <- unique(idat_folders)
}

cat("Found", length(plate_folders), "plate folder(s):\n")
for (folder in plate_folders) {
  idat_count <- length(list.files(folder, pattern = "\\.idat$", ignore.case = TRUE, full.names = TRUE))
  cat("  -", basename(folder), "(", idat_count, "IDAT files)\n")
}

# Create plate manifests
plate_summary <- data.frame(
  Plate = character(),
  Plate_Folder = character(),
  N_IDATs = integer(),
  N_Samples = integer(),
  Chip_IDs = character(),
  stringsAsFactors = FALSE
)

for (plate_folder in plate_folders) {
  plate_name <- basename(plate_folder)
  
  # Find all IDAT files in this plate folder (recursively)
  idat_files <- list.files(plate_folder, pattern = "\\.idat$", 
                          full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  
  if (length(idat_files) == 0) {
    cat("Warning: No IDAT files found in", plate_folder, "\n")
    next
  }
  
  # Write manifest
  manifest_file <- file.path(out_dir, paste0(plate_name, ".txt"))
  writeLines(idat_files, manifest_file)
  
  # Get unique samples and chips
  basenames <- unique(sub("_(Grn|Red)\\.idat$", "", basename(idat_files), ignore.case = TRUE))
  chip_ids <- unique(sapply(basenames, function(bn) sub("_.*$", "", bn)))
  
  cat("Plate", plate_name, ":\n")
  cat("  IDAT files:", length(idat_files), "\n")
  cat("  Samples:", length(basenames), "\n")
  cat("  Chips:", length(chip_ids), "(", paste(chip_ids, collapse = ", "), ")\n")
  cat("  Manifest:", basename(manifest_file), "\n\n")
  
  # Add to summary
  plate_summary <- rbind(plate_summary, data.frame(
    Plate = plate_name,
    Plate_Folder = plate_folder,
    N_IDATs = length(idat_files),
    N_Samples = length(basenames),
    Chip_IDs = paste(chip_ids, collapse = ";"),
    stringsAsFactors = FALSE
  ))
}

# Write summary
write.csv(plate_summary, summary_file, row.names = FALSE)

cat(rep("=", 60), "\n", sep = "")
cat("GROUPING COMPLETED\n")
cat(rep("=", 60), "\n\n")
cat("Summary:\n")
cat("  Total plates:", nrow(plate_summary), "\n")
cat("  Total IDAT files:", sum(plate_summary$N_IDATs), "\n")
cat("  Total samples:", sum(plate_summary$N_Samples), "\n")
cat("  Manifests written to:", out_dir, "\n")
cat("  Summary saved to:", summary_file, "\n")
