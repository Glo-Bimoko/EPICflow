library(meffil)

# ── Find the per-plate samplesheets written by SPLIT_SAMPLESHEET ──────────
ss_dir <- "/home/eiko/Desktop/EPICflow/results/samplesheets"
ss_files <- list.files(ss_dir, pattern="Plate_.*_samplesheet\\.csv$",
                        full.names=TRUE)
# Also try the work dir samplesheets output
if (length(ss_files)==0)
  ss_files <- list.files("/home/eiko/Desktop/EPICflow/samplesheets",
                          pattern="Plate_.*\\.csv$", full.names=TRUE)
cat("Per-plate samplesheets found:", length(ss_files), "\n")
print(ss_files)

# Load and combine all per-plate samplesheets
all_ss <- do.call(rbind, lapply(ss_files, read.csv, stringsAsFactors=FALSE))
cat("\nCombined samplesheet rows:", nrow(all_ss), "\n")
cat("Columns:", colnames(all_ss), "\n")

# Build Sample_Name -> Basename lookup
bn_map <- setNames(all_ss$Basename, all_ss$Sample_Name)
cat("\nBasename sample (first 3):\n")
print(head(bn_map, 3))

# Load QC objects and check males
work    <- "/home/eiko/Desktop/EPICflow/work/1c/dc27dcac93deeb2f38f01731dca6ab/qc_data"
qc_files <- list.files(work, pattern="_qc_objects\\.rds$", full.names=TRUE)
all_qc  <- list()
for (f in qc_files) all_qc <- c(all_qc, readRDS(f))
male_idx <- which(sapply(all_qc, function(x) x$predicted.sex) == "M")

# Verify first 6 males get correct basenames from the map
cat("\nBasename resolution for first 6 males:\n")
for (i in male_idx[1:6]) {
  sname <- all_qc[[i]]$sample.name
  base  <- bn_map[sname]
  grn_ok <- !is.na(base) && file.exists(paste0(base, "_Grn.idat"))
  cat(sprintf("  %-10s -> %s  [Grn exists: %s]\n", sname,
              ifelse(is.na(base), "NOT IN MAP", basename(base)), grn_ok))
}