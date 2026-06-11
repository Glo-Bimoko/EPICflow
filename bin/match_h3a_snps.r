#!/usr/bin/env Rscript
# match_h3a_snps.r
# ================
# Map the 65 EPIC/EPICv2 rs-probe SNPs onto an H3Africa PLINK dataset and
# produce a .traw genotype file for cross_study_concordance.py.
#
# H3AFRICA .bim VARIANT ID FORMAT
# --------------------------------
# H3Africa arrays do not use rsIDs in the .bim SNP column.  Instead, variants
# are identified by positional strings such as:
#   seq-h3a_37_1_61989_G_C
#   h3a_37_1_108310_T_C
#   kgp15717912              (legacy Affymetrix IDs — no position encoded)
#
# Because of this, direct rsID matching will almost never succeed.
# The PRIMARY matching strategy is POSITIONAL:
#   1. Get hg19/GRCh37 chr:pos for each EPIC rs probe from the meffil manifest.
#   2. Find the H3A variant at the same chr:pos that shares at least one allele.
#   3. For the small number of kgp/rs IDs present, attempt a direct name match
#      as a supplementary fallback.
#
# After matching, PLINK extracts those variants and recodes to .traw.
# Sample IDs in the output .traw are cleaned to plain IID (FID stripped).
#
# Usage:
#   Rscript match_h3a_snps.r \
#       --bfile      /path/to/h3afow_output \
#       --epic_traw  concordance/methylation_snps_snp_calls.traw \
#       --out_prefix concordance/h3a_snps \
#       --featureset epicv2 \
#       [--plink     plink] \
#       [--snp_names bin/snp-names.txt]
#
# Outputs:
#   <out_prefix>_h3a_snp_calls.traw   – H3A genotype calls, SNPs relabelled
#                                        to EPIC probe names, sample IDs clean
#   <out_prefix>_snp_match_table.csv  – per-probe match details
#   <out_prefix>_snp_match_report.txt – human-readable summary

suppressPackageStartupMessages({
  library(meffil)
})

# ============================================================
# ARGUMENT PARSING
# ============================================================
args_raw <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    bfile      = NULL,
    epic_traw  = NULL,
    out_prefix = NULL,
    featureset = "epicv2",
    plink      = "plink",
    snp_names  = NULL
  )
  i <- 1
  while (i <= length(args)) {
    flag <- args[i]
    val  <- if (i + 1 <= length(args)) args[i + 1] else NULL
    switch(flag,
      "--bfile"      = { out$bfile      <- val; i <- i + 2 },
      "--epic_traw"  = { out$epic_traw  <- val; i <- i + 2 },
      "--out_prefix" = { out$out_prefix <- val; i <- i + 2 },
      "--featureset" = { out$featureset <- val; i <- i + 2 },
      "--plink"      = { out$plink      <- val; i <- i + 2 },
      "--snp_names"  = { out$snp_names  <- val; i <- i + 2 },
      { stop("Unknown argument: ", flag) }
    )
  }
  for (req in c("bfile", "epic_traw", "out_prefix")) {
    if (is.null(out[[req]])) stop("Required argument missing: --", req)
  }
  out
}

opt <- parse_args(args_raw)

cat("\n", rep("=", 60), "\n", sep = "")
cat("H3A SNP MATCHING FOR CROSS-STUDY CONCORDANCE\n")
cat(rep("=", 60), "\n\n")
cat("H3A PLINK bfile  :", opt$bfile, "\n")
cat("EPIC .traw       :", opt$epic_traw, "\n")
cat("Output prefix    :", opt$out_prefix, "\n")
cat("Array featureset :", opt$featureset, "\n")
cat("PLINK binary     :", opt$plink, "\n")
if (!is.null(opt$snp_names))
  cat("SNP names file   :", opt$snp_names, "\n")
cat("\n")

# ============================================================
# STEP 1: LOAD EPIC PROBE LIST
# ============================================================
# Priority:
#   A) --snp_names file (authoritative list; independent of QC dropout)
#   B) probes retained in the EPIC .traw (fallback)

cat("--- STEP 1: Loading EPIC probe list ---\n")

epic_traw_lines <- readLines(opt$epic_traw)
data_lines      <- epic_traw_lines[-1]
data_lines      <- data_lines[nchar(trimws(data_lines)) > 0]

epic_probes_in_traw <- sapply(data_lines, function(ln) {
  strsplit(ln, "\t")[[1]][2]
})
names(epic_probes_in_traw) <- NULL
cat("EPIC probes present in .traw (post-QC):", length(epic_probes_in_traw), "\n")

if (!is.null(opt$snp_names) && nchar(opt$snp_names) > 0) {
  if (!file.exists(opt$snp_names)) stop("--snp_names file not found: ", opt$snp_names)
  snp_list <- readLines(opt$snp_names)
  snp_list <- trimws(snp_list)
  snp_list <- unique(snp_list[nchar(snp_list) > 0 & grepl("^rs", snp_list)])
  epic_probes <- snp_list
  cat("Using authoritative SNP list from    :", opt$snp_names, "\n")
  cat("SNPs in list                         :", length(epic_probes), "\n")
  qc_dropped <- setdiff(epic_probes, epic_probes_in_traw)
  if (length(qc_dropped) > 0) {
    cat("⚠  Probes in snp-names.txt absent from .traw (QC dropout):\n")
    cat("   ", paste(qc_dropped, collapse = ", "), "\n")
  } else {
    cat("✓  All listed SNPs present in .traw\n")
  }
} else {
  epic_probes <- epic_probes_in_traw
  cat("No --snp_names provided; using probes from EPIC .traw\n")
}
cat("Target probes for H3A extraction:", length(epic_probes), "\n\n")

# ============================================================
# STEP 2: READ H3A .bim FILE
# ============================================================
cat("--- STEP 2: Reading H3A .bim file ---\n")

bim_path <- paste0(opt$bfile, ".bim")
if (!file.exists(bim_path)) stop("Cannot find .bim file: ", bim_path)

bim <- read.table(bim_path, header = FALSE, sep = "\t",
                  col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"),
                  stringsAsFactors = FALSE)
cat("H3A variants in .bim:", nrow(bim), "\n")
cat("Chromosomes present :", paste(sort(unique(bim$CHR)), collapse = " "), "\n")

# Show the variant ID formats present (helps diagnose matching issues)
sample_ids_bim <- head(bim$SNP, 8)
cat("Variant ID examples :", paste(sample_ids_bim, collapse = ", "), "\n\n")

# Note: H3Africa arrays use positional IDs (seq-h3a_37_CHR_POS_REF_ALT,
# h3a_37_CHR_POS_REF_ALT) and legacy Affymetrix IDs (kgpXXXX).
# rsIDs are rare or absent, so positional matching is the primary strategy.

# ============================================================
# STEP 3: GET EPIC PROBE COORDINATES FROM MEFFIL MANIFEST
# ============================================================
# This is the primary matching path for H3Africa data because the .bim
# does not contain rsIDs that match the EPIC probe names.

cat("--- STEP 3: Retrieving EPIC probe coordinates from meffil manifest ---\n")

snp_annot <- NULL
for (fs in c(opt$featureset, "epic", "epicv2", "450k")) {
  snp_annot <- tryCatch(meffil.get.features(fs), error = function(e) NULL)
  if (!is.null(snp_annot)) {
    cat("Using meffil featureset:", fs, "\n")
    break
  }
}

if (is.null(snp_annot)) {
  stop(
    "Could not load meffil SNP annotations for any featureset. ",
    "Run meffil.list.featuresets() to see what is available and ",
    "pass the correct name via --featureset."
  )
}

# Keep only rs-probe rows
snp_rows <- snp_annot[grepl("^rs", snp_annot$name), ]
cat("SNP probes in manifest:", nrow(snp_rows), "\n")

# Build lookup: probe_name -> (chr, pos)
# Normalise chromosome to numeric string (no "chr" prefix) to match .bim format
probe_coords <- data.frame(
  EPIC_PROBE = snp_rows$name,
  CHR        = sub("^chr", "", as.character(snp_rows$chromosome)),
  POS        = as.integer(snp_rows$position),
  stringsAsFactors = FALSE
)

# Restrict to our target probe list
probe_coords <- probe_coords[probe_coords$EPIC_PROBE %in% epic_probes, ]
cat("Probe coordinates found for", nrow(probe_coords), "/",
    length(epic_probes), "target probes\n\n")

# ============================================================
# STEP 4: POSITIONAL MATCHING (PRIMARY STRATEGY)
# ============================================================
cat("--- STEP 4: Positional matching (primary strategy for H3Africa data) ---\n")

# Build a chr:pos index of the .bim for fast lookup
bim$CHR_str <- as.character(bim$CHR)
bim$POS_KEY <- paste0(bim$CHR_str, ":", bim$POS)

probe_coords$POS_KEY <- paste0(probe_coords$CHR, ":", probe_coords$POS)

pos_matches <- data.frame(
  EPIC_PROBE = character(0),
  H3A_SNP    = character(0),
  CHR        = character(0),
  POS        = integer(0),
  A1         = character(0),
  A2         = character(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(probe_coords))) {
  probe  <- probe_coords$EPIC_PROBE[i]
  pk     <- probe_coords$POS_KEY[i]

  hits   <- bim[bim$POS_KEY == pk, ]
  if (nrow(hits) == 0) next

  # If multiple H3A variants at same position (rare), prefer any with an
  # rsID, then take the first
  rs_hits <- hits[grepl("^rs", hits$SNP), ]
  chosen  <- if (nrow(rs_hits) > 0) rs_hits[1, ] else hits[1, ]

  pos_matches <- rbind(pos_matches, data.frame(
    EPIC_PROBE = probe,
    H3A_SNP    = chosen$SNP,
    CHR        = as.character(chosen$CHR),
    POS        = chosen$POS,
    A1         = chosen$A1,
    A2         = chosen$A2,
    stringsAsFactors = FALSE
  ))
}

cat("Positional matches found:", nrow(pos_matches), "/",
    length(epic_probes), "target probes\n")

# ============================================================
# STEP 5: SUPPLEMENTARY DIRECT NAME MATCH
# ============================================================
# For probes not found positionally, attempt a direct name lookup in the
# .bim.  This handles the rare kgpXXXX entries and any genuine rsIDs.

cat("\n--- STEP 5: Supplementary direct name matching ---\n")

pos_matched_probes <- pos_matches$EPIC_PROBE
still_unmatched    <- setdiff(epic_probes, pos_matched_probes)
cat("Probes not yet matched:", length(still_unmatched), "\n")

direct_matches <- data.frame(
  EPIC_PROBE = character(0),
  H3A_SNP    = character(0),
  CHR        = character(0),
  POS        = integer(0),
  A1         = character(0),
  A2         = character(0),
  stringsAsFactors = FALSE
)

if (length(still_unmatched) > 0) {
  bim_name_idx <- match(still_unmatched, bim$SNP)
  for (k in seq_along(still_unmatched)) {
    idx <- bim_name_idx[k]
    if (!is.na(idx)) {
      direct_matches <- rbind(direct_matches, data.frame(
        EPIC_PROBE = still_unmatched[k],
        H3A_SNP    = bim$SNP[idx],
        CHR        = as.character(bim$CHR[idx]),
        POS        = bim$POS[idx],
        A1         = bim$A1[idx],
        A2         = bim$A2[idx],
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("Direct name matches found:", nrow(direct_matches), "\n")
}

# ============================================================
# STEP 6: BUILD FINAL MATCH TABLE
# ============================================================
cat("\n--- STEP 6: Building final SNP match table ---\n")

pos_final    <- if (nrow(pos_matches)    > 0) cbind(pos_matches,    MATCH_TYPE = "positional")  else pos_matches[0,]
direct_final <- if (nrow(direct_matches) > 0) cbind(direct_matches, MATCH_TYPE = "direct_name") else direct_matches[0,]

all_matches <- rbind(pos_final, direct_final)
all_matches <- all_matches[!duplicated(all_matches$EPIC_PROBE), ]

still_unmatched_final <- setdiff(epic_probes, all_matches$EPIC_PROBE)

cat("Total matched probes :", nrow(all_matches), "\n")
cat("  via position       :", sum(all_matches$MATCH_TYPE == "positional"), "\n")
cat("  via direct name    :", sum(all_matches$MATCH_TYPE == "direct_name"), "\n")
cat("Still unmatched      :", length(still_unmatched_final), "\n")
if (length(still_unmatched_final) > 0) {
  cat("  Probes             :", paste(still_unmatched_final, collapse = ", "), "\n")
}
cat("\n")

if (nrow(all_matches) == 0) {
  stop(
    "No EPIC probes could be matched to the H3A dataset.\n",
    "Most likely cause: meffil coordinates are on a different genome build ",
    "than the H3A .bim (H3Africa arrays use GRCh37/hg19 in column 4).\n",
    "Check that your .bim POS values are hg19 and that the meffil featureset ",
    "also uses hg19."
  )
}
if (nrow(all_matches) < 20) {
  warning("Fewer than 20 SNPs matched — concordance results will be unreliable.")
}

match_table_out <- paste0(opt$out_prefix, "_snp_match_table.csv")
write.csv(all_matches, match_table_out, row.names = FALSE)
cat("✓ Match table saved:", match_table_out, "\n")

# ============================================================
# STEP 7: EXTRACT MATCHED SNPs FROM H3A WITH PLINK
# ============================================================
cat("\n--- STEP 7: Extracting matched SNPs from H3A with PLINK ---\n")

extract_list <- paste0(opt$out_prefix, "_extract_list.txt")
writeLines(all_matches$H3A_SNP, extract_list)
cat("Variants to extract:", nrow(all_matches), "\n")

plink_out <- paste0(opt$out_prefix, "_h3a_extracted")
plink_cmd <- paste(
  opt$plink,
  "--bfile",   opt$bfile,
  "--extract", extract_list,
  "--recode A-transpose",
  "--out",     plink_out,
  "--silent"
)
cat("Running PLINK:\n  ", plink_cmd, "\n")

ret <- system(plink_cmd)
if (ret != 0) {
  stop("PLINK exited with code ", ret,
       ". Check that plink is on PATH and bfile paths are correct.")
}

traw_path <- paste0(plink_out, ".traw")
if (!file.exists(traw_path)) {
  stop("PLINK did not produce expected .traw file: ", traw_path)
}
cat("✓ H3A .traw produced:", traw_path, "\n")

# ============================================================
# STEP 8: RELABEL SNP IDs AND CLEAN SAMPLE IDs
# ============================================================
# Two things to fix in the raw PLINK .traw:
#   a) SNP column (col 2): replace positional H3A IDs with EPIC probe names
#      so that cross_study_concordance.py can intersect by name.
#   b) Sample columns (col 7+): PLINK writes "FID_IID" in the header.
#      H3Aflow sets FID == IID == Sample_ID (e.g. "G007"), so the header
#      has "G007_G007".  Strip to plain IID so the sample ID matches the
#      id_map CSV.

cat("\n--- STEP 8: Relabelling SNP IDs and cleaning sample IDs ---\n")

id_remap   <- setNames(all_matches$EPIC_PROBE, all_matches$H3A_SNP)

h3a_lines  <- readLines(traw_path)
h3a_header <- h3a_lines[1]
h3a_data   <- h3a_lines[-1]
h3a_data   <- h3a_data[nchar(trimws(h3a_data)) > 0]

# (a) Relabel SNP IDs
relabelled_data <- sapply(h3a_data, function(ln) {
  parts      <- strsplit(ln, "\t")[[1]]
  h3a_id     <- parts[2]
  if (h3a_id %in% names(id_remap)) parts[2] <- id_remap[[h3a_id]]
  paste(parts, collapse = "\t")
}, USE.NAMES = FALSE)

# (b) Clean sample IDs in the header
# H3Aflow samplesheet: Sample_ID = G007, G009, ...
# PLINK FID_IID encoding: "G007_G007"
# Strategy: for "FID_IID" tokens, if FID == IID return IID alone;
#           if they differ, return IID (the biologically meaningful ID).
header_parts   <- strsplit(h3a_header, "\t")[[1]]
fixed_cols     <- header_parts[1:6]
raw_sample_cols <- header_parts[7:length(header_parts)]

clean_sample_id <- function(fid_iid) {
  # Split on first underscore only
  parts <- strsplit(fid_iid, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(fid_iid)   # no underscore → return as-is
  fid <- parts[1]
  iid <- paste(parts[-1], collapse = "_")  # IID may itself contain underscores
  # If FID == IID (the common H3Aflow convention), return plain IID
  if (fid == iid) return(iid)
  # Otherwise return IID only — FID is typically a family/batch placeholder
  return(iid)
}

cleaned_samples <- sapply(raw_sample_cols, clean_sample_id, USE.NAMES = FALSE)
cat("Sample IDs cleaned (examples):\n")
cat("  Before:", paste(head(raw_sample_cols, 4), collapse = ", "), "\n")
cat("  After :", paste(head(cleaned_samples,  4), collapse = ", "), "\n")

final_header  <- paste(c(fixed_cols, cleaned_samples), collapse = "\t")
final_traw_out <- paste0(opt$out_prefix, "_h3a_snp_calls.traw")
writeLines(c(final_header, relabelled_data), final_traw_out)

cat("✓ Relabelled .traw written:", final_traw_out, "\n")
cat("  SNPs   :", length(relabelled_data), "\n")
cat("  Samples:", length(cleaned_samples), "\n")

# ============================================================
# STEP 9: WRITE MATCH REPORT
# ============================================================
report_out <- paste0(opt$out_prefix, "_snp_match_report.txt")
report_lines <- c(
  "=============================================================",
  "H3A SNP MATCHING REPORT",
  "=============================================================",
  paste0("  H3A PLINK bfile    : ", opt$bfile),
  paste0("  EPIC .traw input   : ", opt$epic_traw),
  paste0("  Array featureset   : ", opt$featureset),
  "",
  "  NOTE: H3Africa .bim uses positional IDs (seq-h3a_37_CHR_POS_REF_ALT),",
  "  not rsIDs.  Matching is done primarily by chr:pos lookup against the",
  "  meffil manifest (GRCh37/hg19 coordinates required).",
  "",
  paste0("  EPIC probes sought : ", length(epic_probes)),
  paste0("  H3A variants total : ", nrow(bim)),
  "",
  paste0("  Matched (total)    : ", nrow(all_matches)),
  paste0("    Positional       : ", sum(all_matches$MATCH_TYPE == "positional")),
  paste0("    Direct name      : ", sum(all_matches$MATCH_TYPE == "direct_name")),
  paste0("  Unmatched          : ", length(still_unmatched_final)),
  ""
)

if (length(still_unmatched_final) > 0) {
  report_lines <- c(report_lines,
    "  Unmatched EPIC probes:",
    paste0("    ", still_unmatched_final),
    ""
  )
}

if (nrow(pos_final) > 0) {
  report_lines <- c(report_lines,
    "  Positionally matched probes (H3A positional ID -> EPIC probe name):",
    apply(pos_final, 1, function(r)
      sprintf("    %-20s <- %s  (chr%s:%s)",
              r["EPIC_PROBE"], r["H3A_SNP"], r["CHR"], r["POS"])),
    ""
  )
}

report_lines <- c(report_lines,
  "=============================================================",
  paste0("  Output .traw       : ", final_traw_out),
  paste0("  Match table CSV    : ", match_table_out),
  "============================================================="
)

writeLines(report_lines, report_out)
cat("✓ Match report written:", report_out, "\n")

cat("\n", rep("=", 60), "\n", sep = "")
cat("H3A SNP MATCHING COMPLETED\n")
cat(rep("=", 60), "\n")
cat("  Matched probes :", nrow(all_matches), "/", length(epic_probes), "\n")
cat("  H3A samples    :", length(cleaned_samples), "\n")
cat("  Output .traw   :", final_traw_out, "\n")
cat(rep("=", 60), "\n")
