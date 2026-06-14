#!/usr/bin/env Rscript
# patch_male_controls.R — final version
#
# Fixes meffil mclapply bug: all male samples receive identical $controls
# vectors, causing meffil.control.matrix(normalize=TRUE) to crash with
# "All control matrix variables have zero variance" on the male subset
# inside meffil.normalize.quantiles().

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript patch_male_controls.R <work_qc_data_dir>")

qc_data_dir <- args[1]
suppressPackageStartupMessages(library(meffil))

cat("\n", rep("=", 60), "\n", sep = "")
cat("PATCH: Re-extracting male control probe values (EPICv2)\n")
cat(rep("=", 60), "\n\n")
cat("QC data dir:", qc_data_dir, "\n\n")

# ── Load all QC objects ───────────────────────────────────────────────────
qc_files <- list.files(qc_data_dir, pattern = "_qc_objects\\.rds$", full.names = TRUE)
if (length(qc_files) == 0) stop("No *_qc_objects.rds files found in: ", qc_data_dir)

cat("Plates found:", length(qc_files), "\n")
all_plates <- list()
for (f in qc_files) {
  pid <- sub("_qc_objects\\.rds$", "", basename(f))
  cat("  Loading", pid, "...")
  obj <- readRDS(f)
  cat(length(obj), "samples\n")
  all_plates[[pid]] <- list(file = f, objects = obj)
}
all_qc   <- do.call(c, lapply(all_plates, function(x) x$objects))
cat("\nTotal QC objects:", length(all_qc), "\n")

male_idx <- which(sapply(all_qc, function(x) x$predicted.sex) == "M")
cat("Male samples:    ", length(male_idx), "\n")
if (length(male_idx) == 0) { cat("No males — nothing to patch.\n"); quit(status=0) }

# ── Confirm corruption ────────────────────────────────────────────────────
ctrl_mat    <- t(sapply(all_qc[male_idx], function(x) x$controls))
unique_rows <- nrow(unique(ctrl_mat))
cat("Unique control rows among males:", unique_rows, "\n")
if (unique_rows > 1) { cat("Controls already vary — no patch needed.\n"); quit(status=0) }
cat("Confirmed: all males share one identical controls vector. Patching...\n\n")

featureset <- all_qc[[male_idx[1]]]$featureset
chip       <- all_qc[[male_idx[1]]]$chip
cat("Featureset:", featureset, "  Chip:", chip, "\n\n")

# ── Reference names from a good female ───────────────────────────────────
female_idx <- which(sapply(all_qc, function(x) x$predicted.sex) == "F")
ref_names  <- names(all_qc[[female_idx[1]]]$controls)
cat("Expected", length(ref_names), "control values\n\n")

# ── Aggregation: lc$values + lc$probes -> 42-element named vector ─────────
# lc$probes has columns: type, target, dye, name  (no address column)
# Addresses ARE in the ref_names as "prefix.ADDRESS" — we use ref_names
# directly to assign names, matched by position within each probe group.
# Groups are ordered identically in lc$probes and ref_names because both
# come from the same meffil.probe.info() table sorted the same way.

make_controls_vector <- function(lc, ref_names) {
  vals   <- lc$values[, 1]
  probes <- lc$probes

  get <- function(tgt, dy, nm=NULL, exact=FALSE) {
    idx <- probes$target == tgt & probes$dye == dy
    if (!is.null(nm)) {
      if (exact) idx <- idx & probes$name %in% nm
      else       idx <- idx & grepl(nm, probes$name, fixed=FALSE)
    }
    vals[idx]
  }
  smean <- function(x) if (length(x)==0) NA_real_ else mean(x, na.rm=TRUE)

  # bisulfite1: mean of paired C-probe sums (C1-C5 present in both channels)
  bs1_G <- get("BISULFITE CONVERSION I", "G", "C[1-5]$")
  bs1_R <- get("BISULFITE CONVERSION I", "R", "C[1-5]$")
  bs1   <- smean(bs1_G + bs1_R)

  bs2   <- smean(get("BISULFITE CONVERSION II", "R"))

  # For address-named elements, extract values in order and name from ref_names
  # ref_names positions (0-indexed from bisulfite1=0):
  #  2-3 : extension.G.*  (2 values)
  #  4-5 : extension.R.*  (2 values)
  #  6-8 : hybe.*         (3 values)
  #  9   : stain.G
  # 10   : stain.R
  # 11-12: nonpoly.G.*    (2 values)
  # 13-14: nonpoly.R.*    (2 values)
  # 15-16: targetrem.*    (2 values)
  # 17-19: spec1.G.*      (3 values)
  # 20-22: spec1.R.*      (3 values)
  # 23-25: spec2.G.*      (3 values)
  # 26-28: spec2.R.*      (3 values)

  ext_G  <- get("EXTENSION",        "G", c("Extension (C)","Extension (G)"), exact=TRUE)
  ext_R  <- get("EXTENSION",        "R", c("Extension (A)","Extension (T)"), exact=TRUE)
  hybe   <- get("HYBRIDIZATION",    "G")
  stainG <- smean(get("STAINING",   "G", "Biotin \\(High\\)"))
  stainR <- smean(get("STAINING",   "R", "DNP \\(High\\)"))
  npG    <- get("NON-POLYMORPHIC",  "G", c("NP (C)","NP (G)"), exact=TRUE)
  npR    <- get("NON-POLYMORPHIC",  "R", c("NP (A)","NP (T)"), exact=TRUE)
  trem   <- get("TARGET REMOVAL",   "G")
  sp1G   <- get("SPECIFICITY I",    "G", sprintf("GT Mismatch %d (PM)", 1:3), exact=TRUE)
  sp1R   <- get("SPECIFICITY I",    "R", sprintf("GT Mismatch %d (PM)", 4:6), exact=TRUE)
  sp2G   <- get("SPECIFICITY II",   "G")
  sp2R   <- get("SPECIFICITY II",   "R")

  # cross-channel values for ratios
  sp1_Rp     <- get("SPECIFICITY I", "R", sprintf("GT Mismatch %d (PM)", 1:3), exact=TRUE)
  sp1_Gp     <- get("SPECIFICITY I", "G", sprintf("GT Mismatch %d (PM)", 4:6), exact=TRUE)
  sp1_ratio1 <- smean(sp1_Rp) / smean(sp1G)
  sp1_ratio2 <- smean(sp1_Gp) / smean(sp1R)
  sp1_ratio  <- (sp1_ratio1 + sp1_ratio2) / 2
  sp2_ratio  <- smean(sp2G)   / smean(sp2R)

  normA  <- smean(get("NORM_A", "R"))
  normT  <- smean(get("NORM_T", "R"))
  normC  <- smean(get("NORM_C", "G"))
  normG  <- smean(get("NORM_G", "G"))
  dye_b  <- (normC + normG) / (normA + normT)

  # OOB: no OOB probes in EPICv2 — use NEGATIVE as surrogate
  oob_G  <- get("NEGATIVE", "G")
  oob_R  <- get("NEGATIVE", "R")
  safe_q <- function(x, p) {
    if (length(x)==0) setNames(rep(NA_real_,length(p)), paste0(p*100,"%"))
    else quantile(x, p, na.rm=TRUE)
  }
  oq_G   <- safe_q(oob_G, c(0.01,0.50,0.99))
  oq_R   <- safe_q(oob_R, c(0.01,0.50,0.99))
  if (!is.na(oq_R["50%"]) && oq_R["50%"] < 1) oq_R["50%"] <- 1
  oob_ratio <- oq_G["50%"] / oq_R["50%"]

  # Assemble in ref_names order, naming address-based entries from ref_names
  result <- c(
    bisulfite1=bs1, bisulfite2=bs2,
    setNames(ext_G,  ref_names[3:4]),
    setNames(ext_R,  ref_names[5:6]),
    setNames(hybe,   ref_names[7:9]),
    stain.G=stainG,  stain.R=stainR,
    setNames(npG,    ref_names[12:13]),
    setNames(npR,    ref_names[14:15]),
    setNames(trem,   ref_names[16:17]),
    setNames(sp1G,   ref_names[18:20]),
    setNames(sp1R,   ref_names[21:23]),
    setNames(sp2G,   ref_names[24:26]),
    setNames(sp2R,   ref_names[27:29]),
    spec1.ratio1=sp1_ratio1, spec1.ratio=sp1_ratio,
    spec2.ratio=sp2_ratio,   spec1.ratio2=sp1_ratio2,
    normA=normA, normC=normC, normT=normT, normG=normG,
    dye.bias=dye_b,
    "oob.G.1%"  = unname(oq_G["1%"]),
    "oob.G.50%" = unname(oq_G["50%"]),
    "oob.G.99%" = unname(oq_G["99%"]),
    oob.ratio   = unname(oob_ratio)
  )

  # Sanity check
  if (length(result) != length(ref_names))
    stop(sprintf("Length mismatch: got %d, expected %d", length(result), length(ref_names)))
  if (!all(names(result) == ref_names))
    stop("Name mismatch:\n  got: ", paste(names(result), collapse=","),
         "\n  exp: ", paste(ref_names, collapse=","))
  result
}

# ── Patch each male ───────────────────────────────────────────────────────
cat("Re-extracting controls for", length(male_idx), "males:\n")
n_patched <- 0L
n_failed  <- 0L

for (i in male_idx) {
  obj   <- all_qc[[i]]
  sname <- obj$sample.name
  cat(sprintf("  [%d/%d] %-10s ... ",
              which(male_idx == i), length(male_idx), sname))

  new_ctrl <- tryCatch({
    lc <- meffil:::meffil.load.controls(
            obj$samplesheet, featureset=featureset,
            chip=chip, verbose=FALSE)
    make_controls_vector(lc, ref_names)
  }, error = function(e) { cat("FAILED (", conditionMessage(e), ")\n"); NULL })

  if (is.null(new_ctrl)) { n_failed <- n_failed + 1L; next }

  old_bs1 <- obj$controls[["bisulfite1"]]
  new_bs1 <- new_ctrl[["bisulfite1"]]
  cat(sprintf("bisulfite1: %.1f -> %.1f\n", old_bs1, new_bs1))
  all_qc[[i]]$controls <- new_ctrl
  n_patched <- n_patched + 1L
}

cat(sprintf("\nPatched: %d  Failed: %d\n\n", n_patched, n_failed))
if (n_patched == 0) stop("No samples patched.")

# ── Verify ────────────────────────────────────────────────────────────────
cat("--- Verifying ---\n")
ctrl_new <- t(sapply(all_qc[male_idx], function(x) x$controls))
vars_new <- apply(ctrl_new, 2, var, na.rm=TRUE)
cat("Unique rows:", nrow(unique(ctrl_new)), "\n")
cat("Zero-var cols:", sum(vars_new==0), "/", ncol(ctrl_new), "\n")
cat("bisulfite1 range:", round(min(ctrl_new[,"bisulfite1"]),1),
    "-", round(max(ctrl_new[,"bisulfite1"]),1), "\n\n")
if (sum(vars_new==0) == ncol(ctrl_new))
  stop("Patch did not restore variance.")

# ── Save ──────────────────────────────────────────────────────────────────
cat("--- Saving patched QC objects ---\n")
for (pid in names(all_plates)) {
  f       <- all_plates[[pid]]$file
  new_obj <- all_plates[[pid]]$objects
  for (s in names(new_obj))
    if (s %in% names(all_qc)) new_obj[[s]] <- all_qc[[s]]
  bak <- paste0(f, ".bak")
  if (!file.exists(bak)) { file.copy(f, bak); cat("  Backed up:", basename(bak), "\n") }
  saveRDS(new_obj, f)
  cat("  Saved:", basename(f), "\n")
}

# ── Normalize test ────────────────────────────────────────────────────────
cat("\n--- Testing meffil.normalize.quantiles ---\n")
all_qc_disk <- list()
for (f in list.files(qc_data_dir, pattern="_qc_objects\\.rds$", full.names=TRUE))
  all_qc_disk <- c(all_qc_disk, readRDS(f))

pass_files <- unique(unlist(lapply(
  c(file.path(dirname(qc_data_dir), "passed_data"), dirname(qc_data_dir)),
  function(d) list.files(d, pattern="_passed_samples\\.txt$", full.names=TRUE))))
pass_files <- pass_files[file.exists(pass_files)]

if (length(pass_files) > 0) {
  all_passed <- unlist(lapply(pass_files, function(f) {
    s <- readLines(f, warn=FALSE); s[nzchar(trimws(s))] }))
  qc_passed <- all_qc_disk[names(all_qc_disk) %in% all_passed]
  cat("Passed samples:", length(qc_passed), "\n")
  norm_test <- tryCatch(
    meffil.normalize.quantiles(qc_passed, number.pcs=16, verbose=FALSE),
    error=function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(norm_test))
    cat("✓ meffil.normalize.quantiles SUCCEEDED -", length(norm_test), "objects\n")
} else {
  cat("NOTE: passed_samples files not found — run -resume to confirm.\n")
}

cat("\n", rep("=", 60), "\n", sep="")
cat("PATCH COMPLETE — patched", n_patched, "of", length(male_idx), "males\n")
cat(rep("=", 60), "\n\n")
cat("Next step:\n")
cat("  nextflow run main_meffil.nf -profile local \\\n")
cat("    --idat_dir ./idats --out_dir ./results \\\n")
cat("    --samplesheet ./samplesheets/all_6_epic.csv \\\n")
cat("    -resume\n\n")
