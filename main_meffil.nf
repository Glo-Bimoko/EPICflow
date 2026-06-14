#!/usr/bin/env nextflow

/*
 * Meffil-based Methylation QC and Normalization Pipeline
 * Processes samples by existing plate folders
 * Based on meffil best practices:
 *   https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets
 *
 * Steps:
 *   1  SPLIT_SAMPLESHEET      – split study CSV into per-plate meffil samplesheets
 *   2  PLATE_QC               – meffil QC per plate (produces QC objects + passed lists)
 *   2B MERGE_SAMPLESHEETS     – combine per-plate samplesheets for joint normalization
 *   4  COMBINED_NORMALIZE     – functional normalization across all plates together
 *   5  GENOTYPE_CONCORDANCE   – all-vs-all identity check on the 65 rs-probe SNPs
 */

params.idat_dir              = "./idats"
params.out_dir               = "./results"
params.qc_thresh             = 0.95
params.conda_env             = "env_meffil_epicv2.yml"
params.rscript               = 'Rscript'
params.conc_threshold        = 0.99   // flag pairs with concordance >= this as duplicates
params.conc_min_snps         = 50     // min jointly-called SNPs for a reliable comparison
params.plate_qc_mem          = 28     // GB RAM reserved for each PLATE_QC process
                                       // (local: keep below machine total; HPC: set to process allocation)
params.norm_mem               = 28     // GB RAM reserved for COMBINED_NORMALIZE
params.norm_n_pcs             = 0      // control-probe PCs for FN: 0 = auto (cross-validated scree)

// ── Study samplesheet ──────────────────────────────────────────────────────
// Required CSV mapping every array position to a study participant.
// The pipeline reads this file, splits it by Plate Number, and feeds
// one samplesheet per plate directly into PLATE_QC — no intermediate
// grouping or samplesheet-creation step needed.
//
// Required columns (case-insensitive, spaces or underscores):
//   Sample ID | BeadChip Barcode | Sentrix Position | Plate Number
// Optional:
//   Collected Gender / Sex / Gender   → used for meffil sex prediction check
//
// --samplesheet and --sample_map are accepted synonyms.
params.samplesheet           = ""     // e.g. "./samplesheets/all_6_epic.csv"
params.sample_map            = ""     // synonym for --samplesheet

// ── H3Africa cross-study concordance ──────────────────────────────────────
// Set --h3a_bfile to enable cross-study ID checking (Steps 5-pre and 6).
// Provide the PLINK bfile prefix (no extension) for the H3Africa genotyping
// dataset produced by H3Aflow, e.g. "--h3a_bfile /data/h3a/h3africa_gwas"
// The dataset must be in bed/bim/fam format with variant IDs in column 2 of
// the .bim file.  rsIDs present in both the EPIC array and the H3A dataset
// are used directly; remaining probes are matched by chromosomal position.
params.h3a_bfile             = ""     // e.g. "/data/h3a/h3africa_gwas"; empty = skip
params.h3a_conc_threshold    = 0.80   // score threshold for reporting (not classification)
params.h3a_conc_min_snps     = 20     // min SNPs for a reliable cross-study pair
params.h3a_featureset        = "epicv2"  // meffil featureset for positional SNP lookup
params.plink                 = "plink"  // path to PLINK 1.9 binary
// Path to snp-names.txt — authoritative list of rs-probe SNPs used for H3A
// extraction.  When provided, ALL 65 SNPs are sought in the H3A .bim even if
// some were dropped during EPIC QC (e.g. low call rate).  Leave empty to fall
// back to probes retained in the EPIC .traw only.
params.snp_names_file        = "${projectDir}/bin/snp-names.txt"
// Path to identity map TSV: epic_sample_id<TAB>h3a_sample_id (one row per
// participant enrolled in both studies).  Enables assigned-vs-best
// classification (correct / swap / outside_dataset / low_confidence).
// Leave empty to run best-match-only mode (no classification).
params.id_map                = ""


// Resolve --samplesheet alias: if --sample_map is not set but --samplesheet
// is, treat them as equivalent so common invocation typos don't silently skip
// Resolve --samplesheet / --sample_map synonym.
def STUDY_SS = params.samplesheet ?: params.sample_map ?: ""
if (!STUDY_SS) { error "Please supply --samplesheet <path/to/all_6_epic.csv>" }

// Resolve IDAT directory to absolute path
def IDAT_DIR = params.idat_dir
if (!IDAT_DIR.startsWith('/')) {
    IDAT_DIR = "${projectDir}/${IDAT_DIR}"
}

log.info """
         =========================================
         Meffil Methylation Pipeline (Plate-based)
         =========================================
         IDAT directory      : ${IDAT_DIR}
         Output directory    : ${params.out_dir}
         QC threshold        : ${params.qc_thresh}
         Study samplesheet   : ${STUDY_SS}
         Normalization       : Combined (all plates together)
         Concordance cutoff  : ${params.conc_threshold}
         Min SNPs (reliable) : ${params.conc_min_snps}
         H3A bfile           : ${params.h3a_bfile ?: '(not provided — cross-study step skipped)'}
         H3A match threshold : ${params.h3a_conc_threshold}
         SNP names file      : ${params.snp_names_file ?: '(not set — using EPIC .traw probes only)'}
         ID map              : ${params.id_map ?: '(not set — best-match only, no classification)'}
         =========================================
         """

// ============================================================
// PROCESS 1: Split study samplesheet into per-plate meffil samplesheets
// ============================================================
// Reads the study-level CSV (--samplesheet), splits by Plate Number, and
// writes one meffil-ready CSV per plate to the work directory.
// Each emitted file is then processed independently by PLATE_QC.
// Rows with a missing Plate Number are written to an "unassigned" file
// and reported as a warning; they are excluded from QC and normalization.
process SPLIT_SAMPLESHEET {
    publishDir "${params.out_dir}/samplesheets", mode: 'copy'
    conda "${params.conda_env}"

    output:
    path "samplesheets/*_samplesheet.csv", emit: per_plate_samplesheets

    script:
    def ss_path = file(STUDY_SS).toAbsolutePath()
    """
    mkdir -p samplesheets
    ${params.rscript} ${projectDir}/bin/split_samplesheet.r \
        ${ss_path} \
        ${IDAT_DIR} \
        samplesheets
    """
}

// ============================================================
// PROCESS 2: QC per plate following meffil best practices
// ============================================================
process PLATE_QC {
    publishDir "${params.out_dir}/qc_results", mode: 'copy'
    memory "${params.plate_qc_mem} GB"
    conda "${params.conda_env}"

    input:
    path samplesheet

    output:
    path "qc_results/*",                    emit: qc_all_files
    path "qc_results/*_qc_objects.rds",     emit: qc_objects
    path "qc_results/*_passed_samples.txt", emit: passed_samples
    path "qc_results/*_qc_summary.rds",     emit: qc_summaries
    path "qc_results/*_qc_metrics.csv",     emit: qc_metrics
    path samplesheet,                       emit: samplesheet_out

    script:
    """
    mkdir -p qc_results
    ${params.rscript} ${projectDir}/bin/plate_qc_meffil.r \
        ${samplesheet} \
        qc_results \
        ${params.qc_thresh} \
        ${params.plate_qc_mem}
    """
}

// ============================================================
// PROCESS 2B: Merge samplesheets for combined normalization
// ============================================================
process MERGE_SAMPLESHEETS {
    publishDir "${params.out_dir}/samplesheets", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path samplesheet_files

    output:
    path "combined_samplesheet.csv", emit: combined_samplesheet

    script:
    """
    #!/usr/bin/env Rscript

    cat("\\n=== MERGING SAMPLESHEETS ===\\n")

    samplesheet_files <- list.files(".", pattern = "_samplesheet\\\\.csv\$",
                                    full.names = TRUE)

    cat("Found", length(samplesheet_files), "samplesheet(s) to merge:\\n")
    for (f in samplesheet_files) cat("  -", basename(f), "\\n")

    if (length(samplesheet_files) == 0) {
        stop("No samplesheets found to merge!")
    }

    all_samples <- lapply(samplesheet_files, function(f) {
        ss <- read.csv(f, stringsAsFactors = FALSE)
        cat(" ", basename(f), ":", nrow(ss), "samples\\n")
        ss
    })

    combined <- do.call(rbind, all_samples)
    cat("\\nTotal samples in combined samplesheet:", nrow(combined), "\\n")

    write.csv(combined, "combined_samplesheet.csv",
              row.names = FALSE, quote = TRUE)
    cat("Combined samplesheet saved\\n")
    """
}

// ============================================================
// PROCESS 4: Combined normalization (all plates together)
// ============================================================
process COMBINED_NORMALIZE {
    publishDir "${params.out_dir}/normalized_combined", mode: 'copy'
    memory "${params.norm_mem} GB"
    conda "${params.conda_env}"

    input:
    path qc_objects_file
    path passed_samples_file
    path qc_summary_file
    path combined_samplesheet

    output:
    path "normalized_combined/*",                  emit: normalized_data
    path "normalized_combined/beta_normalised.gds", emit: gds_file
    // Pass QC objects directory forward for concordance
    path "qc_data",                                emit: qc_objects_dir
    path "passed_data",                            emit: passed_samples_dir

    script:
    """
    mkdir -p qc_data passed_data qc_summaries normalized_combined

    mv ${qc_objects_file}   qc_data/
    mv ${passed_samples_file} passed_data/
    mv ${qc_summary_file}   qc_summaries/

    ${params.rscript} ${projectDir}/bin/combined_normalize_meffil.r \
        qc_data \
        passed_data \
        qc_summaries \
        ${combined_samplesheet} \
        normalized_combined \
        ${params.qc_thresh} \
        ${params.norm_mem} \
        ${params.norm_n_pcs}
    """
}


// ============================================================
// PROCESS 5-PRE: Pre-QC identity check against H3Africa genotype data
// ============================================================
// Runs immediately after PLATE_QC — before any sample is excluded.
// Uses ALL QC objects (passed and failed), matching the approach from
// Sinenhlanhla Mthembu / PURE-SA-NW (methylation_genotype_matching.R):
// a swapped or contaminated sample may fail QC *because* of the swap,
// so identity checking must happen before exclusion.
//
// Only runs when --h3a_bfile is provided.
process IDENTITY_CHECK_PRE_QC {
    publishDir "${params.out_dir}/cross_study_pre_qc", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path qc_objects_files     // all plate *_qc_objects.rds (before filtering)
    path passed_samples_files // passed lists (present but not used in --pre_qc mode)

    output:
    path "pre_qc_check/*",                               emit: all_outputs
    path "pre_qc_check/*_concordance_report.html",       emit: report

    script:
    """
    mkdir -p pre_qc_check qc_data passed_data

    # Stage inputs (collect emits multiple files; stage them into subdirs)
    for f in ${qc_objects_files}; do cp \$f qc_data/; done
    for f in ${passed_samples_files}; do cp \$f passed_data/; done

    # ── Step A: extract SNP betas from ALL samples (--pre_qc skips the
    #            passed-samples filter, so QC failures are included) ────────
    ${params.rscript} ${projectDir}/bin/extract_snp_betas.r \
        qc_data \
        passed_data \
        pre_qc_check/pre_qc_snps \
        --pre_qc

    # ── Step B: map EPIC probes → H3A variants ─────────────────────────────
    ${params.rscript} ${projectDir}/bin/match_h3a_snps.r \
        --bfile      ${params.h3a_bfile} \
        --epic_traw  pre_qc_check/pre_qc_snps_snp_calls.traw \
        --out_prefix pre_qc_check/h3a_snps \
        --featureset ${params.h3a_featureset} \
        --plink      ${params.plink} \
        ${params.snp_names_file ? "--snp_names  ${params.snp_names_file}" : ""}

    # ── Step C: cross-study identity check ────────────────────────────────
    python3 ${projectDir}/bin/cross_study_concordance.py \
        --epic_traw   pre_qc_check/pre_qc_snps_snp_calls.traw \
        --h3a_traw    pre_qc_check/h3a_snps_h3a_snp_calls.traw \
        --out         pre_qc_check \
        --prefix      pre_qc \
        --threshold   ${params.h3a_conc_threshold} \
        --min_snps    ${params.h3a_conc_min_snps} \
        ${params.id_map ? "--id_map  ${params.id_map}" : ""}
    """
}

// ============================================================
// PROCESS 5: Genotype concordance on EPIC rs-probe SNPs
// ============================================================
process GENOTYPE_CONCORDANCE {
    publishDir "${params.out_dir}/concordance", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path qc_objects_dir
    path passed_samples_dir

    output:
    path "concordance/*",                            emit: all_outputs
    path "concordance/*_flagged_duplicates.tsv",     emit: flagged
    path "concordance/*_concordance_report.html",    emit: report

    script:
    """
    mkdir -p concordance

    # ── Step A: extract rs-probe betas → genotype calls → .traw ──────────
    ${params.rscript} ${projectDir}/bin/extract_snp_betas.r \
        ${qc_objects_dir} \
        ${passed_samples_dir} \
        concordance/methylation_snps

    # ── Step B: all-vs-all pairwise concordance ───────────────────────────
    python3 ${projectDir}/bin/pairwise_concordance.py \
        --traw      concordance/methylation_snps_snp_calls.traw \
        --out       concordance \
        --prefix    methylation_snps \
        --threshold ${params.conc_threshold} \
        --min_snps  ${params.conc_min_snps}

    # ── Step C: count SNPs used (for report) ─────────────────────────────
    n_snps=\$(python3 -c "
import csv
with open('concordance/methylation_snps_snp_calls.traw') as f:
    rows = list(csv.reader(f, delimiter='\\t'))
print(len(rows) - 1)
" 2>/dev/null || echo 65)

    # ── Step D: HTML report ───────────────────────────────────────────────
    python3 ${projectDir}/bin/concordance_report.py \
        --matrix    concordance/methylation_snps_concordance_matrix.csv \
        --pairs     concordance/methylation_snps_concordance_pairs.tsv \
        --flagged   concordance/methylation_snps_flagged_duplicates.tsv \
        --summary   concordance/methylation_snps_concordance_summary.txt \
        --out       concordance/methylation_snps_concordance_report.html \
        --threshold ${params.conc_threshold} \
        --n_snps    \${n_snps}
    """
}

// ============================================================
// PROCESS 5B: Methylation profile concordance for flagged pairs
// ============================================================
// Runs after GENOTYPE_CONCORDANCE.  For every pair flagged as a likely
// duplicate (concordance ≥ params.conc_cutoff on the 65 rs-probe SNPs),
// this process loads their genome-wide beta values from the normalised
// GDS file and computes three complementary similarity metrics:
//
//   Pearson r        – correlation of genome-wide beta profiles
//   Median |Δβ|      – median absolute beta difference across all CpGs
//   Frac |Δβ|>0.2    – fraction of CpGs with a large discordant swing
//
// Verdicts per pair:
//   CONFIRMED_DUPLICATE   – all three methylation thresholds met;
//                           these are the same biological sample
//   BORDERLINE            – close but not all criteria satisfied
//   METHYLATION_MISMATCH  – same genotype, divergent methylation;
//                           consistent with a sample swap or MZ twins
//   INSUFFICIENT_DATA     – too few CpGs available to decide
//
// The process exits cleanly with empty placeholder outputs when
// flagged_duplicates.tsv contains no data rows (zero duplicates),
// so it never blocks a clean run.
//
// Parameters (set in nextflow.config or pass on CLI):
//   params.meth_r_threshold      – Pearson r must be ≥ this  (default 0.999)
//   params.meth_delta_threshold  – Median |Δβ| must be ≤ this (default 0.01)
//   params.meth_large_delta_frac – Frac |Δβ|>0.2 must be ≤ this (default 0.01)
//   params.meth_chunk_size       – CpGs per GDS read chunk    (default 5000)
//
// Outputs published to: results/methylation_concordance/
process METHYLATION_CONCORDANCE {
    publishDir "${params.out_dir}/methylation_concordance", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path flagged_tsv   // *_flagged_duplicates.tsv from GENOTYPE_CONCORDANCE
    path gds_file      // beta_normalised.gds from COMBINED_NORMALIZE

    output:
    path "methdup/*_methylation_duplicate_evidence.tsv", emit: evidence
    path "methdup/*_methylation_duplicate_summary.txt",  emit: summary
    path "methdup/*_methylation_duplicate_report.html",  emit: report

    script:
    def r_thr   = params.meth_r_threshold      ?: 0.999
    def d_thr   = params.meth_delta_threshold   ?: 0.01
    def fld_thr = params.meth_large_delta_frac  ?: 0.01
    def chunk   = params.meth_chunk_size        ?: 5000
    """
    mkdir -p methdup

    # Exit cleanly if no pairs were flagged (header-only or empty file)
    n_pairs=\$(tail -n +2 ${flagged_tsv} | wc -l)
    if [ "\${n_pairs}" -eq 0 ]; then
        echo "No flagged pairs — skipping methylation concordance check." | \\
            tee methdup/methylation_snps_methylation_duplicate_summary.txt
        touch methdup/methylation_snps_methylation_duplicate_evidence.tsv
        touch methdup/methylation_snps_methylation_duplicate_report.html
        exit 0
    fi

    python3 ${projectDir}/bin/methylation_duplicate_check.py \\
        --flagged          ${flagged_tsv} \\
        --gds              ${gds_file} \\
        --out              methdup \\
        --prefix           methylation_snps \\
        --r_threshold      ${r_thr} \\
        --delta_threshold  ${d_thr} \\
        --large_delta_frac ${fld_thr} \\
        --chunk_size       ${chunk}
    """
}

// ============================================================
// PROCESS 6: Cross-study concordance against H3Africa genotype data
// ============================================================
// Only runs when --h3a_bfile is provided.
// Requires PLINK 1.9 on PATH (or specify via --plink).
//
// Inputs:
//   qc_objects_dir  – directory of *_qc_objects.rds (from COMBINED_NORMALIZE)
//   passed_samples_dir – directory of *_passed_samples.txt
//   epic_traw       – methylation-derived SNP calls (from GENOTYPE_CONCORDANCE)
//
// What it does:
//   Step A  match_h3a_snps.r
//           – Reads the 65 rs-probe names from the EPIC .traw
//           – Looks them up in the H3A .bim by rsID (direct match)
//           – Falls back to chr:pos matching for any remaining probes
//           – Calls PLINK to extract those variants and write a .traw
//           – Relabels H3A variant IDs to EPIC probe names for comparison
//
//   Step B  cross_study_concordance.py
//           – Intersects SNPs present in both .traw files
//           – Computes all EPIC × H3A pairwise concordance
//           – Reports pairs above --h3a_conc_threshold as ID-confirmed
//
// Outputs published to: results/cross_study_concordance/
process H3A_CONCORDANCE {
    publishDir "${params.out_dir}/cross_study_concordance", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path epic_traw           // methylation_snps_snp_calls.traw from GENOTYPE_CONCORDANCE

    output:
    path "cross_study/*",                            emit: all_outputs
    path "cross_study/*_per_sample_verdict.tsv",     emit: verdicts
    path "cross_study/*_concordance_report.html",      emit: report

    script:
    """
    mkdir -p cross_study

    # ── Step A: map EPIC probes → H3A variants, extract H3A .traw ────────
    ${params.rscript} ${projectDir}/bin/match_h3a_snps.r \
        --bfile      ${params.h3a_bfile} \
        --epic_traw  ${epic_traw} \
        --out_prefix cross_study/h3a_snps \
        --featureset ${params.h3a_featureset} \
        --plink      ${params.plink} \
        ${params.snp_names_file ? "--snp_names  ${params.snp_names_file}" : ""}

    # ── Step B: all EPIC × all H3A identity check ────────────────────────
    python3 ${projectDir}/bin/cross_study_concordance.py \
        --epic_traw   ${epic_traw} \
        --h3a_traw    cross_study/h3a_snps_h3a_snp_calls.traw \
        --out         cross_study \
        --prefix      cross_study \
        --threshold   ${params.h3a_conc_threshold} \
        --min_snps    ${params.h3a_conc_min_snps} \
        ${params.id_map ? "--id_map  ${params.id_map}" : ""}
    """
}

// ============================================================
// Workflow
// ============================================================
workflow {

    // Step 1 — split study CSV into per-plate samplesheets
    SPLIT_SAMPLESHEET()

    // Step 2 — QC each plate independently (runs in parallel)
    per_plate_ch = SPLIT_SAMPLESHEET.out.per_plate_samplesheets.flatten()
    PLATE_QC(per_plate_ch)

    // Step 2B — merge per-plate samplesheets for combined normalization
    all_samplesheets = PLATE_QC.out.samplesheet_out.collect()
    MERGE_SAMPLESHEETS(all_samplesheets)

    // Step 4
    COMBINED_NORMALIZE(
        PLATE_QC.out.qc_objects.collect(),
        PLATE_QC.out.passed_samples.collect(),
        PLATE_QC.out.qc_summaries.collect(),
        MERGE_SAMPLESHEETS.out.combined_samplesheet
    )

    // Step 5 — genotype concordance
    GENOTYPE_CONCORDANCE(
        COMBINED_NORMALIZE.out.qc_objects_dir,
        COMBINED_NORMALIZE.out.passed_samples_dir
    )

    // Step 5B — methylation profile confirmation for genotype-flagged pairs
    flagged_tsv_ch = GENOTYPE_CONCORDANCE.out.flagged
    METHYLATION_CONCORDANCE(
        flagged_tsv_ch,
        COMBINED_NORMALIZE.out.gds_file
    )

    // Step 5-PRE — pre-QC identity check (all samples, before exclusion)
    // Only runs if --h3a_bfile is provided.
    if (params.h3a_bfile) {
        IDENTITY_CHECK_PRE_QC(
            PLATE_QC.out.qc_objects.collect(),
            PLATE_QC.out.passed_samples.collect()
        )
    } else {
        log.info "Skipping IDENTITY_CHECK_PRE_QC: --h3a_bfile not provided."
    }

    // Step 6 — post-QC cross-study concordance (QC-passed samples only)
    // Only runs if --h3a_bfile is provided.
    if (params.h3a_bfile) {
        epic_traw_ch = GENOTYPE_CONCORDANCE.out.all_outputs
            .flatten()
            .filter { it.name.endsWith("_snp_calls.traw") }

        H3A_CONCORDANCE(epic_traw_ch)
    } else {
        log.info "Skipping H3A_CONCORDANCE: --h3a_bfile not provided."
    }
}

workflow.onComplete {
    log.info """
             =========================================
             PIPELINE COMPLETED
             =========================================
             Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
             Duration : ${workflow.duration}
             Results  : ${params.out_dir}

             Key outputs:
               Per-plate samplesheets     : samplesheets/
               QC reports (per plate)     : qc_results/
               Combined normalized data   : normalized_combined/
                 BetaValues_all_plates_combined.csv
                 MValues_all_plates_combined.csv
                 by_plate/
               Genotype concordance       : concordance/
                 methylation_snps_concordance_report.html  ← review this
                 methylation_snps_flagged_duplicates.tsv   ← within-study duplicates
                 methylation_snps_concordance_matrix.csv
               Methylation concordance    : methylation_concordance/
                 methylation_snps_methylation_duplicate_report.html  ← confirms duplicates
                 methylation_snps_methylation_duplicate_evidence.tsv ← per-pair metrics
                 methylation_snps_methylation_duplicate_summary.txt
               Cross-study (pre-QC)       : cross_study_pre_qc/  (if --h3a_bfile set)
                 pre_qc_concordance_report.html        ← ⭐ review FIRST (all samples)
                 pre_qc_per_sample_verdict.tsv         ← per-sample classification
                 pre_qc_classification_summary.tsv     ← category counts
               Cross-study (post-QC)      : cross_study_concordance/
                 cross_study_concordance_report.html     ← review after pre-QC
                 cross_study_per_sample_verdict.tsv     ← QC-passed samples only
                 h3a_snps_snp_match_report.txt             ← SNP mapping summary
             =========================================
             """
}
