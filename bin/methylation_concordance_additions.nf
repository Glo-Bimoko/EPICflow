// ============================================================
// PROCESS 5B: Methylation profile concordance for flagged pairs
// ============================================================
// Runs only when GENOTYPE_CONCORDANCE flags at least one duplicate pair.
// Takes the flagged_duplicates.tsv + the normalised beta GDS file and
// computes three complementary methylation similarity metrics per pair:
//
//   Pearson r          – near 1.0 for true duplicates
//   Median |Δβ|        – near 0 for true duplicates
//   Frac |Δβ| > 0.2   – near 0 for true duplicates
//
// Verdict per pair:
//   CONFIRMED_DUPLICATE   – all three thresholds satisfied
//   BORDERLINE            – close but not all criteria met
//   METHYLATION_MISMATCH  – same genotype, different methylation profile
//                           → possible sample swap or monozygotic twins
//   INSUFFICIENT_DATA     – too few jointly-called CpGs
//
// Parameters (add to nextflow.config or pass on CLI):
//   params.meth_r_threshold      = 0.999   // Pearson r must exceed this
//   params.meth_delta_threshold  = 0.01    // Median |Δβ| must be below this
//   params.meth_large_delta_frac = 0.01    // Fraction |Δβ|>0.2 must be below this
//   params.meth_chunk_size       = 5000    // CpGs per GDS read chunk (memory tuning)
//
// Outputs published to: results/methylation_concordance/
process METHYLATION_CONCORDANCE {
    publishDir "${params.out_dir}/methylation_concordance", mode: 'copy'
    conda "${params.conda_env}"

    input:
    path flagged_tsv   // methylation_snps_flagged_duplicates.tsv from GENOTYPE_CONCORDANCE
    path gds_file      // beta_normalised.gds from COMBINED_NORMALIZE

    output:
    path "methdup/*_methylation_duplicate_evidence.tsv", emit: evidence
    path "methdup/*_methylation_duplicate_summary.txt",  emit: summary
    path "methdup/*_methylation_duplicate_report.html",  emit: report

    script:
    def r_thr    = params.meth_r_threshold      ?: 0.999
    def d_thr    = params.meth_delta_threshold   ?: 0.01
    def fld_thr  = params.meth_large_delta_frac  ?: 0.01
    def chunk    = params.meth_chunk_size        ?: 5000
    """
    mkdir -p methdup

    # Skip gracefully if no pairs were flagged (empty file or header-only)
    n_pairs=\$(tail -n +2 ${flagged_tsv} | wc -l)
    if [ "\${n_pairs}" -eq 0 ]; then
        echo "No flagged pairs — skipping methylation concordance check." | \\
            tee methdup/methylation_snps_methylation_duplicate_summary.txt
        touch methdup/methylation_snps_methylation_duplicate_evidence.tsv
        touch methdup/methylation_snps_methylation_duplicate_report.html
        exit 0
    fi

    python3 ${projectDir}/bin/methylation_duplicate_check.py \\
        --flagged   ${flagged_tsv} \\
        --gds       ${gds_file} \\
        --out       methdup \\
        --prefix    methylation_snps \\
        --r_threshold      ${r_thr} \\
        --delta_threshold  ${d_thr} \\
        --large_delta_frac ${fld_thr} \\
        --chunk_size       ${chunk}
    """
}

// ============================================================
// WORKFLOW ADDITION
// ============================================================
// Insert this block AFTER the existing GENOTYPE_CONCORDANCE call:
//
//   // Step 5B — methylation profile confirmation of genotype-flagged pairs
//   flagged_tsv_ch = GENOTYPE_CONCORDANCE.out.all_outputs
//       .flatten()
//       .filter { it.name.endsWith("_flagged_duplicates.tsv") }
//
//   METHYLATION_CONCORDANCE(
//       flagged_tsv_ch,
//       COMBINED_NORMALIZE.out.gds_file   // ← you need to expose this emit (see note below)
//   )
//
// ──────────────────────────────────────────────────────────────
// NOTE: COMBINED_NORMALIZE currently does not emit the GDS file as a
// named channel.  Add this output stanza to the COMBINED_NORMALIZE process:
//
//   path "normalized_combined/beta_normalised.gds", emit: gds_file
//
// And ensure combined_normalize_meffil.r writes to the expected relative path
// (it already does: gds_file <- file.path(out_dir, "beta_normalised.gds")).
// ──────────────────────────────────────────────────────────────
//
// Also add these params to nextflow.config or params block in main_meffil.nf:
//
//   params.meth_r_threshold      = 0.999
//   params.meth_delta_threshold  = 0.01
//   params.meth_large_delta_frac = 0.01
//   params.meth_chunk_size       = 5000
