#!/usr/bin/env nextflow

/*
 * Meffil-based Methylation QC and Normalization Pipeline
 * Processes samples by existing plate folders
 * Based on meffil best practices: https://github.com/perishky/meffil/wiki/Full-pipeline-for-analysing-massive-datasets
 */

params.idat_dir = "./idats"
params.out_dir = "./results"
params.qc_thresh = 0.95
params.conda_env = "env_meffil_epicv2.yml"
params.rscript = 'Rscript'

// Resolve IDAT directory to absolute path
def IDAT_DIR = params.idat_dir
if (!IDAT_DIR.startsWith('/')) {
    IDAT_DIR = "${projectDir}/${IDAT_DIR}"
}

log.info """
         =========================================
         Meffil Methylation Pipeline (Plate-based)
         =========================================
         IDAT directory : ${IDAT_DIR}
         Output directory : ${params.out_dir}
         QC threshold : ${params.qc_thresh}
         Normalization : Combined (all plates together)
         =========================================
         """

/*
 * PROCESS 1: Group IDATs by existing plate folders
 */
process GROUP_BY_PLATE {
    publishDir "${params.out_dir}/plate_manifests", mode: 'copy'
    conda "${params.conda_env}"
    
    output:
    path "plate_manifests/*.txt", emit: plate_manifests
    path "grouping_summary.txt", emit: summary
    
    script:
    """
    mkdir -p plate_manifests
    ${params.rscript} ${projectDir}/group_by_plate.r \
        ${IDAT_DIR} \
        plate_manifests \
        grouping_summary.txt
    """
}

/*
 * PROCESS 2: Create samplesheets for each plate
 */
process CREATE_SAMPLESHEETS {
    publishDir "${params.out_dir}/samplesheets", mode: 'copy'
    conda "${params.conda_env}"
    
    input:
    path plate_manifest
    
    output:
    tuple path(plate_manifest), path("samplesheets/*.csv"), emit: plate_with_samplesheet
    
    script:
    def plate_id = plate_manifest.baseName
    """
    mkdir -p samplesheets
    ${params.rscript} ${projectDir}/create_samplesheet.r \
        ${plate_manifest} \
        samplesheets/${plate_id}_samplesheet.csv
    """
}

/*
 * PROCESS 3: QC per plate following meffil best practices
 */
process PLATE_QC {
    publishDir "${params.out_dir}/qc_results", mode: 'copy'
    conda "${params.conda_env}"
    
    input:
    tuple path(plate_manifest), path(samplesheet)
    
    output:
    path "qc_results/*", emit: qc_all_files
    path "qc_results/*_qc_objects.rds", emit: qc_objects
    path "qc_results/*_passed_samples.txt", emit: passed_samples
    path "qc_results/*_qc_metrics.csv", emit: qc_metrics
    tuple path(plate_manifest), path(samplesheet), emit: plate_samplesheet
    
    script:
    def plate_id = samplesheet.baseName.replace('_samplesheet', '')
    """
    mkdir -p qc_results
    ${params.rscript} ${projectDir}/plate_qc_meffil.r \
        ${samplesheet} \
        qc_results \
        ${params.qc_thresh}
    """
}

/*
 * PROCESS 3B: Merge samplesheets for combined normalization
 */
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
    
    # Read all samplesheets
    samplesheet_files <- list.files(".", pattern = "_samplesheet\\\\.csv\$", full.names = TRUE)
    
    cat("Found", length(samplesheet_files), "samplesheet(s) to merge:\\n")
    for (f in samplesheet_files) {
      cat("  -", basename(f), "\\n")
    }
    
    if (length(samplesheet_files) == 0) {
      stop("No samplesheets found to merge!")
    }
    
    all_samples <- list()
    for (file in samplesheet_files) {
      ss <- read.csv(file, stringsAsFactors = FALSE)
      cat("  ", basename(file), ":", nrow(ss), "samples\\n")
      all_samples[[length(all_samples) + 1]] <- ss
    }
    
    # Combine all samplesheets
    combined <- do.call(rbind, all_samples)
    
    cat("\\nTotal samples in combined samplesheet:", nrow(combined), "\\n")
    
    # Save combined samplesheet
    write.csv(combined, "combined_samplesheet.csv", row.names = FALSE, quote = TRUE)
    
    cat("Combined samplesheet saved\\n\\n")
    """
}

/*
 * PROCESS 4: Combined normalization (all plates together)
 */
process COMBINED_NORMALIZE {
    publishDir "${params.out_dir}/normalized_combined", mode: 'copy'
    conda "${params.conda_env}"
    
    input:
    path qc_objects_file
    path passed_samples_file
    path combined_samplesheet
    
    output:
    path "normalized_combined/*", emit: normalized_data
    
    script:
    """
    mkdir -p qc_data passed_data normalized_combined
    
    # Organize input files
    mv ${qc_objects_file} qc_data/
    mv ${passed_samples_file} passed_data/
    
    ${params.rscript} ${projectDir}/combined_normalize_meffil.r \
        qc_data \
        passed_data \
        ${combined_samplesheet} \
        normalized_combined \
        ${params.qc_thresh}
    """
}

/*
 * Workflow definition
 */
workflow {
    // Step 1: Group IDATs by existing plate folders
    GROUP_BY_PLATE()
    
    // Step 2: Create samplesheets for each plate
    plate_manifests_ch = GROUP_BY_PLATE.out.plate_manifests.flatten()
    CREATE_SAMPLESHEETS(plate_manifests_ch)
    
    // Step 3: QC each plate independently
    PLATE_QC(CREATE_SAMPLESHEETS.out.plate_with_samplesheet)
    
    // Step 3B: Merge all samplesheets for combined normalization
    all_samplesheets = PLATE_QC.out.plate_samplesheet
        .map { manifest, samplesheet -> samplesheet }
        .collect()
    
    MERGE_SAMPLESHEETS(all_samplesheets)
    
    // Step 4: Combined normalization (all plates together)
    COMBINED_NORMALIZE(
        PLATE_QC.out.qc_objects.collect(),
        PLATE_QC.out.passed_samples.collect(),
        MERGE_SAMPLESHEETS.out.combined_samplesheet
    )
}

workflow.onComplete {
    log.info """
             =========================================
             PIPELINE COMPLETED!
             =========================================
             Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
             Duration: ${workflow.duration}
             Results: ${params.out_dir}
             
             Key outputs:
             - Individual plate QC reports: qc_results/
             - Combined normalized data: normalized_combined/
               * BetaValues_all_plates_combined.csv
               * MValues_all_plates_combined.csv
               * by_plate/ (individual plate subsets)
               * PCA_combined_normalization.png
               * combined_normalization_report.html
             =========================================
             """
}
