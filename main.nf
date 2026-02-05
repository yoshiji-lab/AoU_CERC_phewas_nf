nextflow.enable.dsl = 2

include { DOWNLOAD_PRS }                 from './modules/download_prs'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COUNTS } from './modules/download_gcs_or_local'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COVARS } from './modules/download_gcs_or_local'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_EXCLUDE } from './modules/download_gcs_or_local'
include { PREPARE_PRS_COHORTS }          from './modules/prepare_prs_cohorts'
include { FILTER_COUNTS_TO_COHORT }      from './modules/filter_counts_to_cohort'
include { RUN_PHEWAS }                   from './modules/run_phewas'
include { CUSTOM_PLOT }                  from './modules/custom_plot'

workflow {
  // ========================================
  // 0) Download/stage inputs
  // ========================================
  ch_counts_in = Channel.value(params.phecode_counts)
  ch_covars_in = Channel.value(params.cohort_covariates)
  ch_prs_in    = Channel.value(params.prs_matrix)
  ch_exclude_in = Channel.value(params.exclude_ids)
 
  
  DOWNLOAD_COUNTS(ch_counts_in)
  DOWNLOAD_COVARS(ch_covars_in)
  DOWNLOAD_PRS(ch_prs_in)
  DOWNLOAD_EXCLUDE(ch_exclude_in)
  
  // Rename files to avoid collision in PREPARE_PRS_COHORTS
  ch_counts  = DOWNLOAD_COUNTS.out.file
  ch_covars  = DOWNLOAD_COVARS.out.file
  ch_prs     = DOWNLOAD_PRS.out.prs_matrix
  ch_exclude = DOWNLOAD_EXCLUDE.out.file
  
  // ========================================
  // 1) Prepare cohorts with PRS + apply exclusions
  // ========================================
  PREPARE_PRS_COHORTS(ch_prs, ch_covars, ch_exclude)
  
  // Parse traits manifest
  ch_traits = PREPARE_PRS_COHORTS.out.traits_manifest
    .splitCsv(sep: '\t', header: true)
    .map { row -> 
      tuple(
        row.trait as String, 
        row.trait_tag as String, 
        row.iv_col as String
      ) 
    }
  
  // Parse cohort files into (stratum, file) tuples
  ch_cohorts = PREPARE_PRS_COHORTS.out.cohorts
    .flatten()
    .map { f ->
      def base = f.getName()
      // Extract stratum from filename: cohort_with_covariates.prs.{stratum}.tsv
      def stratum = base
        .replaceFirst(/^cohort_with_covariates\.prs\./, '')
        .replaceFirst(/\.tsv$/, '')
      tuple(stratum, f)
    }
  
  // ========================================
  // 2) Filter counts for EACH cohort (FIXED ISSUE 2/3)
  // ========================================
  // Combine each cohort with raw counts: (stratum, cohort_file, counts_file)
  ch_cohort_counts_inputs = ch_cohorts
    .combine(ch_counts)
    .map { stratum, cohort_file, counts_file ->
      tuple(stratum, cohort_file, counts_file)
    }
  
  // Filter counts to match each cohort's person_ids
  FILTER_COUNTS_TO_COHORT(ch_cohort_counts_inputs)
  
  // Output: (stratum, cohort_file, filtered_counts_file)
  ch_cohorts_with_counts = FILTER_COUNTS_TO_COHORT.out.cohort_counts
  
  // ========================================
  // 3) Cartesian product AFTER filtering (FIXED ISSUE 5)
  // ========================================
  // Now create jobs: traits × (cohort + filtered_counts)
  ch_jobs = ch_traits
    .combine(ch_cohorts_with_counts)
    .map { trait, trait_tag, iv_col, stratum, cohort_file, counts_file ->
      tuple(trait, trait_tag, stratum, cohort_file, iv_col, counts_file)
    }
  
  // Debug output
//  ch_jobs.view { job ->
//    "JOB: trait=${job[0]}, tag=${job[1]}, stratum=${job[2]}, " +
//    "cohort=${job[3].getName()}, iv=${job[4]}, counts=${job[5].getName()}"
//  }
  
  // ========================================
  // 4) Run PheWAS
  // ========================================
  RUN_PHEWAS(ch_jobs)
  
  // ========================================
  // 5) Generate plots
  // ========================================
  CUSTOM_PLOT(RUN_PHEWAS.out.results)
  
  // ========================================
  // Final outputs
  // ========================================
  RUN_PHEWAS.out.results.view { "RESULT: ${it.getName()}" }
  CUSTOM_PLOT.out.pdf.view { "PLOT: ${it.getName()}" }
}