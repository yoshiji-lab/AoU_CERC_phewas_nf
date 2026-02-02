nextflow.enable.dsl = 2
include { DOWNLOAD_PRS }                 from './modules/download_prs'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COUNTS } from './modules/download_gcs_or_local'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COVARS } from './modules/download_gcs_or_local'
include { PREPARE_PRS_COHORTS }          from './modules/prepare_prs_cohorts'
include { FILTER_COUNTS_TO_COHORT }      from './modules/filter_counts_to_cohort'
include { RUN_PHEWAS }                   from './modules/run_phewas'
include { CUSTOM_PLOT }                  from './modules/custom_plot'

workflow {
  // strings (local path or gs://)
  ch_counts_in = Channel.value(params.phecode_counts)
  ch_covars_in = Channel.value(params.cohort_covariates)
  ch_prs_in    = Channel.value(params.prs_matrix)
  
  // 0) Stage inputs locally (call processes)
  DOWNLOAD_COUNTS(ch_counts_in)
  DOWNLOAD_COVARS(ch_covars_in)
  DOWNLOAD_PRS(ch_prs_in)
  
  // Grab emitted file channels
  ch_counts = DOWNLOAD_COUNTS.out.file
  ch_covars = DOWNLOAD_COVARS.out.file
  ch_prs    = DOWNLOAD_PRS.out.prs_matrix
  
  // 1) Prepare cohorts + traits manifest
  PREPARE_PRS_COHORTS(ch_prs, ch_covars)
  
  ch_traits = PREPARE_PRS_COHORTS.out.traits_manifest
    .splitCsv(sep: '\t', header: true)
    .map { row -> tuple(row.trait as String, row.trait_tag as String, row.iv_col as String) }
  
  ch_cohorts = PREPARE_PRS_COHORTS.out.cohorts
    .flatten()
    .map { f ->
        def base = f.getName()
        def stratum = base.replaceFirst(/^cohort_with_covariates\.prs\./, '').replaceFirst(/\.tsv$/, '')
        tuple(stratum, f)
    }
  
  // Filter counts to pooled cohort only
  ch_pooled = ch_cohorts
    .filter { stratum, f -> stratum == 'pooled' }
    .map { stratum, f -> f }
  
  FILTER_COUNTS_TO_COHORT(ch_counts, ch_pooled)
  ch_counts_filtered = FILTER_COUNTS_TO_COHORT.out.file
  
  // Create jobs: every trait × every cohort stratum × filtered counts
  // Use combine() for cartesian product (correct for DSL2)
  ch_runs = ch_traits
    .combine(ch_cohorts)
    .map { trait, trait_tag, iv_col, stratum, cohort_file ->
        tuple(trait, trait_tag, stratum, cohort_file, iv_col)
    }
  
  ch_jobs = ch_runs
    .combine(ch_counts_filtered)
    .map { trait, trait_tag, stratum, cohort_file, iv_col, counts_file ->
        tuple(trait, trait_tag, stratum, cohort_file, iv_col, counts_file)
    }
  
  ch_jobs.view { "JOB\t$it" }
  
  RUN_PHEWAS(ch_jobs)
  CUSTOM_PLOT(RUN_PHEWAS.out.results)
  
  // Final sinks
  RUN_PHEWAS.out.results.view { "RESULT\t$it" }
  CUSTOM_PLOT.out.pdf.view    { "PLOT\t$it" }
}