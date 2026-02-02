nextflow.enable.dsl = 2

include { DOWNLOAD_PRS }                 from './modules/download_prs'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COUNTS } from './modules/download_gcs_or_local'
include { DOWNLOAD_GCS_OR_LOCAL as DOWNLOAD_COVARS } from './modules/download_gcs_or_local'
include { PREPARE_PRS_COHORTS }          from './modules/prepare_prs_cohorts'
include { RUN_PHEWAS }                   from './modules/run_phewas'
include { CUSTOM_PLOT }                  from './modules/custom_plot'

workflow {

  // strings (local path or gs://)
  ch_counts_in = Channel.value(params.phecode_counts)
  ch_covars_in = Channel.value(params.cohort_covariates)
  ch_prs_in    = Channel.value(params.prs_matrix)

  // 0) Stage inputs (call processes)
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
    
  //ch_traits.view { "TRAIT: $it" }  

  ch_cohorts = PREPARE_PRS_COHORTS.out.cohorts
    .flatten()
    .map { f ->
        def base = f.getName()
        def stratum = base.replaceFirst(/^cohort_with_covariates\.prs\./, '').replaceFirst(/\.tsv$/, '')
        tuple(stratum, f)
    }
    
  def sel = (params.ancestry_select ?: 'ALL').toString().trim()
  def selU = sel.toUpperCase()
  def allowed = null

  if( selU != 'ALL' && selU != 'POOLED' ) {
    allowed = sel.split(',')
                 .collect{ it.trim().toLowerCase() }
                 .findAll{ it }
                 .toSet()
    allowed.add('pooled')   // keep pooled by default
  }

  ch_cohorts = ch_cohorts.filter { stratum, f ->
    if( selU == 'POOLED' ) return stratum == 'pooled'
    if( selU == 'ALL' )    return true
    return allowed.contains(stratum.toString().toLowerCase())
  }

    
  //ch_cohorts.view { "COHORT: $it" } 

  //ch_runs = ch_traits
   // .cross(ch_cohorts)
    //.map { t, c ->
    //    def (trait, trait_tag, iv_col) = t
    //    def (stratum, cohort_file) = c
    //    tuple(trait, trait_tag, stratum, cohort_file, iv_col)
   // }

  ch_runs = ch_traits
    .combine(ch_cohorts)
    .map { trait, trait_tag, iv_col, stratum, cohort_file ->
        tuple(trait, trait_tag, stratum, cohort_file, iv_col)
    }
  
  //ch_runs.view { "RUN: $it" }
  
  //ch_counts.view { "COUNTS: $it" }  
    
  ch_jobs = ch_runs
    .combine(ch_counts)
    .map { trait, trait_tag, stratum, cohort_file, iv_col, counts_file ->
        tuple(trait, trait_tag, stratum, cohort_file, iv_col, counts_file)
    }
    
  ch_jobs.view { "JOB: $it" }  
  
  // Attach counts_file to every run (counts has 1 item; cross replicates it)
  //ch_jobs = ch_runs
  //  .cross(ch_counts)
  //  .map { run, counts_file ->
  //      def (trait, trait_tag, stratum, cohort_file, iv_col) = run
  //      tuple(trait, trait_tag, stratum, cohort_file, iv_col, counts_file)
  //  }

  // 2) Run + plot
  RUN_PHEWAS(ch_jobs)
  CUSTOM_PLOT(RUN_PHEWAS.out.results)

  // Final sinks
  RUN_PHEWAS.out.results.view { "RESULT\t$it" }
  CUSTOM_PLOT.out.pdf.view    { "PLOT\t$it" }
}
