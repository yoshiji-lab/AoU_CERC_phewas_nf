process PREPARE_PRS_COHORTS {
  tag "prepare_cohorts"

  input:
    path prs_matrix, stageAs: 'prs_input.tsv'
    path cohort_covariates, stageAs: 'cohort_input.tsv'
    path exclude_ids, stageAs: 'exclude_input.tsv'

  output:
    path "traits_manifest.tsv", emit: traits_manifest
    path "cohort_with_covariates.prs.*.tsv", emit: cohorts
    path "prep.log", emit: log

  script:
    """
    set -euo pipefail
    
    python3 ${projectDir}/bin/prepare_prs_cohorts.py \
      --prs-matrix "prs_input.tsv" \
      --cohort-covariates "cohort_input.tsv" \
      --exclude-ids "exclude_input.tsv" \
      --prs-traits "${params.prs_traits}" \
      --out-prefix "cohort_with_covariates.prs" \
      --stratify-by-ancestry "${params.stratify_by_ancestry}" \
      --ancestry-col "${params.ancestry_col}" \
      --person-id-col "${params.person_id_col}" \
      --ancestry-select "${params.ancestry_select}" \
      --int-seed "${params.int_seed}" \
      --int-jitter-sd "${params.int_jitter_sd}" \
      --log "prep.log"
    """
}