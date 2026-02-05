process FILTER_COUNTS_TO_COHORT {
  tag "${stratum}"

  input:
    tuple val(stratum), path(cohort_tsv), path(counts_tsv)

  output:
    tuple val(stratum), path(cohort_tsv), path("counts.${stratum}.filtered.tsv"), emit: cohort_counts
    path "filter_counts.${stratum}.log", optional: true, emit: log

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/filter_counts_to_cohort.py \
      --counts "${counts_tsv}" \
      --cohort "${cohort_tsv}" \
      --out "counts.${stratum}.filtered.tsv" \
      --counts-person-id-col "person_id" \
      --cohort-person-id-col "${params.person_id_col}" \
      --log "filter_counts.${stratum}.log"
    """
}