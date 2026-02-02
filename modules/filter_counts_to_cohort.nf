process FILTER_COUNTS_TO_COHORT {
  tag "filter_counts"

  input:
    path counts_tsv
    path cohort_tsv

  output:
    path "aou_phecode_counts.filtered.tsv", emit: file
    path "filter_counts.log", optional: true, emit: log

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/filter_counts_to_cohort.py \
      --counts "${counts_tsv}" \
      --cohort "${cohort_tsv}" \
      --out "aou_phecode_counts.filtered.tsv" \
      --counts-person-id-col "person_id" \
      --cohort-person-id-col "${params.person_id_col}" \
      --log "filter_counts.log"
    """
}
