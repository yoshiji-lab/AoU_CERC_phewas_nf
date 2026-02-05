process RUN_PHEWAS {
  tag "${trait_tag}.${stratum}"

  input:
    tuple val(trait), val(trait_tag), val(stratum), path(cohort_file), val(iv_col), path(counts_file)

  output:
    path "phewas.${trait_tag}.${stratum}.tsv", emit: results

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/run_phewas.py \
      --phecode-version "${params.phecode_version}" \
      --phecode-counts "${counts_file}" \
      --cohort "${cohort_file}" \
      --iv-col "${iv_col}" \
      --sex-col "sex_at_birth" \
      --min-cases "${params.min_cases}" \
      --min-phecode-count "${params.min_phecode_count}" \
      --out "phewas.${trait_tag}.${stratum}.tsv"
    """
}