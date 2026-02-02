process DOWNLOAD_PRS {
  tag "prs_matrix"

  input:
    val prs_in

  output:
    path "prs_matrix.tsv", emit: prs_matrix

  script:
    """
    set -euo pipefail

    if [[ "${prs_in}" == gs://* ]]; then
      command -v gsutil >/dev/null 2>&1 || { echo "ERROR: gsutil not found" >&2; exit 2; }
      gsutil cp "${prs_in}" prs_matrix.tsv
    else
      [[ -f "${prs_in}" ]] || { echo "ERROR: local prs_matrix not found: ${prs_in}" >&2; exit 2; }
      ln -sf "${prs_in}" prs_matrix.tsv || cp -f "${prs_in}" prs_matrix.tsv
    fi
    """
}
