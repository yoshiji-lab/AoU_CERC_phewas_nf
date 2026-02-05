process DOWNLOAD_EXCLUDE {
  tag "exclude_ids"

  input:
    val exclude_in

  output:
    path "exclude_ids.tsv", emit: exclude_ids

  script:
    """
    set -euo pipefail

    if [[ "${exclude_in}" == gs://* ]]; then
      command -v gsutil >/dev/null 2>&1 || { echo "ERROR: gsutil not found" >&2; exit 2; }
      gsutil cp "${exclude_in}" exclude_ids.tsv
    else
      [[ -f "${exclude_in}" ]] || { echo "ERROR: local exclude file not found: ${exclude_in}" >&2; exit 2; }
      ln -sf "${exclude_in}" exclude_ids.tsv || cp -f "${exclude_in}" exclude_ids.tsv
    fi
    """
}