process DOWNLOAD_GCS_OR_LOCAL {
  tag "${uri}"

  input:
    val uri

  output:
    path "input", emit: file

  script:
    """
    set -euo pipefail

    if [[ "${uri}" == gs://* ]]; then
      command -v gsutil >/dev/null 2>&1 || { echo "ERROR: gsutil not found" >&2; exit 2; }
      gsutil cp "${uri}" input
    else
      [[ -f "${uri}" ]] || { echo "ERROR: local file not found: ${uri}" >&2; exit 2; }
      ln -sf "${uri}" input || cp -f "${uri}" input
    fi
    """
}

