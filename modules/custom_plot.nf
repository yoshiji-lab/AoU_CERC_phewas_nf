process CUSTOM_PLOT {
  tag "${results_tsv.baseName}"

  input:
    path results_tsv

  output:
    path "*.pdf", emit: pdf

  script:
    """
    set -euo pipefail

    python3 ${projectDir}/bin/custom_phewas_plot.py \
      --results "${results_tsv}" \
      --out "manhattan.${results_tsv.baseName}.pdf" \
      --phecode-version "${params.phecode_version}" \
      --converged-only "${params.converged_only}" \
      --label-count "${params.label_count}" \
      --global-font-size "${params.plot_global_font_size}" \
      --fig-width "${params.plot_fig_width}" \
      --fig-height "${params.plot_fig_height}" \
      --axis-text-size "${params.plot_axis_text_size}" \
      --title-text-size "${params.plot_title_text_size}" \
      --label-size "${params.plot_label_size}" \
      --legend-size "${params.plot_legend_size}" \
      --dpi "${params.plot_dpi}"
    """
}