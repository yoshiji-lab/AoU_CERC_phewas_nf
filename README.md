# AoU PRS → PheWAS (Nextflow DSL2)

This pipeline implements your **Part B** workflow:

1) Stage inputs (local path or `gs://...`):
   - PRS matrix TSV (columns: `FID`, `IID`, PRS trait columns)
   - cohort covariates TSV (from PheTK `Cohort.add_covariates`)
   - phecode counts TSV (from PheTK `Phecode.count_phecode`)
2) Build pooled + ancestry-stratified cohorts with INT-transformed PRS(s)
3) Run logistic PheWAS for each (trait × stratum)
4) Produce a Manhattan PDF per result TSV (custom matplotlib script)

## Requirements

- Nextflow (tested with v24.10.4)
- Python 3
- `phetk` installed in the environment where you run `RUN_PHEWAS`
- `gsutil` on PATH if you use `gs://...` inputs

## Run

```bash
./nextflow run main.nf -c nextflow.config -resume \
  --prs_matrix "gs://.../prs_matrix_ALL_TRAITS.tsv" \
  --cohort_covariates "gs://.../rs2925979_cohort_with_covariates.tsv" \
  --phecode_counts "gs://.../aou_phecode_counts.tsv" \
  --prs_traits "Gynoid_fat" \
  --outdir "analysis/phewas_prs_nf_out"
```

Outputs are published to:

- `${outdir}/inputs/`  (prepared cohort files + traits_manifest.tsv + prep.log)
- `${outdir}/results/` (phewas.*.tsv)
- `${outdir}/plots/`   (manhattan.*.pdf)

## Notes

- `prs_traits` accepts:
  - comma-separated list
  - `ALL`
  - a file path containing one trait per line
- INT is Blom rank-INT with tiny jitter (seeded for reproducibility).
- If you prefer your lab’s custom plotter, swap `bin/custom_manhattan.py`.
