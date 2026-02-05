#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def detect_age_col(cols):
    if "age_at_last_ehr_event" in cols:
        return "age_at_last_ehr_event"
    if "age_at_last_event" in cols:
        return "age_at_last_event"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phecode-version", default="1.2")
    ap.add_argument("--phecode-counts", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--iv-col", required=True)
    ap.add_argument("--sex-col", default="sex_at_birth")
    ap.add_argument("--age-col", default=None)
    ap.add_argument("--min-cases", type=int, default=50)
    ap.add_argument("--min-phecode-count", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from phetk.phewas import PheWAS

    cohort_path = Path(args.cohort)
    counts_path = Path(args.phecode_counts)
    out_path = Path(args.out)

    hdr = pd.read_csv(cohort_path, sep="\t", nrows=0)
    cols = hdr.columns.tolist()

    age_col = args.age_col or detect_age_col(cols)
    if age_col is None:
        raise SystemExit("Could not detect age column; provide --age-col.")

    pc_cols = [f"pc{i}" for i in range(1, 11) if f"pc{i}" in cols]
    covariate_cols = [age_col, args.sex_col] + pc_cols

    phewas = PheWAS(
        phecode_version=args.phecode_version,
        phecode_count_file_path=str(counts_path),
        cohort_file_path=str(cohort_path),
        covariate_cols=covariate_cols,
        independent_variable_of_interest=args.iv_col,
        sex_at_birth_col=args.sex_col,
        male_as_one=True,
        min_cases=int(args.min_cases),
        min_phecode_count=int(args.min_phecode_count),
        method="logit",
        output_file_path=str(out_path),
        suppress_warnings=True,
    )
    phewas.run()

    if not out_path.exists():
        cand = Path(str(out_path) + ".tsv")
        if cand.exists():
            cand.rename(out_path)
        else:
            tsvs = sorted(Path(".").glob("*.tsv"), key=lambda p: p.stat().st_mtime, reverse=True)
            if tsvs:
                tsvs[0].rename(out_path)
            else:
                raise SystemExit("PheWAS finished but no TSV output found.")

    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()