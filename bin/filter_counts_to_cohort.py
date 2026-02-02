#!/usr/bin/env python3
"""
filter_counts_to_cohort.py

Filters a large AoU phecode counts TSV to only those person_ids present in a cohort TSV.

Designed to be fast and memory-safe by streaming the counts TSV line-by-line.

Inputs:
  --counts: counts TSV (must include person_id column)
  --cohort: cohort TSV (must include person_id column)
Output:
  --out: filtered counts TSV (same columns as input counts)

Notes:
  - Keeps header/columns unchanged.
  - Uses a Python set of cohort ids in memory (typically feasible).
"""

import argparse
import csv
from pathlib import Path
from typing import Set


def load_ids(cohort_path: Path, cohort_id_col: str) -> Set[str]:
    ids: Set[str] = set()
    with cohort_path.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: cohort file has no header")
        if cohort_id_col not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: cohort missing column '{cohort_id_col}'. Columns: {reader.fieldnames}"
            )
        for row in reader:
            v = row.get(cohort_id_col)
            if v is None:
                continue
            v = str(v).strip()
            if v:
                ids.add(v)
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="Counts TSV (e.g. aou_phecode_counts.tsv)")
    ap.add_argument("--cohort", required=True, help="Cohort TSV (e.g. cohort_with_covariates.prs.pooled.tsv)")
    ap.add_argument("--out", required=True, help="Output filtered counts TSV")
    ap.add_argument("--counts-person-id-col", default="person_id", help="Person id column in counts file")
    ap.add_argument("--cohort-person-id-col", default="person_id", help="Person id column in cohort file")
    ap.add_argument("--log", default="filter_counts.log", help="Log file path")
    args = ap.parse_args()

    counts_path = Path(args.counts)
    cohort_path = Path(args.cohort)
    out_path = Path(args.out)
    log_path = Path(args.log)

    def log(msg: str):
        with log_path.open("a") as lf:
            lf.write(msg + "\n")

    if not counts_path.exists():
        raise SystemExit(f"ERROR: counts file not found: {counts_path}")
    if not cohort_path.exists():
        raise SystemExit(f"ERROR: cohort file not found: {cohort_path}")

    log("=== FILTER_COUNTS_TO_COHORT ===")
    log(f"counts: {counts_path}")
    log(f"cohort: {cohort_path}")
    log(f"counts_person_id_col: {args.counts_person_id_col}")
    log(f"cohort_person_id_col: {args.cohort_person_id_col}")

    cohort_ids = load_ids(cohort_path, args.cohort_person_id_col)
    log(f"Loaded cohort ids: {len(cohort_ids)}")

    kept = 0
    total = 0

    with counts_path.open("r", newline="") as fin, out_path.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: counts file has no header")
        if args.counts_person_id_col not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: counts missing column '{args.counts_person_id_col}'. Columns: {reader.fieldnames}"
            )

        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            total += 1
            pid = row.get(args.counts_person_id_col)
            if pid is None:
                continue
            pid = str(pid).strip()
            if pid in cohort_ids:
                writer.writerow(row)
                kept += 1

    log(f"Rows read: {total}")
    log(f"Rows kept: {kept}")
    log("Done.")


if __name__ == "__main__":
    main()

