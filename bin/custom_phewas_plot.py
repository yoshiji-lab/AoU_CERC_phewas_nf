#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

def _pick_col(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None

def _to_bool_series(x):
    # Accept bools, 0/1, TRUE/FALSE, etc.
    if x.dtype == bool:
        return x
    s = x.astype(str).str.strip().str.lower()
    return s.isin(["true", "1", "t", "yes", "y", "on"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="PheTK PheWAS results TSV (or CSV)")
    ap.add_argument("--out", required=True, help="Output PDF path")
    ap.add_argument("--phecode-version", default="1.2")
    ap.add_argument("--converged-only", default="true")
    ap.add_argument("--label-count", type=int, default=15)
    ap.add_argument("--global-font-size", type=float, default=14)
    ap.add_argument("--fig-width", type=float, default=12)
    ap.add_argument("--fig-height", type=float, default=4)
    ap.add_argument("--axis-text-size", type=float, default=14)
    ap.add_argument("--title-text-size", type=float, default=18)
    ap.add_argument("--label-size", type=float, default=10)
    args = ap.parse_args()

    results_path = Path(args.results)
    out_pdf = Path(args.out)

    # Read results (PheTK writes TSV by default)
    sep = "\t" if results_path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(results_path, sep=sep, dtype={"phecode": "string"}, low_memory=False)

    # Ensure required columns exist (phewas_plot.py expects these names)
    # phecode
    if "phecode" not in df.columns:
        c = _pick_col(df.columns, ["phecode_id", "phecode_code"])
        if c:
            df["phecode"] = df[c].astype("string")

    # beta
    if "beta" not in df.columns:
        c = _pick_col(df.columns, ["coef", "estimate", "effect", "log_odds", "slope"])
        if c:
            df["beta"] = pd.to_numeric(df[c], errors="coerce")

    # p_value
    if "p_value" not in df.columns:
        c = _pick_col(df.columns, ["pval", "pvalue", "p"])
        if c:
            df["p_value"] = pd.to_numeric(df[c], errors="coerce")

    # neg_log_p_value
    if "neg_log_p_value" not in df.columns and "p_value" in df.columns:
        p = pd.to_numeric(df["p_value"], errors="coerce")
        # protect against 0 -> inf
        p = p.clip(lower=np.nextafter(0, 1), upper=1.0)
        df["neg_log_p_value"] = -np.log10(p)

    # converged
    if "converged" not in df.columns:
        df["converged"] = True
    else:
        df["converged"] = _to_bool_series(df["converged"])

    # phecode_string (label column default)
    if "phecode_string" not in df.columns:
        c = _pick_col(df.columns, ["phenotype", "description", "phecode_name"])
        if c:
            df["phecode_string"] = df[c].astype(str)
        else:
            df["phecode_string"] = df["phecode"].astype(str)

    # phecode_category (used for colors)
    if "phecode_category" not in df.columns:
        c = _pick_col(df.columns, ["category", "phecode_group", "group"])
        if c:
            df["phecode_category"] = df[c].astype(str)
        else:
            # fallback: put everything in one category
            df["phecode_category"] = "all"

    # Write a temp CSV next to output (polars reader in phewas_plot expects CSV)
    csv_path = out_pdf.with_suffix(".plot_input.csv")
    df.to_csv(csv_path, index=False)

    # Import plotting code shipped in pipeline bin/
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from phewas_plot import Plot, set_global_font_size
    except Exception as e:
        raise SystemExit(
            "Failed to import phewas_plot.py. "
            "Make sure bin/phewas_plot.py exists. "
            f"Import error: {e}"
        )

    # Run plot
    set_global_font_size(args.global_font_size)

    try:
        p = Plot(
            str(csv_path),
            converged_only=str(args.converged_only).strip().lower() in ("1","true","t","yes","y","on"),
            phecode_version=str(args.phecode_version),
        )
        p.manhattan(
            fig_width=float(args.fig_width),
            fig_height=float(args.fig_height),
            axis_text_size=float(args.axis_text_size),
            title_text_size=float(args.title_text_size),
            label_size=float(args.label_size),
            label_count=int(args.label_count),
            output_file_name=str(out_pdf),
            output_file_type="pdf",
            save_plot=True,
        )
    except ModuleNotFoundError as e:
        # Most common missing deps: polars, adjustText
        raise SystemExit(
            f"Missing Python dependency: {e}.\n"
            "Install required deps in your notebook environment:\n"
            "  pip install -U polars adjustText\n"
            "Matplotlib/numpy/pandas are usually preinstalled."
        )

    print("Wrote:", out_pdf)

if __name__ == "__main__":
    main()
