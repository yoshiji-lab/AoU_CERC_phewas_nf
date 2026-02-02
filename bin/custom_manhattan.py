#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pick_col(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label-count", type=int, default=12)
    ap.add_argument("--converged-only", default="true")
    args = ap.parse_args()

    results_path = Path(args.results)
    out_path = Path(args.out)

    df = pd.read_csv(results_path, sep="\t", dtype={"phecode": "string"})

    conv_col = pick_col(df.columns, ["converged"])
    if conv_col is not None and str(args.converged_only).strip().lower() in ("1","true","t","yes","y","on"):
        df = df[df[conv_col].astype(str).str.lower().isin(["true", "1"])].copy()

    phecode_col = pick_col(df.columns, ["phecode", "phecode_id"])
    p_col = pick_col(df.columns, ["p_value", "pval", "pvalue", "p"])
    beta_col = pick_col(df.columns, ["beta", "coef", "estimate", "effect", "log_odds", "slope"])

    if phecode_col is None or p_col is None:
        raise SystemExit(f"Could not find required columns. Have: {df.columns.tolist()}")

    df["_p"] = pd.to_numeric(df[p_col], errors="coerce")
    df = df.dropna(subset=["_p"])

    def to_num(x):
        try:
            return float(x)
        except Exception:
            return np.nan

    df["_phe_num"] = df[phecode_col].map(to_num)
    df = df.sort_values(by=["_phe_num", phecode_col], na_position="last").reset_index(drop=True)

    x = np.arange(df.shape[0])
    y = -np.log10(df["_p"].clip(1e-300))

    if beta_col is not None:
        b = pd.to_numeric(df[beta_col], errors="coerce").fillna(0.0).values
    else:
        b = np.zeros(df.shape[0], dtype=float)

    pos = b >= 0
    neg = ~pos

    s = 20 + 80 * (np.abs(b) / (np.nanmax(np.abs(b)) + 1e-12))
    s = np.clip(s, 20, 200)

    plt.figure(figsize=(14, 5))
    plt.scatter(x[pos], y[pos], s=s[pos], marker="^", alpha=0.9)
    plt.scatter(x[neg], y[neg], s=s[neg], marker="v", alpha=0.9)

    m = max(1, df.shape[0])
    thr = -np.log10(0.05 / m)
    plt.axhline(thr, linestyle="--", linewidth=1)

    plt.xlabel("Phecode (sorted)")
    plt.ylabel("-log10(p)")

    k = max(0, int(args.label_count))
    if k > 0:
        top = df.nsmallest(k, "_p").copy()
        for i, r in top.iterrows():
            plt.annotate(str(r[phecode_col]), (x[i], y[i]),
                         textcoords="offset points", xytext=(0, 4),
                         ha="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
