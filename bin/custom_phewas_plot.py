#!/usr/bin/env python3
"""
custom_phewas_plot.py - Publication-ready PheWAS Manhattan plots

Adds an "Infinity" line (dash-dot) when neg_log_p_value contains +inf
(e.g. p_value underflowed to 0 => -log10(0)=inf).
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd

# Font configuration
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def set_global_font_size(size: int) -> None:
    """Set all matplotlib font sizes."""
    keys = [
        "font.size",
        "axes.titlesize",
        "axes.labelsize",
        "xtick.labelsize",
        "ytick.labelsize",
        "legend.fontsize",
        "figure.titlesize",
    ]
    for k in keys:
        mpl.rcParams[k] = size


def _parse_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="PheWAS results TSV")
    ap.add_argument("--out", required=True, help="Output PDF filename")
    ap.add_argument("--phecode-version", default="1.2")
    ap.add_argument("--converged-only", type=_parse_bool, default=True)
    ap.add_argument("--label-count", type=int, default=10)
    ap.add_argument("--global-font-size", type=int, default=10)

    ap.add_argument("--fig-width", type=float, default=14)
    ap.add_argument("--fig-height", type=float, default=10)
    ap.add_argument("--axis-text-size", type=int, default=8)
    ap.add_argument("--title-text-size", type=int, default=10)
    ap.add_argument("--label-size", type=int, default=8)
    ap.add_argument("--legend-size", type=int, default=6)
    ap.add_argument("--dpi", type=int, default=150)

    # New: allow toggling infinity line (default True to match yoshiji-lab behavior)
    ap.add_argument("--infinity-line", type=_parse_bool, default=True)

    args = ap.parse_args()

    set_global_font_size(args.global_font_size)

    # Load results
    df = pd.read_csv(args.results, sep="\t", dtype={"phecode": "string"})

    if len(df) == 0:
        print("WARNING: Results file is empty; nothing to plot.")
        sys.exit(0)

    # Normalize converged column (bool/strings/0-1) then filter if requested
    if "converged" in df.columns:
        df["converged"] = df["converged"].map(
            {
                True: True,
                False: False,
                "True": True,
                "False": False,
                "TRUE": True,
                "FALSE": False,
                1: True,
                0: False,
                "1": True,
                "0": False,
            }
        ).fillna(df["converged"])
        # If still not boolean-like, leave as-is; filter below may do nothing.

    if args.converged_only and "converged" in df.columns:
        df = df[df["converged"] == True].copy()

    if len(df) == 0:
        print("WARNING: No data to plot after converged filtering.")
        sys.exit(0)

    # Ensure neg_log_p_value exists and is numeric
    if "neg_log_p_value" not in df.columns:
        if "p_value" not in df.columns:
            raise ValueError("Missing both 'neg_log_p_value' and 'p_value' columns in results.")
        p = pd.to_numeric(df["p_value"], errors="coerce")
        df["neg_log_p_value"] = -np.log10(p)
    else:
        df["neg_log_p_value"] = pd.to_numeric(df["neg_log_p_value"], errors="coerce")

    # Significance thresholds
    bonferroni = -np.log10(0.05 / max(len(df), 1))
    nominal = -np.log10(0.05)

    # Handle +infinity values and compute inf proxy for plotting
    has_inf = np.isinf(df["neg_log_p_value"]).any()
    inf_proxy = None
    if has_inf:
        finite_mask = np.isfinite(df["neg_log_p_value"])
        if finite_mask.any():
            max_finite = float(df.loc[finite_mask, "neg_log_p_value"].max())
        else:
            # Extremely pathological: all values are inf/NaN
            max_finite = float(max(10.0, bonferroni, nominal))
        inf_proxy = max_finite * 1.2
        df.loc[np.isinf(df["neg_log_p_value"]), "neg_log_p_value"] = inf_proxy

    # Required columns check (fail loudly, since plot assumptions depend on them)
    needed = ["phecode_category", "beta", "phecode"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in results: {missing}")

    # Sort by category and p-value
    df = df.sort_values(["phecode_category", "neg_log_p_value"], ascending=[True, False]).reset_index(drop=True)
    df["phecode_index"] = range(1, len(df) + 1)

    # Color palette (matches phewas_plot.py)
    color_palette = [
        "blue", "indianred", "darkcyan", "goldenrod", "darkblue",
        "magenta", "green", "red", "darkturquoise", "olive",
        "black", "royalblue", "maroon", "darkolivegreen",
        "coral", "purple", "gray"
    ]

    categories = sorted(df["phecode_category"].dropna().astype(str).unique())
    color_dict = {cat: color_palette[i % len(color_palette)] for i, cat in enumerate(categories)}
    df["color"] = df["phecode_category"].astype(str).map(color_dict)

    # Title from filename
    basename = Path(args.results).stem
    title_bits = basename.replace("phewas.", "", 1).split(".")
    trait_name = title_bits[0].replace("_", " ") if title_bits else basename
    stratum = title_bits[1] if len(title_bits) > 1 else "pooled"
    title = f"{trait_name} PheWAS — {stratum} (phecode v{args.phecode_version})"

    # Create figure
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    plt.subplots_adjust(bottom=0.15, left=0.08, right=0.88, top=0.93)
    plt.title(title, weight="bold", size=args.title_text_size)

    # Plot points (black edges)
    pos_beta = df[df["beta"] >= 0]
    neg_beta = df[df["beta"] < 0]

    for data, marker, alpha in [(pos_beta, "^", 0.7), (neg_beta, "v", 0.3)]:
        for cat in categories:
            cat_data = data[data["phecode_category"].astype(str) == cat]
            if len(cat_data) == 0:
                continue
            ax.scatter(
                cat_data["phecode_index"],
                cat_data["neg_log_p_value"],
                c=cat_data["color"],
                marker=marker,
                s=50,
                alpha=alpha,
                edgecolors="black",
                linewidths=0.5,
            )

    # Significance lines
    offset = 9
    x_min = int(df["phecode_index"].min()) - offset - 1
    x_max = int(df["phecode_index"].max()) + offset + 1

    ax.axhline(nominal, color="red", lw=1)
    ax.axhline(bonferroni, color="green", lw=1)

    # NEW: Infinity line (dash-dot), matching yoshiji-lab style
    if args.infinity_line and has_inf and inf_proxy is not None and np.isfinite(inf_proxy):
        ax.hlines(
            inf_proxy * 0.98,
            x_min,
            x_max,
            colors="blue",
            linestyle="dashdot",
            lw=1,
        )

    ax.set_ylabel(r"$-\log_{10}$(p-value)", size=args.axis_text_size)
    ax.set_xlim(x_min, x_max)

    # X-axis category labels
    cat_positions = df.groupby("phecode_category")["phecode_index"].mean()
    ax.set_xticks(cat_positions.values)
    ax.set_xticklabels(cat_positions.index, rotation=45, ha="right", size=args.axis_text_size)

    for tick_label, cat in zip(ax.get_xticklabels(), cat_positions.index):
        tick_label.set_color(color_dict.get(str(cat), "black"))

    # Optional: adjustText for labels (kept as in your original)
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False
        print("WARNING: adjustText not available, labels may overlap")

    # Label top phecodes
    top_phecodes = df.nlargest(args.label_count, "neg_log_p_value")
    texts = []
    for _, row in top_phecodes.iterrows():
        label_text = row.get("phecode_string", str(row["phecode"]))
        if isinstance(label_text, str) and len(label_text) > 30:
            words = label_text.split()
            lines, current, length = [], [], 0
            for word in words:
                if length + len(word) > 30 and current:
                    lines.append(" ".join(current))
                    current, length = [word], len(word)
                else:
                    current.append(word)
                    length += len(word) + 1
            if current:
                lines.append(" ".join(current))
            label_text = "\n".join(lines)

        text = ax.text(
            row["phecode_index"],
            row["neg_log_p_value"],
            label_text,
            color=row["color"],
            size=args.label_size,
            weight="normal",
            bbox=dict(facecolor="white", edgecolor="none", boxstyle="round", alpha=0.3, pad=0.1),
            ha="center",
        )
        texts.append(text)

    if use_adjust_text and texts:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.5),
            expand_points=(1.5, 1.5),
            expand_text=(1.2, 1.2),
            force_points=0.5,
            force_text=0.5,
        )

    # Legend (explicit handles)
    legend_elements = [
        Line2D([0], [0], color="green", lw=1, label="Bonferroni\nCorrection"),
        Line2D([0], [0], color="red", lw=1, label="Nominal\nSignificance"),
    ]

    if args.infinity_line and has_inf and inf_proxy is not None and np.isfinite(inf_proxy):
        legend_elements.append(
            Line2D([0], [0], color="blue", lw=1, linestyle="dashdot", label="Infinity")
        )

    legend_elements.extend([
        Line2D([0], [0], marker="^", color="white",
               markerfacecolor="blue",
               markeredgecolor="black", markeredgewidth=0.5,
               markersize=args.legend_size, alpha=0.7,
               label="Increased\nRisk"),
        Line2D([0], [0], marker="v", color="white",
               markerfacecolor="blue",
               markeredgecolor="black", markeredgewidth=0.5,
               markersize=args.legend_size, alpha=0.3,
               label="Decreased\nRisk"),
    ])

    ax.legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize=args.legend_size,
        frameon=True,
    )

    plt.savefig(args.out, dpi=args.dpi)
    print(f"Plot saved: {args.out}")
    plt.close()


if __name__ == "__main__":
    main()
