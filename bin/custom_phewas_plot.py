#!/usr/bin/env python3
"""
custom_phewas_plot.py - Publication-ready PheWAS Manhattan plots
"""

import argparse
import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

# Font configuration
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

def set_global_font_size(size):
    """Set all matplotlib font sizes"""
    keys = ["font.size", "axes.titlesize", "axes.labelsize",
            "xtick.labelsize", "ytick.labelsize",
            "legend.fontsize", "figure.titlesize"]
    for k in keys:
        mpl.rcParams[k] = size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="PheWAS results TSV")
    ap.add_argument("--out", required=True, help="Output PDF filename")
    ap.add_argument("--phecode-version", default="1.2")
    ap.add_argument("--converged-only", type=lambda x: str(x).lower() == "true", default=True)
    ap.add_argument("--label-count", type=int, default=10)  # Increased from 10
    ap.add_argument("--global-font-size", type=int, default=10)  # Increased from 10
    ap.add_argument("--fig-width", type=float, default=14)
    ap.add_argument("--fig-height", type=float, default=10)  # Updated to 8
    ap.add_argument("--axis-text-size", type=int, default=8)  # Increased from 8
    ap.add_argument("--title-text-size", type=int, default=10)  # Increased from 10
    ap.add_argument("--label-size", type=int, default=8)  # Increased from 8
    ap.add_argument("--legend-size", type=int, default=6)  # Increased from 6
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    # Set global font size
    set_global_font_size(args.global_font_size)

    # Load results
    df = pd.read_csv(args.results, sep="\t", dtype={"phecode": str})
    
    # Filter converged if requested
    if args.converged_only and "converged" in df.columns:
        df = df[df["converged"] == True].copy()
    
    if len(df) == 0:
        print("WARNING: No data to plot after filtering")
        sys.exit(0)
    
    # Calculate significance thresholds
    bonferroni = -np.log10(0.05 / len(df))
    nominal = -np.log10(0.05)
    
    # Handle infinite p-values
    max_finite = df[df["neg_log_p_value"] != np.inf]["neg_log_p_value"].max()
    if (df["neg_log_p_value"] == np.inf).any():
        df.loc[df["neg_log_p_value"] == np.inf, "neg_log_p_value"] = max_finite * 1.2
    
    # Sort by category and p-value
    df = df.sort_values(["phecode_category", "neg_log_p_value"], 
                        ascending=[True, False]).reset_index(drop=True)
    df["phecode_index"] = range(1, len(df) + 1)
    
    # Color palette (matches phewas_plot.py)
    color_palette = [
        "blue", "indianred", "darkcyan", "goldenrod", "darkblue",
        "magenta", "green", "red", "darkturquoise", "olive",
        "black", "royalblue", "maroon", "darkolivegreen",
        "coral", "purple", "gray"
    ]
    
    categories = sorted(df["phecode_category"].unique())
    color_dict = {cat: color_palette[i % len(color_palette)] 
                  for i, cat in enumerate(categories)}
    df["color"] = df["phecode_category"].map(color_dict)
    
    # FIX: Extract title from filename - remove duplicate "pooled"
    basename = Path(args.results).stem
    # Remove "phewas." prefix
    title = basename.replace("phewas.", "", 1)
    # Split and clean
    parts = title.split(".")
    trait_name = parts[0].replace("_", " ") if parts else title
    stratum = parts[1] if len(parts) > 1 else "pooled"
    
    # Create clean title - avoid duplicate "pooled"
    title = f"{trait_name} PheWAS — {stratum} (phecode v{args.phecode_version})"
    
    # Create figure with proper height
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    
    # Add margins to prevent label cutoff
    plt.subplots_adjust(bottom=0.15, left=0.08, right=0.88, top=0.93)
    
    # Set title
    plt.title(title, weight="bold", size=args.title_text_size)
    
    # FIX: Plot data points with BLACK EDGES on markers
    pos_beta = df[df["beta"] >= 0]
    neg_beta = df[df["beta"] < 0]
    
    for data, marker, alpha in [(pos_beta, "^", 0.7), (neg_beta, "v", 0.3)]:
        for cat in categories:
            cat_data = data[data["phecode_category"] == cat]
            if len(cat_data) > 0:
                ax.scatter(
                    cat_data["phecode_index"],
                    cat_data["neg_log_p_value"],
                    c=cat_data["color"],
                    marker=marker,
                    s=50,  # Slightly larger markers
                    alpha=alpha,
                    edgecolors="black",  # BLACK EDGE
                    linewidths=0.5       # Edge width
                )
    
    # Add significance lines
    offset = 9
    x_min = df["phecode_index"].min() - offset - 1
    x_max = df["phecode_index"].max() + offset + 1
    
    ax.axhline(nominal, color="red", lw=1, label="Nominal Significance")
    ax.axhline(bonferroni, color="green", lw=1, label="Bonferroni Correction")
    
    # Set axis labels
    ax.set_ylabel(r"$-\log_{10}$(p-value)", size=args.axis_text_size)
    ax.set_xlim(x_min, x_max)
    
    # FIX: X-axis category labels with MORE SPACING
    cat_positions = df.groupby("phecode_category")["phecode_index"].mean()
    ax.set_xticks(cat_positions.values)
    ax.set_xticklabels(
        cat_positions.index,
        rotation=45,
        ha="right",
        size=args.axis_text_size
    )
    
    # Color the x-tick labels
    for tick_label, cat in zip(ax.get_xticklabels(), cat_positions.index):
        tick_label.set_color(color_dict[cat])
    
    # FIX: Use adjustText for better label positioning (prevents overlap)
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False
        print("WARNING: adjustText not available, labels may overlap")
    
    # Add labels for top significant phecodes
    top_phecodes = df.nlargest(args.label_count, "neg_log_p_value")
    
    texts = []
    for _, row in top_phecodes.iterrows():
        label_text = row.get("phecode_string", str(row["phecode"]))
        # Split long labels
        if len(label_text) > 30:
            words = label_text.split()
            lines = []
            current = []
            length = 0
            for word in words:
                if length + len(word) > 30 and current:
                    lines.append(" ".join(current))
                    current = [word]
                    length = len(word)
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
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                boxstyle="round",
                alpha=0.3,
                pad=0.1
            ),
            ha="center"
        )
        texts.append(text)
    
    # Adjust text positions to avoid overlap
    if use_adjust_text and texts:
        adjust_text(
            texts,
            arrowprops=dict(
                arrowstyle="->",
                color="gray",
                lw=0.5
            ),
            expand_points=(1.5, 1.5),  # More spacing
            expand_text=(1.2, 1.2),
            force_points=0.5,
            force_text=0.5
        )
    
    # Legend with black-edged markers
    legend_elements = [
        Line2D([0], [0], color="green", lw=1, label="Bonferroni\nCorrection"),
        Line2D([0], [0], color="red", lw=1, label="Nominal\nSignificance"),
        Line2D([0], [0], marker="^", color="white",
               markerfacecolor="blue", 
               markeredgecolor="black",  # Black edge
               markeredgewidth=0.5,
               markersize=args.legend_size,
               alpha=0.7, 
               label="Increased\nRisk"),
        Line2D([0], [0], marker="v", color="white",
               markerfacecolor="blue",
               markeredgecolor="black",  # Black edge
               markeredgewidth=0.5,
               markersize=args.legend_size,
               alpha=0.3, 
               label="Decreased\nRisk"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize=args.legend_size,
        frameon=True
    )
    
    # Save without tight layout (prevents clipping)
    plt.savefig(args.out, dpi=args.dpi)
    print(f"Plot saved: {args.out}")
    plt.close()


if __name__ == "__main__":
    main()