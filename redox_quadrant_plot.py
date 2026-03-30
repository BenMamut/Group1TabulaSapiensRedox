#!/usr/bin/env python3
"""
Redox Quadrant Scatter Plot
============================
Plots cell-type pseudobulk observations on a 2D plane:
  X-axis  →  ROS Production  (z-scored NOX_SYSTEM + MITO_ROS average)
  Y-axis  →  Antioxidant Defense  (z-scored SOD + CAT/PRXSM + GSH + TXN average)

Points are colored by cell type (top N highlighted, rest gray).
Quadrant lines drawn at the median of each composite axis.
"""

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

# ── Module groupings ────────────────────────────────────────────────
ROS_MODULES = ["score_NOX_SYSTEM", "score_MITO_ROS"]
ANTIOXIDANT_MODULES = [
    "score_SOD_SYSTEM",
    "score_CATALASE_PEROXISOME",
    "score_GLUTATHIONE_SYSTEM",
    "score_THIOREDOXIN_SYSTEM",
]

TOP_N = 15  # number of cell types to label individually


def load_and_score(parquet_path: str) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)

    # Z-score each module across all rows
    for col in ROS_MODULES + ANTIOXIDANT_MODULES:
        m, s = df[col].mean(), df[col].std()
        df[f"z_{col}"] = (df[col] - m) / s

    # Composite scores (mean of z-scored components)
    df["ROS_production"] = df[[f"z_{c}" for c in ROS_MODULES]].mean(axis=1)
    df["Antioxidant_defense"] = df[[f"z_{c}" for c in ANTIOXIDANT_MODULES]].mean(axis=1)

    return df


def assign_colors(df: pd.DataFrame):
    """Return (color_col, palette_dict, top_types)."""
    top_types = df["cell_type"].value_counts().head(TOP_N).index.tolist()
    df["ct_group"] = df["cell_type"].astype(str).where(df["cell_type"].isin(top_types), "Other")

    # Qualitative palette
    tab20 = plt.cm.tab20.colors
    palette = {}
    for i, ct in enumerate(top_types):
        palette[ct] = tab20[i % 20]
    palette["Other"] = (0.82, 0.82, 0.82, 0.35)

    return palette, top_types


def plot_quadrant_scatter(df, palette, top_types, outdir, tag="all_rows"):
    """Main scatter: one dot per pseudobulk row."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw "Other" first (background)
    other = df[df["ct_group"] == "Other"]
    ax.scatter(
        other["ROS_production"], other["Antioxidant_defense"],
        c=[palette["Other"]], s=18, alpha=0.30, edgecolors="none", label=None, zorder=1,
    )

    # Highlighted cell types on top
    for ct in top_types:
        sub = df[df["ct_group"] == ct]
        ax.scatter(
            sub["ROS_production"], sub["Antioxidant_defense"],
            c=[palette[ct]], s=32, alpha=0.80, edgecolors="white", linewidths=0.3,
            label=ct, zorder=2,
        )

    # Quadrant lines at median
    med_ros = df["ROS_production"].median()
    med_aox = df["Antioxidant_defense"].median()
    ax.axvline(med_ros, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(med_aox, color="k", ls="--", lw=0.8, alpha=0.5)

    # Quadrant labels
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    txt_kw = dict(fontsize=9, fontstyle="italic", alpha=0.6,
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.text(xlim[1] * 0.85, ylim[1] * 0.90, "High ROS\nHigh Antioxidant",
            ha="center", va="center", color="darkgreen", **txt_kw)
    ax.text(xlim[0] * 0.85, ylim[1] * 0.90, "Low ROS\nHigh Antioxidant",
            ha="center", va="center", color="steelblue", **txt_kw)
    ax.text(xlim[1] * 0.85, ylim[0] * 0.85, "High ROS\nLow Antioxidant",
            ha="center", va="center", color="firebrick", **txt_kw)
    ax.text(xlim[0] * 0.85, ylim[0] * 0.85, "Low ROS\nLow Antioxidant",
            ha="center", va="center", color="gray", **txt_kw)

    ax.set_xlabel("ROS Production  (z-scored NOX + MITO_ROS)", fontsize=12)
    ax.set_ylabel("Antioxidant Defense  (z-scored SOD + CAT + GSH + TXN)", fontsize=12)
    ax.set_title("Redox Balance: ROS Production vs Antioxidant Defense\n"
                 "(Tabula Sapiens pseudobulk, colored by cell type)", fontsize=13)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, framealpha=0.9, title="Cell Type", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / f"redox_quadrant_{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"redox_quadrant_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved redox_quadrant_{tag}.png / .pdf")


def plot_celltype_means(df, palette, top_types, outdir):
    """One dot per cell type (averaged across donors/tissues), with labels."""
    ct_mean = (
        df.groupby("cell_type")[["ROS_production", "Antioxidant_defense"]]
        .mean()
        .reset_index()
    )
    ct_mean["ct_group"] = ct_mean["cell_type"].astype(str).where(
        ct_mean["cell_type"].isin(top_types), "Other"
    )
    ct_mean["n"] = ct_mean["cell_type"].map(df["cell_type"].value_counts())

    fig, ax = plt.subplots(figsize=(12, 10))

    other = ct_mean[ct_mean["ct_group"] == "Other"]
    ax.scatter(
        other["ROS_production"], other["Antioxidant_defense"],
        c=[palette["Other"]], s=15, alpha=0.25, edgecolors="none", zorder=1,
    )

    for ct in top_types:
        row = ct_mean[ct_mean["cell_type"] == ct]
        if row.empty:
            continue
        ax.scatter(
            row["ROS_production"], row["Antioxidant_defense"],
            c=[palette[ct]], s=80, alpha=0.90, edgecolors="black", linewidths=0.5,
            label=ct, zorder=3,
        )
        ax.annotate(
            ct, (row["ROS_production"].values[0], row["Antioxidant_defense"].values[0]),
            fontsize=6.5, alpha=0.85, textcoords="offset points", xytext=(6, 4),
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    med_ros = ct_mean["ROS_production"].median()
    med_aox = ct_mean["Antioxidant_defense"].median()
    ax.axvline(med_ros, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(med_aox, color="k", ls="--", lw=0.8, alpha=0.5)

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    txt_kw = dict(fontsize=9, fontstyle="italic", alpha=0.6,
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.text(xlim[1] * 0.85, ylim[1] * 0.90, "High ROS\nHigh Antioxidant",
            ha="center", va="center", color="darkgreen", **txt_kw)
    ax.text(xlim[0] * 0.85, ylim[1] * 0.90, "Low ROS\nHigh Antioxidant",
            ha="center", va="center", color="steelblue", **txt_kw)
    ax.text(xlim[1] * 0.85, ylim[0] * 0.85, "High ROS\nLow Antioxidant",
            ha="center", va="center", color="firebrick", **txt_kw)
    ax.text(xlim[0] * 0.85, ylim[0] * 0.85, "Low ROS\nLow Antioxidant",
            ha="center", va="center", color="gray", **txt_kw)

    ax.set_xlabel("ROS Production  (z-scored NOX + MITO_ROS)", fontsize=12)
    ax.set_ylabel("Antioxidant Defense  (z-scored SOD + CAT + GSH + TXN)", fontsize=12)
    ax.set_title("Redox Balance by Cell Type (mean across donors & tissues)\n"
                 "Tabula Sapiens, top 15 cell types labeled", fontsize=13)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, framealpha=0.9, title="Cell Type", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "redox_quadrant_celltype_means.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / "redox_quadrant_celltype_means.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Saved redox_quadrant_celltype_means.png / .pdf")


def print_quadrant_summary(df: pd.DataFrame):
    """Print which cell types land in each quadrant."""
    med_ros = df["ROS_production"].median()
    med_aox = df["Antioxidant_defense"].median()

    ct_mean = (
        df.groupby("cell_type")[["ROS_production", "Antioxidant_defense"]]
        .mean()
        .reset_index()
    )

    quadrants = {
        "High ROS / High Antioxidant": ct_mean[
            (ct_mean["ROS_production"] > med_ros) & (ct_mean["Antioxidant_defense"] > med_aox)
        ],
        "High ROS / Low Antioxidant": ct_mean[
            (ct_mean["ROS_production"] > med_ros) & (ct_mean["Antioxidant_defense"] <= med_aox)
        ],
        "Low ROS / High Antioxidant": ct_mean[
            (ct_mean["ROS_production"] <= med_ros) & (ct_mean["Antioxidant_defense"] > med_aox)
        ],
        "Low ROS / Low Antioxidant": ct_mean[
            (ct_mean["ROS_production"] <= med_ros) & (ct_mean["Antioxidant_defense"] <= med_aox)
        ],
    }

    freq = df["cell_type"].value_counts()
    top30 = set(freq.head(30).index)

    print("\n── Quadrant Summary (cell types with ≥5 observations, top-30 highlighted) ──")
    for qname, qdf in quadrants.items():
        qdf_filt = qdf[qdf["cell_type"].isin(freq[freq >= 5].index)]
        qdf_top = qdf_filt[qdf_filt["cell_type"].isin(top30)].sort_values(
            "ROS_production", ascending=False
        )
        print(f"\n  {qname}  ({len(qdf)} total, {len(qdf_filt)} with n≥5)")
        if not qdf_top.empty:
            for _, r in qdf_top.iterrows():
                print(f"    • {r['cell_type']:50s}  ROS={r['ROS_production']:+.2f}  AOX={r['Antioxidant_defense']:+.2f}")


def main():
    parser = argparse.ArgumentParser(description="Redox ROS vs Antioxidant quadrant plot")
    parser.add_argument(
        "--input", default="results/tabulasapiens_combined_pseudobulk_full.parquet",
        help="Path to combined pseudobulk parquet",
    )
    parser.add_argument(
        "--outdir", default="results/redox_quadrant_plots",
        help="Output directory",
    )
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading pseudobulk data...")
    df = load_and_score(args.input)
    print(f"  {len(df)} rows, {df['cell_type'].nunique()} cell types")

    palette, top_types = assign_colors(df)

    print("\nPlot 1: All pseudobulk rows scatter...")
    plot_quadrant_scatter(df, palette, top_types, outdir)

    print("Plot 2: Cell-type means scatter...")
    plot_celltype_means(df, palette, top_types, outdir)

    print_quadrant_summary(df)

    # Save quadrant assignments
    ct_mean = (
        df.groupby("cell_type")[["ROS_production", "Antioxidant_defense"]]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    ct_mean.columns = ["cell_type", "ROS_mean", "ROS_std", "ROS_n",
                        "AOX_mean", "AOX_std", "AOX_n"]
    med_ros = df["ROS_production"].median()
    med_aox = df["Antioxidant_defense"].median()
    ct_mean["quadrant"] = np.where(
        ct_mean["ROS_mean"] > med_ros,
        np.where(ct_mean["AOX_mean"] > med_aox, "High ROS / High AOX", "High ROS / Low AOX"),
        np.where(ct_mean["AOX_mean"] > med_aox, "Low ROS / High AOX", "Low ROS / Low AOX"),
    )
    ct_mean = ct_mean.sort_values("ROS_mean", ascending=False)
    ct_mean.to_csv(outdir / "quadrant_assignments.csv", index=False)
    print(f"\n  Saved quadrant_assignments.csv ({len(ct_mean)} cell types)")

    print("\nDone!")


if __name__ == "__main__":
    main()
