#!/usr/bin/env python3
"""
Redox PCA Scatter Plots (improved)
===================================
Addresses the issues with the original 2-axis ROS/AOX plot:

1. All 12 modules are used via PCA (not manually split into 2 groups)
2. PC1 = general redox metabolic activity, PC2 = NOX/innate-immune axis
3. Donor-corrected version removes per-donor mean shifts
4. Cell-type-averaged version shows clean cluster positions

Outputs:
  - PCA scatter, all rows (colored by cell type)
  - PCA scatter, donor-corrected (same, tighter clusters)
  - PCA scatter, cell-type means (one dot per cell type, labeled)
  - PCA scatter, colored by tissue
  - PCA loading biplot overlay
  - Quadrant summary CSV
"""

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

TOP_N = 15  # number of cell types to highlight

# Broad cell-lineage grouping for cleaner coloring
LINEAGE_MAP = {
    "macrophage": "Myeloid",
    "monocyte": "Myeloid",
    "classical monocyte": "Myeloid",
    "intermediate monocyte": "Myeloid",
    "non-classical monocyte": "Myeloid",
    "neutrophil": "Myeloid",
    "myeloid dendritic cell": "Myeloid",
    "mast cell": "Myeloid",
    "B cell": "Lymphoid",
    "plasma cell": "Lymphoid",
    "CD4-positive, alpha-beta T cell": "Lymphoid",
    "CD8-positive, alpha-beta T cell": "Lymphoid",
    "T cell": "Lymphoid",
    "natural killer cell": "Lymphoid",
    "mature NK T cell": "Lymphoid",
    "regulatory T cell": "Lymphoid",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "Lymphoid",
    "hematopoietic precursor cell": "Lymphoid",
    "erythrocyte": "Erythroid",
    "fibroblast": "Stromal",
    "myofibroblast cell": "Stromal",
    "smooth muscle cell": "Stromal",
    "vascular associated smooth muscle cell": "Stromal",
    "pericyte": "Stromal",
    "endothelial cell": "Endothelial",
    "endothelial cell of artery": "Endothelial",
    "endothelial cell of lymphatic vessel": "Endothelial",
    "capillary endothelial cell": "Endothelial",
    "vein endothelial cell": "Endothelial",
    "basal cell": "Epithelial",
}

LINEAGE_COLORS = {
    "Myeloid": "#E41A1C",
    "Lymphoid": "#377EB8",
    "Erythroid": "#984EA3",
    "Stromal": "#4DAF4A",
    "Endothelial": "#FF7F00",
    "Epithelial": "#A65628",
    "Other": "#BBBBBB",
}


def load_data(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def get_score_cols(df: pd.DataFrame) -> list:
    return [c for c in df.columns if c.startswith("score_")]


def run_pca(df: pd.DataFrame, score_cols: list, suffix: str = ""):
    """Z-score → PCA. Returns (pca_model, PC_dataframe, scaler)."""
    scaler = StandardScaler()
    X = scaler.fit_transform(df[score_cols])
    pca = PCA(n_components=min(6, len(score_cols)))
    pcs = pca.fit_transform(X)
    pc_df = pd.DataFrame(
        pcs,
        columns=[f"PC{i+1}{suffix}" for i in range(pcs.shape[1])],
        index=df.index,
    )
    return pca, pc_df, scaler


def donor_correct(df: pd.DataFrame, score_cols: list) -> pd.DataFrame:
    """Subtract per-donor mean from each module score."""
    df2 = df.copy()
    donor_means = df2.groupby("donor_id")[score_cols].transform("mean")
    grand_means = df2[score_cols].mean()
    df2[score_cols] = df2[score_cols] - donor_means + grand_means
    return df2


def assign_lineage(df: pd.DataFrame) -> pd.Series:
    return df["cell_type"].map(LINEAGE_MAP).fillna("Other")


def assign_ct_group(df: pd.DataFrame, top_types: list) -> pd.Series:
    return df["cell_type"].astype(str).where(
        df["cell_type"].isin(top_types), "Other"
    )


def make_ct_palette(top_types: list) -> dict:
    tab20 = plt.cm.tab20.colors
    palette = {}
    for i, ct in enumerate(top_types):
        palette[ct] = tab20[i % 20]
    palette["Other"] = (0.82, 0.82, 0.82, 0.35)
    return palette


# ── Plot functions ──────────────────────────────────────────────────

def plot_pca_scatter(df, pc_df, pca, palette, group_col, top_groups,
                     outdir, tag, title_extra="",
                     pc_x="PC1", pc_y="PC2", label_points=False):
    """Generic PCA scatter: color by group_col."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Background: "Other" group
    mask_other = df[group_col] == "Other"
    if mask_other.any():
        ax.scatter(
            pc_df.loc[mask_other, pc_x], pc_df.loc[mask_other, pc_y],
            c=[palette.get("Other", (0.82, 0.82, 0.82, 0.35))],
            s=12, alpha=0.20, edgecolors="none", zorder=1,
        )

    # Highlighted groups
    for grp in top_groups:
        mask = df[group_col] == grp
        if not mask.any():
            continue
        ax.scatter(
            pc_df.loc[mask, pc_x], pc_df.loc[mask, pc_y],
            c=[palette[grp]], s=28, alpha=0.75,
            edgecolors="white", linewidths=0.3,
            label=grp, zorder=2,
        )
        if label_points:
            cx = pc_df.loc[mask, pc_x].mean()
            cy = pc_df.loc[mask, pc_y].mean()
            ax.annotate(
                grp, (cx, cy), fontsize=6.5, alpha=0.85, fontweight="bold",
                textcoords="offset points", xytext=(6, 4),
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
            )

    ev = pca.explained_variance_ratio_
    pc_x_idx = int(pc_x.replace("PC", "").rstrip("_dc")) - 1
    pc_y_idx = int(pc_y.replace("PC", "").rstrip("_dc")) - 1
    ax.set_xlabel(f"{pc_x} ({ev[pc_x_idx]:.1%} variance)", fontsize=12)
    ax.set_ylabel(f"{pc_y} ({ev[pc_y_idx]:.1%} variance)", fontsize=12)
    ax.set_title(f"Redox Module PCA — {title_extra}\n"
                 f"Tabula Sapiens pseudobulk ({len(df)} observations)", fontsize=13)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, framealpha=0.9, title=group_col.replace("_", " ").title(),
              title_fontsize=9)
    ax.axhline(0, color="k", lw=0.4, alpha=0.3)
    ax.axvline(0, color="k", lw=0.4, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"{tag}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {tag}.png / .pdf")


def plot_celltype_mean_pca(df, pc_df, pca, palette, top_types, outdir, tag,
                           title_extra="", pc_x="PC1", pc_y="PC2"):
    """Average PCs per cell type, one dot per cell type."""
    merged = pc_df[[pc_x, pc_y]].copy()
    merged["cell_type"] = df["cell_type"].values
    ct_mean = merged.groupby("cell_type").mean().reset_index()
    ct_mean["n"] = ct_mean["cell_type"].map(df["cell_type"].value_counts())
    ct_mean["ct_group"] = ct_mean["cell_type"].astype(str).where(
        ct_mean["cell_type"].isin(top_types), "Other"
    )

    fig, ax = plt.subplots(figsize=(13, 10))

    other = ct_mean[ct_mean["ct_group"] == "Other"]
    ax.scatter(
        other[pc_x], other[pc_y],
        c=[(0.82, 0.82, 0.82, 0.35)], s=15, alpha=0.25, edgecolors="none", zorder=1,
    )

    for ct in top_types:
        row = ct_mean[ct_mean["cell_type"] == ct]
        if row.empty:
            continue
        ax.scatter(
            row[pc_x], row[pc_y],
            c=[palette[ct]], s=90, alpha=0.90, edgecolors="black", linewidths=0.5,
            label=ct, zorder=3,
        )
        ax.annotate(
            ct, (row[pc_x].values[0], row[pc_y].values[0]),
            fontsize=6.5, alpha=0.85, fontweight="bold",
            textcoords="offset points", xytext=(7, 5),
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
        )

    ev = pca.explained_variance_ratio_
    pc_x_idx = int(pc_x.replace("PC", "").rstrip("_dc")) - 1
    pc_y_idx = int(pc_y.replace("PC", "").rstrip("_dc")) - 1
    ax.set_xlabel(f"{pc_x} ({ev[pc_x_idx]:.1%} variance)", fontsize=12)
    ax.set_ylabel(f"{pc_y} ({ev[pc_y_idx]:.1%} variance)", fontsize=12)
    ax.set_title(f"Redox PCA — Cell Type Means — {title_extra}\n"
                 f"({ct_mean['cell_type'].nunique()} cell types, top {TOP_N} labeled)",
                 fontsize=13)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, framealpha=0.9, title="Cell Type", title_fontsize=9)
    ax.axhline(0, color="k", lw=0.4, alpha=0.3)
    ax.axvline(0, color="k", lw=0.4, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"{tag}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {tag}.png / .pdf")


def plot_biplot(pca, score_cols, pc_df, outdir, tag, pc_x="PC1", pc_y="PC2"):
    """Loading biplot showing how each module contributes to PC1/PC2."""
    fig, ax = plt.subplots(figsize=(10, 8))

    pc_x_idx = int(pc_x.replace("PC", "").rstrip("_dc")) - 1
    pc_y_idx = int(pc_y.replace("PC", "").rstrip("_dc")) - 1

    # Scatter the PC scores in light gray for context
    ax.scatter(pc_df[pc_x], pc_df[pc_y], c="lightgray", s=5, alpha=0.15, zorder=1)

    # Scale arrows so they're visible relative to data spread
    spread = max(pc_df[pc_x].std(), pc_df[pc_y].std()) * 2.5
    loadings = pca.components_

    colors_12 = plt.cm.Set3.colors[:12]
    for i, col in enumerate(score_cols):
        name = col.replace("score_", "")
        lx = loadings[pc_x_idx, i] * spread
        ly = loadings[pc_y_idx, i] * spread
        ax.annotate(
            "", xy=(lx, ly), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=colors_12[i % 12], lw=2),
        )
        ax.text(
            lx * 1.12, ly * 1.12, name, fontsize=7.5, color=colors_12[i % 12],
            fontweight="bold", ha="center", va="center",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ev = pca.explained_variance_ratio_
    ax.set_xlabel(f"{pc_x} ({ev[pc_x_idx]:.1%})", fontsize=12)
    ax.set_ylabel(f"{pc_y} ({ev[pc_y_idx]:.1%})", fontsize=12)
    ax.set_title("PCA Loading Biplot — all 12 redox modules", fontsize=13)
    ax.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax.axvline(0, color="k", lw=0.5, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"{tag}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {tag}.png / .pdf")


def plot_lineage_pca(df, pc_df, pca, outdir, tag, pc_x="PC1", pc_y="PC2",
                     title_extra=""):
    """Color by broad cell lineage for cleaner grouping."""
    lineage = assign_lineage(df)
    df_copy = df.copy()
    df_copy["lineage"] = lineage

    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw each lineage
    lineage_order = ["Myeloid", "Lymphoid", "Erythroid", "Stromal",
                     "Endothelial", "Epithelial", "Other"]
    for lin in lineage_order:
        mask = df_copy["lineage"] == lin
        if not mask.any():
            continue
        color = LINEAGE_COLORS[lin]
        alpha = 0.20 if lin == "Other" else 0.70
        s = 12 if lin == "Other" else 30
        ax.scatter(
            pc_df.loc[mask, pc_x], pc_df.loc[mask, pc_y],
            c=color, s=s, alpha=alpha, edgecolors="white", linewidths=0.2,
            label=f"{lin} (n={mask.sum()})", zorder=1 if lin == "Other" else 2,
        )

    ev = pca.explained_variance_ratio_
    pc_x_idx = int(pc_x.replace("PC", "").rstrip("_dc")) - 1
    pc_y_idx = int(pc_y.replace("PC", "").rstrip("_dc")) - 1
    ax.set_xlabel(f"{pc_x} ({ev[pc_x_idx]:.1%} variance)", fontsize=12)
    ax.set_ylabel(f"{pc_y} ({ev[pc_y_idx]:.1%} variance)", fontsize=12)
    ax.set_title(f"Redox PCA by Cell Lineage — {title_extra}", fontsize=13)
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, framealpha=0.9, title="Lineage", title_fontsize=10)
    ax.axhline(0, color="k", lw=0.4, alpha=0.3)
    ax.axvline(0, color="k", lw=0.4, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"{tag}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {tag}.png / .pdf")


def plot_donor_effect(df, pc_df_raw, pc_df_dc, pca_raw, pca_dc, outdir,
                      top_types, palette):
    """Side-by-side: raw vs donor-corrected PCA for same cell types."""
    fig, axes = plt.subplots(1, 2, figsize=(22, 9))
    ev_r = pca_raw.explained_variance_ratio_
    ev_d = pca_dc.explained_variance_ratio_

    for ax_idx, (ax, pcdf, pca_m, label) in enumerate(zip(
        axes,
        [pc_df_raw, pc_df_dc],
        [pca_raw, pca_dc],
        ["Raw", "Donor-corrected"],
    )):
        ev = pca_m.explained_variance_ratio_
        pcx = pcdf.columns[0]
        pcy = pcdf.columns[1]

        mask_other = df["ct_group"] == "Other"
        if mask_other.any():
            ax.scatter(
                pcdf.loc[mask_other, pcx], pcdf.loc[mask_other, pcy],
                c=[(0.82, 0.82, 0.82, 0.25)], s=10, alpha=0.15,
                edgecolors="none", zorder=1,
            )
        for ct in top_types:
            mask = df["ct_group"] == ct
            if not mask.any():
                continue
            ax.scatter(
                pcdf.loc[mask, pcx], pcdf.loc[mask, pcy],
                c=[palette[ct]], s=24, alpha=0.70,
                edgecolors="white", linewidths=0.2,
                label=ct if ax_idx == 0 else None, zorder=2,
            )

        ax.set_xlabel(f"PC1 ({ev[0]:.1%})", fontsize=11)
        ax.set_ylabel(f"PC2 ({ev[1]:.1%})", fontsize=11)
        ax.set_title(f"{label}", fontsize=13)
        ax.axhline(0, color="k", lw=0.4, alpha=0.3)
        ax.axvline(0, color="k", lw=0.4, alpha=0.3)

    axes[0].legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 1),
                   frameon=True, framealpha=0.9, ncol=1, title="Cell Type",
                   title_fontsize=8)
    fig.suptitle("Donor Effect: Raw vs Donor-Corrected PCA\n"
                 "(same cell types should cluster tighter after correction)",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "pca_raw_vs_donor_corrected.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(outdir / "pca_raw_vs_donor_corrected.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Saved pca_raw_vs_donor_corrected.png / .pdf")


def compute_cluster_tightness(df, pc_df, group_col, pc_x, pc_y):
    """Measure mean within-group spread in PCA space."""
    merged = pd.DataFrame({
        pc_x: pc_df[pc_x].values,
        pc_y: pc_df[pc_y].values,
        group_col: df[group_col].values,
    })
    stats = []
    for grp, sub in merged.groupby(group_col):
        if len(sub) < 3:
            continue
        spread = np.sqrt(sub[pc_x].var() + sub[pc_y].var())
        stats.append({"group": grp, "n": len(sub), "spread": spread})
    return pd.DataFrame(stats).sort_values("spread")


def main():
    parser = argparse.ArgumentParser(description="Improved Redox PCA scatter plots")
    parser.add_argument(
        "--input", default="results/tabulasapiens_combined_pseudobulk_full.parquet",
    )
    parser.add_argument("--outdir", default="results/redox_pca_plots")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    df = load_data(args.input)
    score_cols = get_score_cols(df)
    print(f"  {len(df)} rows, {df['cell_type'].nunique()} cell types, "
          f"{df['donor_id'].nunique()} donors, {len(score_cols)} modules")

    # ── 1. Raw PCA ──────────────────────────────────────────────────
    print("\n── Raw PCA (all 12 modules) ──")
    pca_raw, pc_raw, _ = run_pca(df, score_cols)
    for i, ev in enumerate(pca_raw.explained_variance_ratio_[:4]):
        print(f"  PC{i+1}: {ev:.1%}")

    top_types = df["cell_type"].value_counts().head(TOP_N).index.tolist()
    df["ct_group"] = assign_ct_group(df, top_types)
    palette = make_ct_palette(top_types)

    print("\nPlot 1: Raw PCA — colored by cell type...")
    plot_pca_scatter(df, pc_raw, pca_raw, palette, "ct_group", top_types,
                     outdir, "pca_raw_celltype", "Raw (all 12 modules)")

    print("Plot 2: Raw PCA — cell type means...")
    plot_celltype_mean_pca(df, pc_raw, pca_raw, palette, top_types,
                           outdir, "pca_raw_celltype_means", "Raw")

    print("Plot 3: Raw PCA — by cell lineage...")
    plot_lineage_pca(df, pc_raw, pca_raw, outdir, "pca_raw_lineage", title_extra="Raw")

    print("Plot 4: Loading biplot...")
    plot_biplot(pca_raw, score_cols, pc_raw, outdir, "pca_biplot")

    # ── 2. Donor-corrected PCA ──────────────────────────────────────
    print("\n── Donor-corrected PCA ──")
    df_dc = donor_correct(df, score_cols)
    pca_dc, pc_dc, _ = run_pca(df_dc, score_cols, suffix="_dc")
    for i, ev in enumerate(pca_dc.explained_variance_ratio_[:4]):
        print(f"  PC{i+1}: {ev:.1%}")

    df_dc["ct_group"] = df["ct_group"]  # reuse same groups

    print("\nPlot 5: Donor-corrected PCA — colored by cell type...")
    plot_pca_scatter(df_dc, pc_dc, pca_dc, palette, "ct_group", top_types,
                     outdir, "pca_dc_celltype", "Donor-corrected",
                     pc_x="PC1_dc", pc_y="PC2_dc")

    print("Plot 6: Donor-corrected PCA — cell type means...")
    plot_celltype_mean_pca(df_dc, pc_dc, pca_dc, palette, top_types,
                           outdir, "pca_dc_celltype_means", "Donor-corrected",
                           pc_x="PC1_dc", pc_y="PC2_dc")

    print("Plot 7: Donor-corrected PCA — by lineage...")
    plot_lineage_pca(df_dc, pc_dc, pca_dc, outdir, "pca_dc_lineage",
                     pc_x="PC1_dc", pc_y="PC2_dc", title_extra="Donor-corrected")

    # ── 3. Side-by-side comparison ──────────────────────────────────
    print("\nPlot 8: Raw vs Donor-corrected side-by-side...")
    plot_donor_effect(df, pc_raw, pc_dc, pca_raw, pca_dc, outdir,
                      top_types, palette)

    # ── 4. Quantify donor effect on cluster tightness ───────────────
    print("\n── Cluster tightness (within-cell-type spread in PC space) ──")
    tight_raw = compute_cluster_tightness(df, pc_raw, "cell_type", "PC1", "PC2")
    tight_dc = compute_cluster_tightness(df_dc, pc_dc, "cell_type", "PC1_dc", "PC2_dc")
    merged_tight = tight_raw.merge(tight_dc, on="group", suffixes=("_raw", "_dc"))
    merged_tight["spread_reduction"] = (
        1 - merged_tight["spread_dc"] / merged_tight["spread_raw"]
    )
    merged_tight = merged_tight.sort_values("spread_reduction", ascending=False)

    freq = df["cell_type"].value_counts()
    top30 = set(freq.head(30).index)
    top_tight = merged_tight[merged_tight["group"].isin(top30)]

    print(f"\n  Mean spread reduction (top-30 cell types): "
          f"{top_tight['spread_reduction'].mean():.1%}")
    print(f"  Mean spread reduction (all): "
          f"{merged_tight['spread_reduction'].mean():.1%}")
    print("\n  Top cell types by spread reduction after donor correction:")
    for _, r in top_tight.head(10).iterrows():
        print(f"    {r['group']:50s}  raw={r['spread_raw']:.2f}  "
              f"dc={r['spread_dc']:.2f}  reduction={r['spread_reduction']:+.1%}")

    merged_tight.to_csv(outdir / "cluster_tightness.csv", index=False)
    print(f"\n  Saved cluster_tightness.csv ({len(merged_tight)} cell types)")

    # Save PCA coordinates + metadata
    export = df[["donor_id", "cell_type", "tissue_general", "n_cells"]].copy()
    export["PC1_raw"] = pc_raw["PC1"].values
    export["PC2_raw"] = pc_raw["PC2"].values
    export["PC1_dc"] = pc_dc["PC1_dc"].values
    export["PC2_dc"] = pc_dc["PC2_dc"].values
    export["lineage"] = assign_lineage(df)
    export.to_csv(outdir / "pca_coordinates.csv", index=False)
    print(f"  Saved pca_coordinates.csv")

    print("\nDone! All plots in", outdir)


if __name__ == "__main__":
    main()
