#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from redox_census_pipeline import (
    REDOX_MODULES,
    TABULA_SAPIENS_COLLECTION_ID,
    add_dataset_filter,
    build_gene_list,
    cpm_log1p,
    dataset_ids_for_collection,
    fetch_subset,
    module_means,
    pseudobulk_scores,
)


def pca1_module_scores(X, var_names, modules):
    var_to_idx = {g: i for i, g in enumerate(var_names)}
    out = {}

    for module, genes in modules.items():
        idx = [var_to_idx[g] for g in genes if g in var_to_idx]
        col = f"score_{module}"
        out[col] = np.full(X.shape[0], np.nan, dtype=float)

        if len(idx) < 2:
            continue

        sub = X[:, idx]
        if hasattr(sub, "toarray"):
            sub = sub.toarray()
        sub = np.asarray(sub, dtype=float)

        means = sub.mean(axis=0, keepdims=True)
        stds = sub.std(axis=0, keepdims=True)
        stds[stds == 0] = 1.0
        Z = (sub - means) / stds

        try:
            U, S, _ = np.linalg.svd(Z, full_matrices=False)
            pc1 = U[:, 0] * S[0]
        except np.linalg.LinAlgError:
            continue

        # Make sign comparable with mean-score direction.
        mean_vec = sub.mean(axis=1)
        if np.corrcoef(pc1, mean_vec)[0, 1] < 0:
            pc1 *= -1

        out[col] = pc1

    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description="Compare mean vs PCA1 redox module scoring")
    ap.add_argument("--census-version", default="2025-11-08")
    ap.add_argument("--obs-filter", default="is_primary_data == True and tissue_general == 'blood'")
    ap.add_argument("--max-cells", type=int, default=12000)
    ap.add_argument("--random-state", type=int, default=33)
    ap.add_argument("--outdir", default="results/scoring_method_experiment")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dataset_ids = dataset_ids_for_collection(args.census_version, TABULA_SAPIENS_COLLECTION_ID)
    obs_filter = add_dataset_filter(args.obs_filter, dataset_ids)

    # Sample joinids manually for reproducibility.
    from redox_census_pipeline import sample_joinids, filter_from_joinids

    joinids = sample_joinids(args.census_version, obs_filter, args.max_cells, args.random_state)
    if joinids is not None:
        obs_filter = filter_from_joinids(obs_filter, joinids)

    genes = build_gene_list(REDOX_MODULES)
    adata = fetch_subset(args.census_version, obs_filter, genes)

    X = cpm_log1p(adata.X)
    var_names = adata.var["feature_name"].tolist()

    mean_scores = module_means(X, var_names, REDOX_MODULES)
    pca_scores = pca1_module_scores(X, var_names, REDOX_MODULES)

    cmp_rows = []
    for c in [c for c in mean_scores.columns if c.startswith("score_")]:
        s1 = mean_scores[c]
        s2 = pca_scores[c]
        ok = s1.notna() & s2.notna()
        if ok.sum() < 5:
            continue
        pearson = float(np.corrcoef(s1[ok], s2[ok])[0, 1])
        spearman = float(pd.Series(s1[ok]).corr(pd.Series(s2[ok]), method="spearman"))
        cmp_rows.append({"module": c, "pearson": pearson, "spearman": spearman, "n": int(ok.sum())})

    cmp = pd.DataFrame(cmp_rows).sort_values("spearman", ascending=False)
    cmp.to_csv(outdir / "mean_vs_pca1_module_correlations.csv", index=False)

    # Plot module-level correlation comparison.
    if not cmp.empty:
        plt.figure(figsize=(11, 5))
        x = np.arange(cmp.shape[0])
        plt.bar(x - 0.2, cmp["pearson"], width=0.4, label="Pearson")
        plt.bar(x + 0.2, cmp["spearman"], width=0.4, label="Spearman")
        plt.xticks(x, cmp["module"], rotation=60, ha="right")
        plt.ylim(-1, 1)
        plt.ylabel("Correlation")
        plt.title("Module Score Agreement: Mean vs PCA1")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "mean_vs_pca1_module_correlations.png", dpi=180)
        plt.close()

    # Build pseudobulk outputs for each method.
    base_obs = adata.obs.reset_index(drop=True)
    mean_df = pd.concat([base_obs, mean_scores], axis=1)
    pca_df = pd.concat([base_obs, pca_scores], axis=1)

    mean_pb = pseudobulk_scores(mean_df, ["donor_id", "cell_type", "tissue_general"], min_cells=20)
    pca_pb = pseudobulk_scores(pca_df, ["donor_id", "cell_type", "tissue_general"], min_cells=20)

    mean_pb.to_parquet(outdir / "mean_method_pseudobulk.parquet", index=False)
    pca_pb.to_parquet(outdir / "pca1_method_pseudobulk.parquet", index=False)

    print(f"Wrote scoring-method experiment artifacts to {outdir}")


if __name__ == "__main__":
    main()
