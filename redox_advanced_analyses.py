#!/usr/bin/env python3
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import spearmanr

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
    sample_joinids,
    filter_from_joinids,
    variance_partition,
)


def _score_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c.startswith("score_") and "_z_in_" not in c]


def donor_coherence(df: pd.DataFrame, min_shared_donors: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_cols = _score_columns(df)
    pair_rows = []

    for score_col in score_cols:
        for cell_type, sub_ct in df.groupby("cell_type", dropna=False):
            pivot = sub_ct.pivot_table(index="donor_id", columns="tissue_general", values=score_col, aggfunc="mean")
            tissues = [t for t in pivot.columns.tolist() if pd.notna(t)]

            for t1, t2 in combinations(tissues, 2):
                pair = pivot[[t1, t2]].dropna()
                n = pair.shape[0]
                if n < min_shared_donors:
                    continue
                corr = spearmanr(pair[t1], pair[t2]).correlation
                pair_rows.append(
                    {
                        "module": score_col,
                        "cell_type": cell_type,
                        "tissue_a": t1,
                        "tissue_b": t2,
                        "n_shared_donors": int(n),
                        "spearman": float(corr) if corr == corr else np.nan,
                    }
                )

    pairs_df = pd.DataFrame(pair_rows)
    if pairs_df.empty:
        return pairs_df, pd.DataFrame()

    summary = (
        pairs_df.groupby("module", dropna=False)
        .agg(
            mean_spearman=("spearman", "mean"),
            median_spearman=("spearman", "median"),
            n_pairs=("spearman", "size"),
            mean_shared_donors=("n_shared_donors", "mean"),
        )
        .reset_index()
        .sort_values("mean_spearman", ascending=False)
    )
    return pairs_df, summary


def aucell_style_scores(
    X: sparse.spmatrix | np.ndarray,
    var_names: List[str],
    modules: Dict[str, List[str]],
    max_rank: int,
) -> pd.DataFrame:
    """
    AUCell-like score on the fetched gene universe.
    This ranks genes within each cell and measures module recovery in top ranks.
    """
    if sparse.issparse(X):
        A = X.toarray()
    else:
        A = np.asarray(X)

    n_cells, n_genes = A.shape
    max_rank = min(max_rank, n_genes)

    # ranks is 1-based rank per gene per cell (1 means highest expression in that cell).
    order = np.argsort(-A, axis=1)
    ranks = np.empty_like(order)
    row_idx = np.arange(n_cells)[:, None]
    ranks[row_idx, order] = np.arange(1, n_genes + 1)

    var_to_idx = {g: i for i, g in enumerate(var_names)}
    out = {}
    for module, genes in modules.items():
        idx = [var_to_idx[g] for g in genes if g in var_to_idx]
        col = f"score_{module}"
        if not idx:
            out[col] = np.full(n_cells, np.nan, dtype=float)
            continue

        r = ranks[:, idx]
        contrib = np.clip(max_rank - r + 1, 0, None)
        score = contrib.sum(axis=1) / (len(idx) * max_rank)
        out[col] = score.astype(float)

    return pd.DataFrame(out)


def background_adjusted_scores(
    X: sparse.spmatrix | np.ndarray,
    var_names: List[str],
    modules: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Module mean minus per-cell global background over the fetched redox-gene universe.
    """
    mean_df = module_means(X, var_names, modules)
    if sparse.issparse(X):
        bg = np.asarray(X.mean(axis=1)).ravel()
    else:
        bg = np.asarray(X.mean(axis=1)).ravel()

    out = pd.DataFrame(index=mean_df.index)
    for c in mean_df.columns:
        if c.startswith("score_"):
            out[c] = mean_df[c] - bg
    return out


def _method_agreement_table(method_scores: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    methods = list(method_scores.keys())
    rows = []

    modules = _score_columns(next(iter(method_scores.values())))
    for module in modules:
        for i, m1 in enumerate(methods):
            for m2 in methods[i + 1 :]:
                s1 = method_scores[m1][module]
                s2 = method_scores[m2][module]
                ok = s1.notna() & s2.notna()
                if ok.sum() < 5:
                    continue
                rows.append(
                    {
                        "module": module,
                        "method_a": m1,
                        "method_b": m2,
                        "pearson": float(np.corrcoef(s1[ok], s2[ok])[0, 1]),
                        "spearman": float(pd.Series(s1[ok]).corr(pd.Series(s2[ok]), method="spearman")),
                        "n": int(ok.sum()),
                    }
                )
    return pd.DataFrame(rows)


def run_donor_coherence(args: argparse.Namespace) -> None:
    in_file = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if in_file.suffix.lower() == ".csv":
        df = pd.read_csv(in_file)
    else:
        df = pd.read_parquet(in_file)

    pairs_df, summary = donor_coherence(df, min_shared_donors=args.min_shared_donors)
    pairs_df.to_csv(outdir / "donor_coherence_pairs.csv", index=False)
    summary.to_csv(outdir / "donor_coherence_summary.csv", index=False)

    if not summary.empty:
        plt.figure(figsize=(11, 5))
        s = summary.sort_values("mean_spearman", ascending=False)
        sns.barplot(data=s, x="module", y="mean_spearman", color="#2E86AB")
        plt.xticks(rotation=60, ha="right")
        plt.ylabel("Mean Spearman across tissue pairs")
        plt.title("Cross-Tissue Donor Coherence by Module")
        plt.tight_layout()
        plt.savefig(outdir / "donor_coherence_module_bar.png", dpi=180)
        plt.close()

    if not pairs_df.empty:
        ct = (
            pairs_df.groupby("cell_type", dropna=False)
            .agg(mean_spearman=("spearman", "mean"), n_pairs=("spearman", "size"))
            .sort_values("mean_spearman", ascending=False)
            .head(25)
            .reset_index()
        )
        ct.to_csv(outdir / "donor_coherence_top_celltypes.csv", index=False)
        plt.figure(figsize=(10, 7))
        sns.barplot(data=ct, x="mean_spearman", y="cell_type", color="#4C956C")
        plt.xlabel("Mean Spearman")
        plt.ylabel("Cell type")
        plt.title("Top Cell Types by Cross-Tissue Donor Coherence")
        plt.tight_layout()
        plt.savefig(outdir / "donor_coherence_top_celltypes.png", dpi=180)
        plt.close()

    print(f"Wrote donor coherence artifacts to {outdir}")


def run_scoring_benchmark(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dataset_ids = dataset_ids_for_collection(args.census_version, TABULA_SAPIENS_COLLECTION_ID)
    obs_filter = add_dataset_filter(args.obs_filter, dataset_ids)

    joinids = sample_joinids(args.census_version, obs_filter, args.max_cells, args.random_state)
    if joinids is not None:
        obs_filter = filter_from_joinids(obs_filter, joinids)

    genes = build_gene_list(REDOX_MODULES)
    adata = fetch_subset(args.census_version, obs_filter, genes)

    X = cpm_log1p(adata.X)
    var_names = adata.var["feature_name"].tolist()

    mean_scores = module_means(X, var_names, REDOX_MODULES)
    aucell_scores = aucell_style_scores(X, var_names, REDOX_MODULES, max_rank=args.aucell_max_rank)
    bg_scores = background_adjusted_scores(X, var_names, REDOX_MODULES)

    methods = {
        "mean": mean_scores,
        "aucell_style": aucell_scores,
        "bg_adjusted": bg_scores,
    }

    agreement = _method_agreement_table(methods)
    agreement.to_csv(outdir / "method_agreement.csv", index=False)

    if not agreement.empty:
        agg = (
            agreement.groupby(["method_a", "method_b"], dropna=False)
            .agg(mean_spearman=("spearman", "mean"), mean_pearson=("pearson", "mean"))
            .reset_index()
        )
        agg.to_csv(outdir / "method_agreement_summary.csv", index=False)

        plt.figure(figsize=(9, 5))
        sns.barplot(data=agg, x="method_a", y="mean_spearman", hue="method_b")
        plt.ylim(-1, 1)
        plt.ylabel("Mean module Spearman")
        plt.title("Average Agreement Between Scoring Methods")
        plt.tight_layout()
        plt.savefig(outdir / "method_agreement_summary.png", dpi=180)
        plt.close()

    base_obs = adata.obs.reset_index(drop=True)
    var_rows = []
    for name, score_df in methods.items():
        method_df = pd.concat([base_obs, score_df], axis=1)
        pb = pseudobulk_scores(method_df, ["donor_id", "cell_type", "tissue_general"], min_cells=args.min_cells)
        pb.to_parquet(outdir / f"{name}_pseudobulk.parquet", index=False)

        vp = variance_partition(pb, factors=["cell_type", "tissue_general", "donor_id"])
        vp.to_csv(outdir / f"variance_partition_{name}.csv", index=False)

        if not vp.empty:
            var_rows.append(
                {
                    "method": name,
                    "avg_frac_cell_type": float(vp["frac_cell_type"].mean()),
                    "avg_frac_tissue_general": float(vp["frac_tissue_general"].mean()),
                    "avg_frac_donor_id": float(vp["frac_donor_id"].mean()),
                }
            )

    if var_rows:
        var_summary = pd.DataFrame(var_rows)
        var_summary.to_csv(outdir / "variance_partition_method_summary.csv", index=False)

        plt.figure(figsize=(8, 5))
        x = np.arange(var_summary.shape[0])
        w = 0.25
        plt.bar(x - w, var_summary["avg_frac_cell_type"], width=w, label="cell_type")
        plt.bar(x, var_summary["avg_frac_tissue_general"], width=w, label="tissue")
        plt.bar(x + w, var_summary["avg_frac_donor_id"], width=w, label="donor")
        plt.xticks(x, var_summary["method"])
        plt.ylabel("Average variance fraction")
        plt.title("Variance Partition Sensitivity by Scoring Method")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "variance_partition_method_summary.png", dpi=180)
        plt.close()

    print(f"Wrote scoring benchmark artifacts to {outdir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Advanced redox analyses")
    p.add_argument("--census-version", default="2025-11-08")

    sub = p.add_subparsers(dest="command", required=True)

    p_coh = sub.add_parser("donor-coherence", help="Cross-tissue donor coherence analysis")
    p_coh.add_argument("--input", required=True, help="Combined pseudobulk file (.parquet or .csv)")
    p_coh.add_argument("--min-shared-donors", type=int, default=4)
    p_coh.add_argument("--outdir", required=True)
    p_coh.set_defaults(func=run_donor_coherence)

    p_bench = sub.add_parser("scoring-benchmark", help="Benchmark mean vs AUCell-style vs background-adjusted module scoring")
    p_bench.add_argument("--obs-filter", default="is_primary_data == True and tissue_general == 'blood'")
    p_bench.add_argument("--max-cells", type=int, default=12000)
    p_bench.add_argument("--random-state", type=int, default=33)
    p_bench.add_argument("--aucell-max-rank", type=int, default=50)
    p_bench.add_argument("--min-cells", type=int, default=20)
    p_bench.add_argument("--outdir", required=True)
    p_bench.set_defaults(func=run_scoring_benchmark)

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
