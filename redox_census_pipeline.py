#!/usr/bin/env python3
"""
Memory-safe redox module exploration using CELLxGENE Census.

This script is designed to avoid loading all Tabula Sapiens cells at once.
It lets you:
1) Preview cell metadata counts for a filtered subset.
2) Compute redox module scores on a filtered subset.
3) Run tissue-by-tissue batches and write each result to disk.

Example usage:
    python redox_census_pipeline.py preview \
      --obs-filter "is_primary_data == True" \
      --groupby tissue_general

    python redox_census_pipeline.py score \
      --obs-filter "is_primary_data == True and tissue_general == 'blood'" \
      --out results/blood_scores.parquet

    python redox_census_pipeline.py batch \
      --tissues blood liver lung \
      --outdir results/by_tissue
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse

import cellxgene_census


REDOX_MODULES: Dict[str, List[str]] = {
    "NOX_SYSTEM": ["NOX1", "NOX2", "CYBB", "NOX3", "NOX4", "NOX5", "DUOX1", "DUOX2", "CYBA", "NCF1", "NCF2", "NCF4", "RAC1", "RAC2", "NOXA1", "NOXO1", "POLDIP2"],
    "MITO_ROS": ["NDUFA1", "NDUFA2", "NDUFA9", "NDUFS1", "NDUFS2", "NDUFS3", "SDHA", "SDHB", "SDHC", "SDHD", "UQCRC1", "UQCRC2", "UQCRFS1", "COX4I1", "COX5A", "COX6C", "ATP5F1A", "ATP5F1B", "FDX1", "ETFDH"],
    "SOD_SYSTEM": ["SOD1", "SOD2", "SOD3", "CCS"],
    "CATALASE_PEROXISOME": ["CAT", "PRDX5", "ACOX1", "ACOX2", "ACOX3", "HSD17B4", "EHHADH", "PEX1", "PEX2", "PEX3", "PEX5", "PEX6", "PEX10", "PEX11A", "PEX11B", "PEX13", "PEX14", "SCP2"],
    "GLUTATHIONE_SYSTEM": ["GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX3", "GPX4", "GPX7", "GPX8", "GSTP1", "GSTA1", "GSTA2", "GSTM1", "GSTM2", "GSTM3", "MGST1", "MGST2", "MGST3", "GLRX", "GLRX2", "GLRX3", "TXNDC12"],
    "THIOREDOXIN_SYSTEM": ["TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6", "SRXN1", "SELENOW", "SELENOH", "TXNDC5", "TXNDC9", "TXNDC17"],
    "PPP_MODULE": ["G6PD", "PGLS", "PGD", "RPIA", "TKT", "TALDO1", "RPE", "PRPS1", "PRPS2", "H6PD"],
    "NADPH_DEHYDROGENASES": ["ME1", "ME2", "ME3", "IDH1", "IDH2", "MTHFD1", "MTHFD1L", "MTHFD2", "MTHFD2L", "ALDH1L1", "ALDH1L2", "NNT"],
    "NRF2_TARGETS": ["NFE2L2", "KEAP1", "NQO1", "HMOX1", "GCLC", "GCLM", "TXNRD1", "SRXN1", "PRDX1", "G6PD", "PGD", "ME1", "FTH1", "FTL", "SLC7A11", "ABCC1", "ABCC2", "GSTP1", "MGST1"],
    "FOXO_OXSTRESS": ["FOXO1", "FOXO3", "FOXO4", "SOD2", "CAT", "PRDX3", "SESN1", "SESN2", "SESN3", "GADD45A", "GADD45B", "CDKN1A", "BNIP3", "TXNIP", "LC3B", "ATG5", "ATG7"],
    "OXPHOS_MODULE": ["NDUFS1", "NDUFS2", "NDUFS3", "SDHA", "SDHB", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "COX6B1", "ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5PD", "CYCS", "TFAM"],
    "GLYCOLYSIS_MODULE": ["HK1", "HK2", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB", "SLC2A1", "SLC2A3"],
}

TABULA_SAPIENS_COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"

OBS_COLUMNS = [
    "cell_type",
    "tissue",
    "tissue_general",
    "donor_id",
    "sex",
    "assay",
    "soma_joinid",
]


def build_gene_list(modules: Dict[str, Iterable[str]]) -> List[str]:
    return sorted({g for genes in modules.values() for g in genes})


def census_soma_uri(census_version: str) -> str:
    return f"s3://cellxgene-census-public-us-west-2/cell-census/{census_version}/soma/"


def open_census(census_version: str):
    return cellxgene_census.open_soma(uri=census_soma_uri(census_version))


def cpm_log1p(X: sparse.spmatrix | np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """Compute log1p(CPM-like) normalized expression without densifying full matrix."""
    if sparse.issparse(X):
        X = X.tocsr(copy=True)
        cell_sums = np.asarray(X.sum(axis=1)).ravel()
        scale = np.divide(target_sum, cell_sums, out=np.zeros_like(cell_sums, dtype=float), where=cell_sums > 0)
        X = sparse.diags(scale) @ X
        X.data = np.log1p(X.data)
        return X

    # Dense fallback (small subsets)
    cell_sums = X.sum(axis=1)
    scale = np.divide(target_sum, cell_sums, out=np.zeros_like(cell_sums, dtype=float), where=cell_sums > 0)
    X_norm = X * scale[:, None]
    return np.log1p(X_norm)


def module_means(X: sparse.spmatrix | np.ndarray, var_names: List[str], modules: Dict[str, List[str]]) -> pd.DataFrame:
    var_to_idx = {g: i for i, g in enumerate(var_names)}
    out = {}

    for module_name, genes in modules.items():
        idx = [var_to_idx[g] for g in genes if g in var_to_idx]
        out[f"score_{module_name}"] = np.full(X.shape[0], np.nan, dtype=float)
        out[f"n_genes_{module_name}"] = len(idx)

        if not idx:
            continue

        sub = X[:, idx]
        if sparse.issparse(sub):
            score = np.asarray(sub.mean(axis=1)).ravel()
        else:
            score = sub.mean(axis=1)
        out[f"score_{module_name}"] = score

    return pd.DataFrame(out)


def zscore_within_group(df: pd.DataFrame, score_cols: List[str], group_col: str) -> pd.DataFrame:
    out = df.copy()
    group_stats = out.groupby(group_col)[score_cols].agg(["mean", "std"])

    for col in score_cols:
        means = out[group_col].map(group_stats[(col, "mean")])
        stds = out[group_col].map(group_stats[(col, "std")])
        out[f"{col}_z_in_{group_col}"] = (out[col] - means) / stds.replace(0, np.nan)

    return out


def obs_table(
    census_version: str,
    obs_filter: str,
    columns: List[str],
) -> pd.DataFrame:
    with open_census(census_version) as census:
        obs = (
            census["census_data"]["homo_sapiens"].obs.read(
                value_filter=obs_filter,
                column_names=columns,
            )
        ).concat().to_pandas()
    return obs


def dataset_ids_for_collection(census_version: str, collection_id: str) -> List[str]:
    with open_census(census_version) as census:
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    hits = datasets[datasets["collection_id"] == collection_id]
    return hits["dataset_id"].astype(str).tolist()


def add_dataset_filter(obs_filter: str, dataset_ids: List[str]) -> str:
    if not dataset_ids:
        raise ValueError("Collection filter matched zero datasets.")
    return f"({obs_filter}) and dataset_id in {dataset_ids}"


def filter_from_joinids(base_filter: str, joinids: List[int]) -> str:
    return f"({base_filter}) and soma_joinid in {joinids}"


def sample_joinids(
    census_version: str,
    obs_filter: str,
    max_cells: int | None,
    random_state: int,
) -> List[int] | None:
    if max_cells is None:
        return None

    obs = obs_table(
        census_version=census_version,
        obs_filter=obs_filter,
        columns=["soma_joinid"],
    )
    if obs.empty:
        return []

    n = min(max_cells, obs.shape[0])
    sampled = obs.sample(n=n, random_state=random_state, replace=False)
    return sampled["soma_joinid"].astype(int).tolist()


def fetch_subset(
    census_version: str,
    obs_filter: str,
    genes: List[str],
) -> "anndata.AnnData":
    gene_filter = f"feature_name in {genes}"

    with open_census(census_version) as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            measurement_name="RNA",
            obs_value_filter=obs_filter,
            obs_column_names=OBS_COLUMNS,
            var_value_filter=gene_filter,
        )
    return adata


def preview_cells(census_version: str, obs_filter: str, groupby: str, out_csv: Path | None) -> pd.DataFrame:
    obs = obs_table(census_version=census_version, obs_filter=obs_filter, columns=[groupby])
    counts = obs.groupby(groupby, dropna=False).size().rename("n_cells").sort_values(ascending=False)
    df = counts.reset_index()

    if out_csv is not None:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)

    return df


def distinct_obs_values(census_version: str, obs_filter: str, column: str) -> List[str]:
    df = obs_table(census_version=census_version, obs_filter=obs_filter, columns=[column])
    values = df[column].dropna().astype(str).sort_values().unique().tolist()
    return values


def compute_scores(
    census_version: str,
    obs_filter: str,
    modules: Dict[str, List[str]],
    zscore_group: str | None,
    max_cells: int | None,
    random_state: int,
) -> pd.DataFrame:
    joinids = sample_joinids(
        census_version=census_version,
        obs_filter=obs_filter,
        max_cells=max_cells,
        random_state=random_state,
    )
    if joinids is not None:
        if len(joinids) == 0:
            score_cols = [f"score_{m}" for m in modules]
            n_gene_cols = [f"n_genes_{m}" for m in modules]
            empty = pd.DataFrame(columns=OBS_COLUMNS + score_cols + n_gene_cols)
            if zscore_group is not None:
                for c in score_cols:
                    empty[f"{c}_z_in_{zscore_group}"] = pd.Series(dtype=float)
            return empty
        obs_filter = filter_from_joinids(obs_filter, joinids)

    genes = build_gene_list(modules)
    adata = fetch_subset(census_version=census_version, obs_filter=obs_filter, genes=genes)

    X = cpm_log1p(adata.X)
    var_names = adata.var["feature_name"].tolist()
    scores = module_means(X, var_names, modules)

    result = pd.concat([adata.obs.reset_index(drop=True), scores], axis=1)

    score_cols = [c for c in result.columns if c.startswith("score_")]
    if zscore_group is not None and not result.empty:
        result = zscore_within_group(result, score_cols=score_cols, group_col=zscore_group)

    return result


def pseudobulk_scores(
    df: pd.DataFrame,
    group_cols: List[str],
    min_cells: int,
) -> pd.DataFrame:
    score_cols = [c for c in df.columns if c.startswith("score_") and "_z_in_" not in c]
    if not score_cols:
        raise ValueError("No raw score_ columns found for pseudobulk aggregation.")

    if df.empty:
        columns = group_cols + ["n_cells"] + score_cols
        return pd.DataFrame(columns=columns)

    missing = [c for c in group_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing group columns for pseudobulk: {missing}")

    grouped = df.groupby(group_cols, dropna=False)
    counts = grouped.size().rename("n_cells")
    means = grouped[score_cols].mean()
    out = counts.to_frame().join(means).reset_index()
    out = out[out["n_cells"] >= min_cells].reset_index(drop=True)
    return out


def combine_tables(paths: List[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
        else:
            df = pd.read_parquet(path)
        if df.empty:
            continue
        df = df.copy()
        df["source_file"] = path.name
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def variance_partition(df: pd.DataFrame, factors: List[str]) -> pd.DataFrame:
    score_cols = [c for c in df.columns if c.startswith("score_") and "_z_in_" not in c]
    rows = []

    def sse(y: np.ndarray, X: np.ndarray) -> float:
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - (X @ coef)
        return float(np.sum(resid * resid))

    for score_col in score_cols:
        usable_factors = [factor for factor in factors if factor in df.columns and df[factor].nunique(dropna=False) > 1]
        if len(usable_factors) == 0:
            continue

        model_df = df[[score_col] + usable_factors].dropna().copy()
        if model_df.shape[0] <= len(usable_factors) + 1:
            continue

        y = model_df[score_col].to_numpy(dtype=float)
        y_mean = float(np.mean(y))
        sst = float(np.sum((y - y_mean) ** 2))
        if sst <= 0:
            continue

        n = y.shape[0]
        X_prev = np.ones((n, 1), dtype=float)
        sse_prev = sse(y, X_prev)
        row = {"module": score_col, "n_rows": n}

        for factor in factors:
            if factor not in usable_factors:
                row[f"frac_{factor}"] = 0.0
                continue

            dummies = pd.get_dummies(model_df[factor].astype(str), drop_first=True, dtype=float)
            if dummies.shape[1] == 0:
                row[f"frac_{factor}"] = 0.0
                continue

            X_new = np.column_stack([X_prev, dummies.to_numpy(dtype=float)])
            sse_new = sse(y, X_new)
            explained = max(0.0, sse_prev - sse_new)
            row[f"frac_{factor}"] = explained / sst
            X_prev = X_new
            sse_prev = sse_new

        row["r_squared"] = max(0.0, min(1.0, 1.0 - (sse_prev / sst)))
        row["frac_residual"] = sse_prev / sst
        rows.append(row)

    return pd.DataFrame(rows)


def run_preview(args: argparse.Namespace) -> None:
    obs_filter = args.obs_filter
    if args.collection_id:
        dataset_ids = dataset_ids_for_collection(args.census_version, args.collection_id)
        obs_filter = add_dataset_filter(obs_filter, dataset_ids)

    df = preview_cells(
        census_version=args.census_version,
        obs_filter=obs_filter,
        groupby=args.groupby,
        out_csv=Path(args.out) if args.out else None,
    )
    print(df.head(args.head).to_string(index=False))
    print(f"\nTotal groups: {df.shape[0]}")


def run_score(args: argparse.Namespace) -> None:
    obs_filter = args.obs_filter
    if args.collection_id:
        dataset_ids = dataset_ids_for_collection(args.census_version, args.collection_id)
        obs_filter = add_dataset_filter(obs_filter, dataset_ids)

    df = compute_scores(
        census_version=args.census_version,
        obs_filter=obs_filter,
        modules=REDOX_MODULES,
        zscore_group=args.zscore_group,
        max_cells=args.max_cells,
        random_state=args.random_state,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".csv":
        df.to_csv(out, index=False)
    else:
        df.to_parquet(out, index=False)

    print(f"Wrote {df.shape[0]} rows x {df.shape[1]} columns to {out}")


def run_pseudobulk(args: argparse.Namespace) -> None:
    in_file = Path(args.input)
    if in_file.suffix.lower() == ".csv":
        df = pd.read_csv(in_file)
    else:
        df = pd.read_parquet(in_file)

    out_df = pseudobulk_scores(df, group_cols=args.group_cols, min_cells=args.min_cells)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".csv":
        out_df.to_csv(out, index=False)
    else:
        out_df.to_parquet(out, index=False)

    print(f"Wrote {out_df.shape[0]} pseudobulk rows x {out_df.shape[1]} columns to {out}")


def run_score_pseudobulk(args: argparse.Namespace) -> None:
    obs_filter = args.obs_filter
    if args.collection_id:
        dataset_ids = dataset_ids_for_collection(args.census_version, args.collection_id)
        obs_filter = add_dataset_filter(obs_filter, dataset_ids)

    df = compute_scores(
        census_version=args.census_version,
        obs_filter=obs_filter,
        modules=REDOX_MODULES,
        zscore_group=args.zscore_group,
        max_cells=args.max_cells,
        random_state=args.random_state,
    )
    out_df = pseudobulk_scores(df, group_cols=args.group_cols, min_cells=args.min_cells)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".csv":
        out_df.to_csv(out, index=False)
    else:
        out_df.to_parquet(out, index=False)

    print(f"Wrote {out_df.shape[0]} pseudobulk rows x {out_df.shape[1]} columns to {out}")


def run_batch_score_pseudobulk(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "manifest.csv"

    base_obs_filter = args.base_obs_filter
    if args.collection_id:
        dataset_ids = dataset_ids_for_collection(args.census_version, args.collection_id)
        base_obs_filter = add_dataset_filter(base_obs_filter, dataset_ids)

    tissues = args.tissues
    if not tissues:
        tissues = distinct_obs_values(args.census_version, base_obs_filter, "tissue_general")

    manifest_rows = []
    if args.resume and manifest_path.exists():
        manifest_rows = pd.read_csv(manifest_path).to_dict("records")
    completed_tissues = {row["tissue_general"] for row in manifest_rows}

    for tissue in tissues:
        if args.resume and tissue in completed_tissues:
            print(f"Skipping tissue_general={tissue}; already present in manifest")
            continue

        obs_filter = f"{base_obs_filter} and tissue_general == '{tissue}'"
        print(f"Scoring+pseudobulking tissue_general={tissue} ...")
        df = compute_scores(
            census_version=args.census_version,
            obs_filter=obs_filter,
            modules=REDOX_MODULES,
            zscore_group=args.zscore_group,
            max_cells=args.max_cells,
            random_state=args.random_state,
        )
        out_df = pseudobulk_scores(df, group_cols=args.group_cols, min_cells=args.min_cells)

        out_file = outdir / f"pseudobulk_tissue_general_{tissue}.parquet"
        out_df.to_parquet(out_file, index=False)
        manifest_rows.append(
            {
                "tissue_general": tissue,
                "n_scored_cells": int(df.shape[0]),
                "n_pseudobulk_rows": int(out_df.shape[0]),
                "output_file": out_file.name,
            }
        )
        pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

        if out_df.empty:
            print(f"  no qualifying pseudobulk rows; wrote empty table to {out_file}")
        else:
            print(f"  wrote {out_df.shape[0]} pseudobulk rows to {out_file}")

    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    print(f"Wrote manifest to {manifest_path}")


def run_combine(args: argparse.Namespace) -> None:
    input_paths = []
    if args.input_dir:
        input_paths.extend(sorted(Path(args.input_dir).glob(args.pattern)))
    if args.inputs:
        input_paths.extend([Path(p) for p in args.inputs])

    # Deduplicate while preserving order.
    seen = set()
    unique_paths = []
    for path in input_paths:
        key = str(path.resolve())
        if key not in seen:
            seen.add(key)
            unique_paths.append(path)

    df = combine_tables(unique_paths)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".csv":
        df.to_csv(out, index=False)
    else:
        df.to_parquet(out, index=False)

    print(f"Wrote combined table with {df.shape[0]} rows x {df.shape[1]} columns to {out}")


def run_variance_partition(args: argparse.Namespace) -> None:
    in_file = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if in_file.suffix.lower() == ".csv":
        df = pd.read_csv(in_file)
    else:
        df = pd.read_parquet(in_file)

    result = variance_partition(df, factors=args.factors)
    result.to_csv(outdir / "variance_partition.csv", index=False)

    long_df = result.melt(
        id_vars=["module", "n_rows", "r_squared"],
        value_vars=[f"frac_{factor}" for factor in args.factors if f"frac_{factor}" in result.columns] + ["frac_residual"],
        var_name="component",
        value_name="fraction",
    )
    if not long_df.empty:
        plt.figure(figsize=(12, 6))
        pivot = long_df.pivot(index="module", columns="component", values="fraction").fillna(0)
        pivot.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="tab20")
        plt.ylabel("Fraction of variance")
        plt.title("Variance Partition by Redox Module")
        plt.xticks(rotation=60, ha="right")
        plt.tight_layout()
        plt.savefig(outdir / "variance_partition_stacked.png", dpi=180)
        plt.close()

    print(f"Wrote variance partition artifacts to {outdir}")


def run_analyze(args: argparse.Namespace) -> None:
    in_file = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if in_file.suffix.lower() == ".csv":
        df = pd.read_csv(in_file)
    else:
        df = pd.read_parquet(in_file)

    score_cols = [c for c in df.columns if c.startswith("score_") and "_z_in_" not in c]
    if not score_cols:
        raise ValueError("No score_ columns found in input file.")

    # 1) Module correlation heatmap
    corr = df[score_cols].corr(method="spearman")
    corr.to_csv(outdir / "module_correlation_spearman.csv")

    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, cmap="coolwarm", center=0, square=True)
    plt.title("Redox Module Correlation (Spearman)")
    plt.tight_layout()
    plt.savefig(outdir / "module_correlation_heatmap.png", dpi=180)
    plt.close()

    # 2) Donor x module mean heatmap (for one selected cell type)
    if args.cell_type_for_donor and "cell_type" in df.columns and "donor_id" in df.columns:
        sub = df[df["cell_type"] == args.cell_type_for_donor].copy()
        if not sub.empty:
            donor_means = sub.groupby("donor_id")[score_cols].mean()
            donor_means.to_csv(outdir / "donor_module_means.csv")

            plt.figure(figsize=(12, 8))
            sns.heatmap(donor_means, cmap="vlag", center=0)
            plt.title(f"Donor Mean Module Scores in {args.cell_type_for_donor}")
            plt.tight_layout()
            plt.savefig(outdir / "donor_module_heatmap.png", dpi=180)
            plt.close()

    # 3) Top varying modules barplot by donor-level variance
    if "donor_id" in df.columns:
        donor_summary = df.groupby("donor_id")[score_cols].mean()
        donor_var = donor_summary.var(axis=0).sort_values(ascending=False)
        donor_var.to_csv(outdir / "donor_level_module_variance.csv", header=["variance"])

        top = donor_var.head(min(10, donor_var.shape[0]))
        plt.figure(figsize=(10, 5))
        sns.barplot(x=top.index, y=top.values, color="#4C78A8")
        plt.xticks(rotation=60, ha="right")
        plt.ylabel("Variance of donor means")
        plt.title("Top Modules by Donor-Level Variability")
        plt.tight_layout()
        plt.savefig(outdir / "top_donor_variability_modules.png", dpi=180)
        plt.close()

    print(f"Wrote analysis artifacts to {outdir}")


def run_batch(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base_obs_filter = args.base_obs_filter
    if args.collection_id:
        dataset_ids = dataset_ids_for_collection(args.census_version, args.collection_id)
        base_obs_filter = add_dataset_filter(base_obs_filter, dataset_ids)

    for tissue in args.tissues:
        obs_filter = f"{base_obs_filter} and tissue_general == '{tissue}'"
        print(f"Scoring tissue_general={tissue} ...")
        df = compute_scores(
            census_version=args.census_version,
            obs_filter=obs_filter,
            modules=REDOX_MODULES,
            zscore_group=args.zscore_group,
            max_cells=args.max_cells,
            random_state=args.random_state,
        )

        out_file = outdir / f"scores_tissue_general_{tissue}.parquet"
        df.to_parquet(out_file, index=False)
        if df.empty:
            print(f"  no cells found; wrote empty table to {out_file}")
        else:
            print(f"  wrote {df.shape[0]} rows to {out_file}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Memory-safe redox module scoring with CELLxGENE Census")
    p.add_argument("--census-version", default="2025-11-08", help="CELLxGENE Census release version")

    sub = p.add_subparsers(dest="command", required=True)

    p_preview = sub.add_parser("preview", help="Preview cell counts by metadata group for a filtered subset")
    p_preview.add_argument("--obs-filter", required=True, help="Census obs_value_filter")
    p_preview.add_argument("--collection-id", default=None, help="Optional collection_id to restrict matching cells")
    p_preview.add_argument("--groupby", default="tissue_general", help="obs column to group counts by")
    p_preview.add_argument("--out", help="Optional CSV path for group counts")
    p_preview.add_argument("--head", type=int, default=20, help="Rows to print")
    p_preview.set_defaults(func=run_preview)

    p_score = sub.add_parser("score", help="Compute redox scores for one filtered subset")
    p_score.add_argument("--obs-filter", required=True, help="Census obs_value_filter")
    p_score.add_argument("--collection-id", default=None, help="Optional collection_id to restrict matching cells")
    p_score.add_argument("--zscore-group", default="cell_type", help="Group for within-group z-scoring; use empty string to skip")
    p_score.add_argument("--max-cells", type=int, default=None, help="Optional random cell cap for fast experiments")
    p_score.add_argument("--random-state", type=int, default=7, help="Random seed for --max-cells sampling")
    p_score.add_argument("--out", required=True, help="Output path (.parquet recommended, .csv also supported)")
    p_score.set_defaults(func=run_score)

    p_pseudobulk = sub.add_parser("pseudobulk", help="Aggregate scored cells into donor/cell-type/tissue summaries")
    p_pseudobulk.add_argument("--input", required=True, help="Scored table input (.parquet or .csv)")
    p_pseudobulk.add_argument("--group-cols", nargs="+", default=["donor_id", "cell_type", "tissue_general"], help="Grouping columns for pseudobulk aggregation")
    p_pseudobulk.add_argument("--min-cells", type=int, default=20, help="Minimum cells required per pseudobulk group")
    p_pseudobulk.add_argument("--out", required=True, help="Output path (.parquet recommended, .csv also supported)")
    p_pseudobulk.set_defaults(func=run_pseudobulk)

    p_score_pseudobulk = sub.add_parser("score-pseudobulk", help="Compute scores for a filtered subset and immediately pseudobulk them")
    p_score_pseudobulk.add_argument("--obs-filter", required=True, help="Census obs_value_filter")
    p_score_pseudobulk.add_argument("--collection-id", default=None, help="Optional collection_id to restrict matching cells")
    p_score_pseudobulk.add_argument("--zscore-group", default="cell_type", help="Group for within-group z-scoring before aggregation; use empty string to skip")
    p_score_pseudobulk.add_argument("--max-cells", type=int, default=None, help="Optional random cell cap for fast experiments")
    p_score_pseudobulk.add_argument("--random-state", type=int, default=7, help="Random seed for --max-cells sampling")
    p_score_pseudobulk.add_argument("--group-cols", nargs="+", default=["donor_id", "cell_type", "tissue_general"], help="Grouping columns for pseudobulk aggregation")
    p_score_pseudobulk.add_argument("--min-cells", type=int, default=20, help="Minimum cells required per pseudobulk group")
    p_score_pseudobulk.add_argument("--out", required=True, help="Output path (.parquet recommended, .csv also supported)")
    p_score_pseudobulk.set_defaults(func=run_score_pseudobulk)

    p_batch_score_pseudobulk = sub.add_parser("batch-score-pseudobulk", help="Loop across tissues, score cells, and write one pseudobulk file per tissue")
    p_batch_score_pseudobulk.add_argument("--base-obs-filter", default="is_primary_data == True", help="Base filter shared across tissues")
    p_batch_score_pseudobulk.add_argument("--collection-id", default=None, help="Optional collection_id to restrict matching cells")
    p_batch_score_pseudobulk.add_argument("--tissues", nargs="*", default=None, help="Optional list of tissue_general values; if omitted, discover automatically")
    p_batch_score_pseudobulk.add_argument("--zscore-group", default="cell_type", help="Group for within-group z-scoring before aggregation; use empty string to skip")
    p_batch_score_pseudobulk.add_argument("--max-cells", type=int, default=None, help="Optional random cell cap per tissue")
    p_batch_score_pseudobulk.add_argument("--random-state", type=int, default=7, help="Random seed for per-tissue sampling")
    p_batch_score_pseudobulk.add_argument("--group-cols", nargs="+", default=["donor_id", "cell_type", "tissue_general"], help="Grouping columns for pseudobulk aggregation")
    p_batch_score_pseudobulk.add_argument("--min-cells", type=int, default=20, help="Minimum cells required per pseudobulk group")
    p_batch_score_pseudobulk.add_argument("--resume", action="store_true", help="Skip tissues already recorded in manifest.csv")
    p_batch_score_pseudobulk.add_argument("--outdir", required=True, help="Output directory for per-tissue pseudobulk tables")
    p_batch_score_pseudobulk.set_defaults(func=run_batch_score_pseudobulk)

    p_combine = sub.add_parser("combine", help="Combine per-tissue pseudobulk files into one analysis table")
    p_combine.add_argument("--input-dir", default=None, help="Directory containing tables to combine")
    p_combine.add_argument("--pattern", default="*.parquet", help="Glob pattern for files inside --input-dir")
    p_combine.add_argument("--inputs", nargs="*", default=None, help="Optional explicit list of files to combine")
    p_combine.add_argument("--out", required=True, help="Output path (.parquet recommended, .csv also supported)")
    p_combine.set_defaults(func=run_combine)

    p_varpart = sub.add_parser("variance-partition", help="Estimate variance fractions explained by donor, cell type, and tissue")
    p_varpart.add_argument("--input", required=True, help="Combined pseudobulk input (.parquet or .csv)")
    p_varpart.add_argument("--factors", nargs="+", default=["cell_type", "tissue_general", "donor_id"], help="Categorical factors for variance partitioning")
    p_varpart.add_argument("--outdir", required=True, help="Directory for variance partition outputs")
    p_varpart.set_defaults(func=run_variance_partition)

    p_analyze = sub.add_parser("analyze", help="Generate visual analysis from a scored file")
    p_analyze.add_argument("--input", required=True, help="Scored table input (.parquet or .csv)")
    p_analyze.add_argument("--outdir", required=True, help="Directory for plots and summary tables")
    p_analyze.add_argument("--cell-type-for-donor", default=None, help="Optional cell type for donor heatmap")
    p_analyze.set_defaults(func=run_analyze)

    p_batch = sub.add_parser("batch", help="Compute redox scores tissue-by-tissue")
    p_batch.add_argument("--base-obs-filter", default="is_primary_data == True", help="Base filter shared across tissues")
    p_batch.add_argument("--collection-id", default=None, help="Optional collection_id to restrict matching cells")
    p_batch.add_argument("--tissues", nargs="+", required=True, help="List of tissue_general values")
    p_batch.add_argument("--zscore-group", default="cell_type", help="Group for within-group z-scoring; use empty string to skip")
    p_batch.add_argument("--max-cells", type=int, default=None, help="Optional random cell cap per tissue")
    p_batch.add_argument("--random-state", type=int, default=7, help="Random seed for per-tissue sampling")
    p_batch.add_argument("--outdir", required=True, help="Output directory for per-tissue parquet files")
    p_batch.set_defaults(func=run_batch)

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if hasattr(args, "zscore_group") and args.zscore_group == "":
        args.zscore_group = None

    args.func(args)


if __name__ == "__main__":
    main()
