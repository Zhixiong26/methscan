#!/usr/bin/env python3
"""Summarize raw-layer log-normalized marker expression by cell type."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from workflow_config import EXCLUDED_CELL_TYPES, RESULT_ROOT, RNA_H5AD, ordered_cell_types


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=RNA_H5AD)
    parser.add_argument(
        "--markers", type=Path, default=RESULT_ROOT / "02_markers" / "marker_genes.tsv"
    )
    parser.add_argument(
        "--crosswalk",
        type=Path,
        default=RESULT_ROOT / "01_audit" / "cell_type_crosswalk.tsv",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "07_expression"
    )
    return parser.parse_args()


def bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def main() -> None:
    args = parse_args()
    for path in (args.h5ad, args.markers, args.crosswalk):
        if not path.is_file():
            raise FileNotFoundError(path)
    markers = pd.read_csv(args.markers, sep="\t")
    crosswalk = pd.read_csv(args.crosswalk, sep="\t")
    use_values = bool_series(crosswalk["use_in_joint_plot"])
    common_types = ordered_cell_types(
        crosswalk.loc[use_values, "cell_type"].astype(str).tolist()
    )

    adata = sc.read_h5ad(args.h5ad)
    if adata.raw is None:
        raise ValueError("RNA h5ad has no adata.raw")
    required_obs = {"cell_type_integrated", "exclude_from_main_analysis"}
    missing = required_obs.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"RNA obs lacks: {sorted(missing)}")
    cell_type_values = adata.obs["cell_type_integrated"].astype(str)
    keep = ~bool_series(adata.obs["exclude_from_main_analysis"])
    keep &= cell_type_values.isin(common_types)
    keep &= ~cell_type_values.isin(EXCLUDED_CELL_TYPES)
    adata = adata[keep].copy()
    cell_types = adata.obs["cell_type_integrated"].astype(str).to_numpy()

    feature_names = markers["feature_name"].astype(str).tolist()
    missing_features = sorted(set(feature_names).difference(adata.raw.var_names.astype(str)))
    if missing_features:
        raise ValueError(f"Marker feature absent from adata.raw: {missing_features[0]}")
    matrix = adata.raw[:, feature_names].X
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    else:
        matrix = np.asarray(matrix)

    rows: list[dict[str, object]] = []
    for cell_type in common_types:
        mask = cell_types == cell_type
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        subset = matrix[mask]
        if sparse.issparse(subset):
            means = np.asarray(subset.mean(axis=0)).ravel()
            fractions = np.asarray((subset > 0).mean(axis=0)).ravel()
        else:
            means = np.asarray(subset.mean(axis=0)).ravel()
            fractions = np.asarray((subset > 0).mean(axis=0)).ravel()
        for index, marker in markers.reset_index(drop=True).iterrows():
            rows.append(
                {
                    "gene_symbol": marker["gene_symbol"],
                    "feature_name": marker["feature_name"],
                    "marker_cell_type": marker["marker_cell_type"],
                    "marker_rank": int(marker["rank"]),
                    "plot_row": int(marker["plot_row"]),
                    "cell_type": cell_type,
                    "rna_cells": n_cells,
                    "mean_normalized_expression": float(means[index]),
                    "expression_fraction": float(fractions[index]),
                }
            )
    output = pd.DataFrame(rows)

    output["scaled_mean_expression"] = np.nan
    for _, indices in output.groupby("gene_symbol", sort=False).groups.items():
        values = output.loc[indices, "mean_normalized_expression"].to_numpy(float)
        standard_deviation = float(np.std(values, ddof=0))
        if np.isfinite(standard_deviation) and standard_deviation > 0:
            output.loc[indices, "scaled_mean_expression"] = (
                values - float(np.mean(values))
            ) / standard_deviation

    output = output.sort_values(["plot_row", "cell_type"], kind="stable")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "celltype_marker_gene_expression.tsv"
    output.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    parameters = {
        "h5ad": str(args.h5ad),
        "expression_source": "adata.raw.X",
        "mean_expression": "arithmetic mean including zero values",
        "expression_fraction": "fraction of cells with raw-layer normalized expression > 0",
        "scaled_expression": "per-gene Z-score across plotted cell types (population SD)",
        "rna_cells": int(adata.n_obs),
        "marker_genes": int(markers.shape[0]),
        "cell_types": common_types,
    }
    (args.output_dir / "expression_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )
    print(f"[OK] {output_path}")


if __name__ == "__main__":
    main()
