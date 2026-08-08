#!/usr/bin/env python3
"""Select globally unique positive Top5 marker genes for common cell types."""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from workflow_config import (
    EXCLUDED_CELL_TYPES,
    MARKER_ADJUSTED_P_MAX,
    MARKER_EXPRESSION_FRACTION_MIN,
    MARKER_LOG2FC_MIN,
    MARKERS_PER_CELL_TYPE,
    RESULT_ROOT,
    RNA_H5AD,
    ordered_cell_types,
    strip_ensembl_version,
)


def gene_identity_table(adata: sc.AnnData) -> pd.DataFrame:
    raw_var = adata.raw.var.copy()
    raw_var.index = raw_var.index.astype(str)
    names = pd.Index(raw_var.index)
    ensembl_like = names.str.startswith("ENSG").mean() > 0.5
    symbol_column = next(
        (
            column
            for column in ("gene_symbol", "gene_symbols", "gene_name", "symbol")
            if column in raw_var.columns
        ),
        None,
    )
    id_column = next(
        (
            column
            for column in ("gene_id", "gene_ids", "ensembl_id")
            if column in raw_var.columns
        ),
        None,
    )
    symbols = (
        raw_var[symbol_column].astype(str)
        if ensembl_like and symbol_column is not None
        else pd.Series(names, index=raw_var.index, dtype=str)
    )
    gene_ids = (
        raw_var[id_column].astype(str)
        if id_column is not None
        else (
            pd.Series(names, index=raw_var.index, dtype=str)
            if ensembl_like
            else pd.Series("", index=raw_var.index, dtype=str)
        )
    )
    result = pd.DataFrame(
        {
            "feature_name": names,
            "gene_symbol": symbols.to_numpy(),
            "gene_id": [strip_ensembl_version(value) for value in gene_ids],
        }
    )
    if result["gene_symbol"].isin({"", "nan", "None"}).any():
        raise ValueError("RNA raw.var contains missing gene symbols")
    return result


def bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def main() -> int:
    audit_crosswalk = RESULT_ROOT / "01_audit" / "cell_type_crosswalk.tsv"
    if not audit_crosswalk.is_file():
        raise FileNotFoundError("Run 01_audit_joint_inputs.py first")
    crosswalk = pd.read_csv(audit_crosswalk, sep="\t")
    common_types = crosswalk.loc[
        bool_series(crosswalk["use_in_joint_plot"]), "cell_type"
    ].astype(str).tolist()
    common_types = ordered_cell_types(common_types)

    output_dir = RESULT_ROOT / "02_markers"
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(RNA_H5AD)
    if adata.raw is None:
        raise ValueError("RNA h5ad has no adata.raw")
    required = {"cell_type_integrated", "exclude_from_main_analysis"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"RNA obs lacks columns: {sorted(missing)}")
    cell_types = adata.obs["cell_type_integrated"].astype(str)
    keep = ~bool_series(adata.obs["exclude_from_main_analysis"])
    keep &= cell_types.isin(common_types)
    keep &= ~cell_types.isin(EXCLUDED_CELL_TYPES)
    adata = adata[keep].copy()
    adata.obs["cell_type_integrated"] = pd.Categorical(
        adata.obs["cell_type_integrated"].astype(str), categories=common_types
    )
    adata.obs["cell_type_integrated"] = (
        adata.obs["cell_type_integrated"].cat.remove_unused_categories()
    )
    observed_types = list(adata.obs["cell_type_integrated"].cat.categories)
    if len(observed_types) < 2:
        raise ValueError("Fewer than two common RNA cell types")

    feature_table = gene_identity_table(adata)
    feature_to_symbol = dict(
        zip(feature_table["feature_name"], feature_table["gene_symbol"])
    )
    feature_to_id = dict(zip(feature_table["feature_name"], feature_table["gene_id"]))

    key = "meth_expr_bubble_markers"
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type_integrated",
        method="wilcoxon",
        use_raw=True,
        pts=True,
        n_genes=adata.raw.n_vars,
        key_added=key,
    )

    frames: list[pd.DataFrame] = []
    for cell_type in observed_types:
        frame = sc.get.rank_genes_groups_df(adata, group=cell_type, key=key)
        frame["marker_cell_type"] = cell_type
        frames.append(frame)
    candidates = pd.concat(frames, ignore_index=True)
    required_columns = {"names", "logfoldchanges", "pvals_adj"}
    missing_columns = required_columns.difference(candidates.columns)
    if missing_columns:
        raise ValueError(f"Scanpy marker output lacks: {sorted(missing_columns)}")
    if "pct_nz_group" not in candidates.columns:
        raise ValueError("Scanpy marker output lacks pct_nz_group; pts=True is required")

    candidates["feature_name"] = candidates["names"].astype(str)
    candidates["gene_symbol"] = candidates["feature_name"].map(feature_to_symbol)
    candidates["gene_id"] = candidates["feature_name"].map(feature_to_id).fillna("")
    candidates = candidates.rename(
        columns={
            "logfoldchanges": "log2FC",
            "pvals_adj": "adjusted_pvalue",
            "pvals": "pvalue",
            "pct_nz_group": "expression_fraction",
            "pct_nz_reference": "rest_expression_fraction",
        }
    )
    numeric_columns = ["log2FC", "adjusted_pvalue", "expression_fraction"]
    for column in numeric_columns:
        candidates[column] = pd.to_numeric(candidates[column], errors="coerce")
    passing = candidates.loc[
        np.isfinite(candidates["log2FC"])
        & np.isfinite(candidates["adjusted_pvalue"])
        & np.isfinite(candidates["expression_fraction"])
        & (candidates["adjusted_pvalue"] < MARKER_ADJUSTED_P_MAX)
        & (candidates["log2FC"] > MARKER_LOG2FC_MIN)
        & (candidates["expression_fraction"] >= MARKER_EXPRESSION_FRACTION_MIN)
        & candidates["gene_symbol"].notna()
    ].copy()
    passing = passing.sort_values(
        ["gene_symbol", "log2FC", "adjusted_pvalue", "marker_cell_type"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    assigned = passing.drop_duplicates("gene_symbol", keep="first").copy()
    order_map = {value: index for index, value in enumerate(observed_types)}
    assigned["_cell_type_order"] = assigned["marker_cell_type"].map(order_map)
    assigned = assigned.sort_values(
        ["_cell_type_order", "log2FC", "adjusted_pvalue", "gene_symbol"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    selected = assigned.groupby("marker_cell_type", sort=False, group_keys=False).head(
        MARKERS_PER_CELL_TYPE
    )
    selected = selected.copy()
    selected["rank"] = selected.groupby("marker_cell_type", sort=False).cumcount() + 1
    selected["plot_row"] = np.arange(1, len(selected) + 1)

    output_columns = [
        "gene_symbol",
        "gene_id",
        "feature_name",
        "marker_cell_type",
        "log2FC",
        "adjusted_pvalue",
        "pvalue",
        "expression_fraction",
        "rest_expression_fraction",
        "rank",
        "plot_row",
    ]
    for column in output_columns:
        if column not in selected.columns:
            selected[column] = np.nan
    selected[output_columns].to_csv(
        output_dir / "marker_genes.tsv", sep="\t", index=False
    )
    passing.to_csv(output_dir / "marker_candidates_passing.tsv.gz", sep="\t", index=False)
    pd.DataFrame(
        {
            "cell_type": observed_types,
            "rna_cells": [
                int((adata.obs["cell_type_integrated"] == value).sum())
                for value in observed_types
            ],
            "selected_unique_markers": [
                int((selected["marker_cell_type"] == value).sum())
                for value in observed_types
            ],
        }
    ).to_csv(output_dir / "marker_summary.tsv", sep="\t", index=False)
    print(f"Selected {len(selected)} globally unique marker genes")
    print(f"Completed: {output_dir / 'marker_genes.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
