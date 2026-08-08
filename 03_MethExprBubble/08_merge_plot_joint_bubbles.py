#!/usr/bin/env python3
"""Merge summaries and draw aligned gene-body methylation/expression bubbles."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

from workflow_config import RESULT_ROOT, ordered_cell_types


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--markers", type=Path, default=RESULT_ROOT / "02_markers" / "marker_genes.tsv"
    )
    parser.add_argument(
        "--methylation",
        type=Path,
        default=RESULT_ROOT
        / "06_region_methylation"
        / "celltype_gene_region_methylation.tsv",
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=RESULT_ROOT / "07_expression" / "celltype_marker_gene_expression.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "08_joint_plot")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--expression-z-clip", type=float, default=2.0)
    parser.add_argument("--min-dot-size", type=float, default=8.0)
    parser.add_argument("--max-dot-size", type=float, default=170.0)
    return parser.parse_args()


def parse_numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")


def build_plot_table(
    markers: pd.DataFrame, methylation: pd.DataFrame, expression: pd.DataFrame
) -> pd.DataFrame:
    markers = markers.sort_values("plot_row", kind="stable").copy()
    cell_types = ordered_cell_types(expression["cell_type"].astype(str).unique().tolist())
    base = pd.MultiIndex.from_product(
        [markers["gene_symbol"].astype(str), cell_types, ["gene_body"]],
        names=["gene_symbol", "cell_type", "region"],
    ).to_frame(index=False)
    marker_meta = markers[
        [
            "gene_symbol",
            "gene_id",
            "marker_cell_type",
            "log2FC",
            "adjusted_pvalue",
            "expression_fraction",
            "rank",
            "plot_row",
        ]
    ].rename(
        columns={
            "expression_fraction": "marker_expression_fraction",
            "rank": "marker_rank",
        }
    )
    base = base.merge(marker_meta, on="gene_symbol", how="left", validate="many_to_one")
    methyl_columns = [
        "gene_symbol",
        "cell_type",
        "region",
        "mean_methylation_probability_raw",
        "mean_methylation_probability",
        "valid_cells",
        "total_cells",
        "valid_cell_fraction",
        "passes_coverage_filter",
        "mapped_sample_dmr_count",
        "effective_cpg_weight_sum",
        "mean_contributing_dmrs_per_valid_cell",
        "region_dmr_covered_bp",
        "region_length_bp",
        "region_dmr_coverage_fraction",
    ]
    base = base.merge(
        methylation[methyl_columns],
        on=["gene_symbol", "cell_type", "region"],
        how="left",
        validate="one_to_one",
    )
    expression_columns = [
        "gene_symbol",
        "cell_type",
        "rna_cells",
        "mean_normalized_expression",
        "scaled_mean_expression",
        "expression_fraction",
    ]
    base = base.merge(
        expression[expression_columns],
        on=["gene_symbol", "cell_type"],
        how="left",
        validate="many_to_one",
    )
    base = base.rename(
        columns={"mean_methylation_probability": "methylation_probability"}
    )
    cell_order = {value: index for index, value in enumerate(cell_types)}
    region_order = {"gene_body": 0}
    base["_cell_order"] = base["cell_type"].map(cell_order)
    base["_region_order"] = base["region"].map(region_order)
    return base.sort_values(
        ["plot_row", "_cell_order", "_region_order"], kind="stable"
    ).drop(columns=["_cell_order", "_region_order"])


def dot_sizes(fractions: np.ndarray, minimum: float, maximum: float) -> np.ndarray:
    fractions = np.clip(fractions, 0, 1)
    return minimum + fractions * (maximum - minimum)


def group_boundaries(markers: pd.DataFrame):
    ordered = markers.sort_values("plot_row", kind="stable").reset_index(drop=True)
    groups = []
    start = 0
    for index in range(1, len(ordered) + 1):
        if index == len(ordered) or ordered.loc[index, "marker_cell_type"] != ordered.loc[start, "marker_cell_type"]:
            groups.append((ordered.loc[start, "marker_cell_type"], start, index - 1))
            start = index
    return groups


def draw_panel(
    ax: plt.Axes,
    table: pd.DataFrame,
    value_column: str,
    fraction_column: str,
    title: str,
    cmap: str,
    norm: Normalize,
    genes: list[str],
    cell_types: list[str],
    minimum_size: float,
    maximum_size: float,
) -> None:
    gene_index = {gene: index for index, gene in enumerate(genes)}
    cell_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    values = pd.to_numeric(table[value_column], errors="coerce").to_numpy(float)
    fractions = pd.to_numeric(table[fraction_column], errors="coerce").to_numpy(float)
    valid = np.isfinite(values) & np.isfinite(fractions)
    x = table["cell_type"].map(cell_index).to_numpy(int)
    y = table["gene_symbol"].map(gene_index).to_numpy(int)
    ax.scatter(
        x[valid],
        y[valid],
        c=values[valid],
        s=dot_sizes(fractions[valid], minimum_size, maximum_size),
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        rasterized=True,
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlim(-0.6, len(cell_types) - 0.4)
    ax.set_ylim(len(genes) - 0.4, -0.6)
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=90, fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.grid(color="#e5e5e5", linewidth=0.45)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.6)


def main() -> None:
    args = parse_args()
    if args.dpi < 72 or args.expression_z_clip <= 0:
        raise ValueError("Invalid plotting parameters")
    for path in (args.markers, args.methylation, args.expression):
        if not path.is_file():
            raise FileNotFoundError(path)
    markers = pd.read_csv(args.markers, sep="\t")
    methylation = pd.read_csv(args.methylation, sep="\t", na_values=["NA"])
    expression = pd.read_csv(args.expression, sep="\t", na_values=["NA"])
    parse_numeric(
        methylation,
        [
            "mean_methylation_probability_raw",
            "mean_methylation_probability",
            "valid_cell_fraction",
        ],
    )
    parse_numeric(
        expression,
        ["mean_normalized_expression", "scaled_mean_expression", "expression_fraction"],
    )
    plot_table = build_plot_table(markers, methylation, expression)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "joint_methylation_expression_bubble_data.tsv.gz"
    plot_table.to_csv(table_path, sep="\t", index=False, na_rep="NA", compression="gzip")

    markers = markers.sort_values("plot_row", kind="stable").reset_index(drop=True)
    genes = markers["gene_symbol"].astype(str).tolist()
    cell_types = ordered_cell_types(expression["cell_type"].astype(str).unique().tolist())
    gene_body = plot_table.loc[plot_table["region"] == "gene_body"]
    expression_once = gene_body

    height = max(10.0, 0.24 * len(genes) + 3.3)
    width = max(17.0, 1.05 * len(cell_types) + 5.0)
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(width, height),
        sharey=True,
        constrained_layout=True,
    )
    methyl_norm = Normalize(0, 1)
    expression_norm = Normalize(-args.expression_z_clip, args.expression_z_clip)
    draw_panel(
        axes[0],
        gene_body,
        "methylation_probability",
        "valid_cell_fraction",
        "Gene-body methylation",
        "Reds",
        methyl_norm,
        genes,
        cell_types,
        args.min_dot_size,
        args.max_dot_size,
    )
    draw_panel(
        axes[1],
        expression_once,
        "scaled_mean_expression",
        "expression_fraction",
        "Gene expression",
        "RdBu_r",
        expression_norm,
        genes,
        cell_types,
        args.min_dot_size,
        args.max_dot_size,
    )
    axes[0].set_yticklabels(genes, fontsize=8)
    axes[0].set_ylabel("Top marker genes (grouped by marker cell type)", fontsize=11)
    axes[1].tick_params(axis="y", labelleft=False)
    for ax in axes:
        ax.set_xlabel("Cell type", fontsize=10)

    groups = group_boundaries(markers)
    palette = plt.get_cmap("tab20")
    for group_index, (cell_type, start, end) in enumerate(groups):
        for ax in axes:
            if end < len(genes) - 1:
                ax.axhline(end + 0.5, color="#444444", linewidth=0.8)
        midpoint = (start + end) / 2
        axes[0].text(
            -0.28,
            midpoint,
            cell_type,
            transform=axes[0].get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=8,
            fontweight="bold",
            color=palette(group_index % 20),
            clip_on=False,
        )

    methyl_mappable = plt.cm.ScalarMappable(norm=methyl_norm, cmap="Reds")
    expression_mappable = plt.cm.ScalarMappable(norm=expression_norm, cmap="RdBu_r")
    methyl_bar = figure.colorbar(methyl_mappable, ax=axes[0], shrink=0.55, pad=0.015)
    methyl_bar.set_label("Mean methylation probability", fontsize=9)
    expression_bar = figure.colorbar(expression_mappable, ax=axes[1], shrink=0.55, pad=0.015)
    expression_bar.set_label("Gene-wise Z-score of mean expression", fontsize=9)

    legend_handles = [
        axes[1].scatter([], [], s=dot_sizes(np.array([fraction]), args.min_dot_size, args.max_dot_size)[0], color="#777777", label=f"{fraction:.0%}")
        for fraction in (0.2, 0.5, 1.0)
    ]
    axes[1].legend(
        handles=legend_handles,
        title="Cell fraction\n(size)",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.98),
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    figure.suptitle(
        "DMR-supported gene-body methylation and marker-gene expression",
        fontsize=16,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.002,
        "Methylation represents only DMR-covered portions of gene bodies; blank = insufficient or missing methylation coverage.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#555555",
    )
    png_path = args.output_dir / "gene_body_methylation_expression_top5_marker_bubble_plot.png"
    pdf_path = args.output_dir / "gene_body_methylation_expression_top5_marker_bubble_plot.pdf"
    figure.savefig(png_path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    qc = {
        "marker_genes": len(genes),
        "cell_types": len(cell_types),
        "gene_body_plotted_values": int(gene_body["methylation_probability"].notna().sum()),
        "expression_plotted_values": int(expression_once["scaled_mean_expression"].notna().sum()),
        "methylation_color_range": [0, 1],
        "expression_zscore_clip": args.expression_z_clip,
        "png": str(png_path),
        "pdf": str(pdf_path),
        "plot_data": str(table_path),
    }
    (args.output_dir / "plot_summary.json").write_text(json.dumps(qc, indent=2) + "\n")
    print(f"[OK] {png_path}")


if __name__ == "__main__":
    main()
