#!/usr/bin/env python3
"""Aggregate DMR probabilities into single-cell promoter/gene-body values."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from workflow_config import (
    METHYLATION_MIN_VALID_CELLS,
    METHYLATION_MIN_VALID_FRACTION,
    RESULT_ROOT,
    SAMPLE_NAMES,
    cpg_count_matrix_path,
    ordered_cell_types,
    ratio_matrix_dir,
    ratio_matrix_path,
    sample_short,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overlaps",
        type=Path,
        default=RESULT_ROOT
        / "05_dmr_gene_region_map"
        / "dmr_marker_gene_region_overlaps.tsv.gz",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "06_region_methylation"
    )
    parser.add_argument("--min-valid-cells", type=int, default=METHYLATION_MIN_VALID_CELLS)
    parser.add_argument(
        "--min-valid-fraction", type=float, default=METHYLATION_MIN_VALID_FRACTION
    )
    return parser.parse_args()


def read_wide_matrix(
    path: Path, value_kind: str
) -> tuple[list[str], list[tuple[str, int, int, str]], np.ndarray]:
    coords: list[tuple[str, int, int, str]] = []
    rows: list[np.ndarray] = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:4] != ["chrom", "start", "end", "dmr_id"]:
            raise ValueError(f"Unexpected matrix header: {path}")
        cells = header[4:]
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: {len(fields)} columns, expected {len(header)}"
                )
            coords.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
            if value_kind == "ratio":
                row = np.array(
                    [np.nan if x in {"NA", "", "nan"} else float(x) for x in fields[4:]],
                    dtype=np.float64,
                )
            else:
                row = np.array([int(x) for x in fields[4:]], dtype=np.float64)
            rows.append(row)
    return cells, coords, np.vstack(rows)


def read_annotations(sample: str, expected_cells: list[str]) -> dict[str, str]:
    path = ratio_matrix_dir(sample) / "cell_annotations.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    mapping: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"cell", "cell_type"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing cell/cell_type in {path}")
        for row in reader:
            mapping[row["cell"]] = row["cell_type"]
    missing = [cell for cell in expected_cells if cell not in mapping]
    if missing:
        raise ValueError(f"{path}: missing annotation for {missing[0]}")
    return mapping


def read_overlaps(path: Path) -> tuple[dict[str, list[dict[str, str]]], list[dict[str, str]]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    all_rows: list[dict[str, str]] = []
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_sample[row["sample"]].append(row)
            all_rows.append(row)
    return by_sample, all_rows


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    records = sorted(set(intervals))
    total = 0
    current_start, current_end = records[0]
    for start, end in records[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def compute_region_coverage(rows: list[dict[str, str]]):
    loci: dict[tuple[str, str], dict[tuple[str, int, int], list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in rows:
        key = (row["gene_symbol"], row["region"])
        locus = (row["chrom"], int(row["region_start"]), int(row["region_end"]))
        loci[key][locus].append((int(row["overlap_start"]), int(row["overlap_end"])))
    output = {}
    for key, locus_map in loci.items():
        covered = sum(interval_union_length(intervals) for intervals in locus_map.values())
        length = sum(end - start for _, start, end in locus_map)
        output[key] = (covered, length, covered / length if length else np.nan)
    return output


def main() -> None:
    args = parse_args()
    if args.min_valid_cells < 1:
        raise ValueError("--min-valid-cells must be >= 1")
    if not 0 <= args.min_valid_fraction <= 1:
        raise ValueError("--min-valid-fraction must be between 0 and 1")
    if not args.overlaps.is_file():
        raise FileNotFoundError(args.overlaps)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    overlaps_by_sample, overlap_rows = read_overlaps(args.overlaps)
    region_coverage = compute_region_coverage(overlap_rows)
    mapped_dmr_counts: dict[tuple[str, str], int] = defaultdict(int)
    for key in {
        (r["sample"], r["dmr_id"], r["gene_symbol"], r["region"])
        for r in overlap_rows
    }:
        mapped_dmr_counts[(key[2], key[3])] += 1

    pooled_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    pooled_weights: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    pooled_dmr_contrib: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    total_cells_by_type: dict[str, int] = defaultdict(int)
    sample_summaries: list[dict[str, object]] = []

    cell_output = args.output_dir / "single_cell_gene_region_methylation.tsv.gz"
    cell_fields = [
        "sample",
        "sample_short",
        "cell",
        "cell_type",
        "gene_symbol",
        "region",
        "methylation_probability",
        "effective_cpg_weight",
        "contributing_dmrs",
    ]
    with gzip.open(cell_output, "wt", newline="") as cell_handle:
        cell_writer = csv.DictWriter(cell_handle, fieldnames=cell_fields, delimiter="\t")
        cell_writer.writeheader()

        for sample in SAMPLE_NAMES:
            ratio_cells, ratio_coords, ratios = read_wide_matrix(
                ratio_matrix_path(sample), "ratio"
            )
            count_cells, count_coords, counts = read_wide_matrix(
                cpg_count_matrix_path(sample), "count"
            )
            if ratio_cells != count_cells or ratio_coords != count_coords:
                raise ValueError(f"Ratio/count matrix mismatch for {sample}")
            annotations = read_annotations(sample, ratio_cells)
            for cell in ratio_cells:
                total_cells_by_type[annotations[cell]] += 1

            dmr_index = {coord[3]: index for index, coord in enumerate(ratio_coords)}
            grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
            for row in overlaps_by_sample.get(sample, []):
                grouped[(row["gene_symbol"], row["region"])].append(row)

            sample_valid_values = 0
            for (gene_symbol, region), mapped_rows in grouped.items():
                numerator = np.zeros(len(ratio_cells), dtype=np.float64)
                denominator = np.zeros(len(ratio_cells), dtype=np.float64)
                contributing = np.zeros(len(ratio_cells), dtype=np.int32)
                for mapped in mapped_rows:
                    index = dmr_index.get(mapped["dmr_id"])
                    if index is None:
                        raise ValueError(
                            f"{sample}: mapped DMR absent from ratio matrix: {mapped['dmr_id']}"
                        )
                    probability = ratios[index]
                    weight = counts[index] * float(mapped["overlap_fraction"])
                    valid = np.isfinite(probability) & (weight > 0)
                    numerator[valid] += probability[valid] * weight[valid]
                    denominator[valid] += weight[valid]
                    contributing[valid] += 1
                values = np.full(len(ratio_cells), np.nan, dtype=np.float64)
                valid = denominator > 0
                values[valid] = numerator[valid] / denominator[valid]
                sample_valid_values += int(valid.sum())

                for column, cell in enumerate(ratio_cells):
                    if not valid[column]:
                        continue
                    cell_type = annotations[cell]
                    key = (gene_symbol, region, cell_type)
                    value = float(values[column])
                    weight_value = float(denominator[column])
                    pooled_values[key].append(value)
                    pooled_weights[key].append(weight_value)
                    pooled_dmr_contrib[key].append(int(contributing[column]))
                    cell_writer.writerow(
                        {
                            "sample": sample,
                            "sample_short": sample_short(sample),
                            "cell": cell,
                            "cell_type": cell_type,
                            "gene_symbol": gene_symbol,
                            "region": region,
                            "methylation_probability": f"{value:.8g}",
                            "effective_cpg_weight": f"{weight_value:.8g}",
                            "contributing_dmrs": int(contributing[column]),
                        }
                    )
            sample_summaries.append(
                {
                    "sample": sample,
                    "cells": len(ratio_cells),
                    "dmrs": len(ratio_coords),
                    "mapped_gene_regions": len(grouped),
                    "valid_cell_gene_region_values": sample_valid_values,
                }
            )
            print(
                f"[{sample_short(sample)} OK] cells={len(ratio_cells)} mapped_gene_regions={len(grouped)}",
                flush=True,
            )

    cell_types = ordered_cell_types(list(total_cells_by_type))
    gene_regions = sorted(
        {(r["gene_symbol"], r["region"]) for r in overlap_rows},
        key=lambda x: (x[0].casefold(), x[1]),
    )
    summary_fields = [
        "gene_symbol",
        "region",
        "cell_type",
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
    summary_path = args.output_dir / "celltype_gene_region_methylation.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for gene_symbol, region in gene_regions:
            covered_bp, region_length, coverage_fraction = region_coverage[(gene_symbol, region)]
            for cell_type in cell_types:
                key = (gene_symbol, region, cell_type)
                values = pooled_values.get(key, [])
                weights = pooled_weights.get(key, [])
                contrib = pooled_dmr_contrib.get(key, [])
                valid_cells = len(values)
                total_cells = total_cells_by_type[cell_type]
                fraction = valid_cells / total_cells if total_cells else 0.0
                raw_mean = float(np.mean(values)) if values else np.nan
                passes = (
                    valid_cells >= args.min_valid_cells
                    and fraction >= args.min_valid_fraction
                )
                writer.writerow(
                    {
                        "gene_symbol": gene_symbol,
                        "region": region,
                        "cell_type": cell_type,
                        "mean_methylation_probability_raw": "NA" if not np.isfinite(raw_mean) else f"{raw_mean:.8g}",
                        "mean_methylation_probability": "NA" if not passes else f"{raw_mean:.8g}",
                        "valid_cells": valid_cells,
                        "total_cells": total_cells,
                        "valid_cell_fraction": f"{fraction:.8g}",
                        "passes_coverage_filter": "yes" if passes else "no",
                        "mapped_sample_dmr_count": mapped_dmr_counts[(gene_symbol, region)],
                        "effective_cpg_weight_sum": f"{sum(weights):.8g}",
                        "mean_contributing_dmrs_per_valid_cell": "NA" if not contrib else f"{np.mean(contrib):.8g}",
                        "region_dmr_covered_bp": covered_bp,
                        "region_length_bp": region_length,
                        "region_dmr_coverage_fraction": f"{coverage_fraction:.8g}",
                    }
                )

    with (args.output_dir / "sample_methylation_summary.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(sample_summaries[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sample_summaries)
    parameters = {
        "region_probability": "sum(DMR_probability * unique_CpG_count * overlap_fraction) / sum(unique_CpG_count * overlap_fraction)",
        "celltype_mean": "unweighted arithmetic mean of valid single-cell region probabilities",
        "min_valid_cells": args.min_valid_cells,
        "min_valid_fraction": args.min_valid_fraction,
        "missing_rule": "NA; never replace with zero",
        "interpretation": "DMR-covered portion of promoter/gene body, not all regional CpGs",
    }
    (args.output_dir / "methylation_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )
    print(f"[OK] {summary_path}")


if __name__ == "__main__":
    main()
