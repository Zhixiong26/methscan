#!/usr/bin/env python3
"""Count unique covered CpG loci per single cell and Top200 DMR."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from workflow_config import (
    RESULT_ROOT,
    SAMPLE_NAMES,
    cov_dir,
    cpg_count_dir,
    cpg_count_matrix_path,
    ratio_matrix_path,
    sample_short,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "samples",
        nargs="*",
        default=["all"],
        help="Sample names or short names (IR01); default: all",
    )
    parser.add_argument("--cell-jobs", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_samples(values: list[str]) -> list[str]:
    if not values or values == ["all"]:
        return list(SAMPLE_NAMES)
    by_short = {sample_short(s): s for s in SAMPLE_NAMES}
    selected: list[str] = []
    for value in values:
        sample = by_short.get(value, value)
        if sample not in SAMPLE_NAMES:
            raise ValueError(f"Unknown sample: {value}")
        if sample not in selected:
            selected.append(sample)
    return selected


def read_ratio_header_and_dmrs(
    path: Path,
) -> tuple[list[str], list[tuple[str, int, int, str]]]:
    dmrs: list[tuple[str, int, int, str]] = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) < 5 or header[:4] != ["chrom", "start", "end", "dmr_id"]:
            raise ValueError(f"Unexpected ratio matrix header: {path}")
        cells = header[4:]
        if len(cells) != len(set(cells)):
            raise ValueError(f"Duplicate cell names in {path}")
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t", 4)
            if len(fields) < 4:
                raise ValueError(f"Malformed row {line_number}: {path}")
            chrom, start, end, dmr_id = fields[:4]
            start_i, end_i = int(start), int(end)
            if start_i >= end_i:
                raise ValueError(f"Invalid DMR interval at row {line_number}: {path}")
            dmrs.append((chrom, start_i, end_i, dmr_id))
    if not dmrs:
        raise ValueError(f"No DMRs in {path}")
    if len({row[3] for row in dmrs}) != len(dmrs):
        raise ValueError(f"Duplicate DMR IDs in {path}")
    return cells, dmrs


def build_interval_index(dmrs: list[tuple[str, int, int, str]]):
    grouped: dict[str, list[tuple[int, int, int]]] = {}
    for index, (chrom, start, end, _) in enumerate(dmrs):
        grouped.setdefault(chrom, []).append((start, end, index))
    compact: dict[str, tuple[list[int], list[int], list[int]]] = {}
    for chrom, records in grouped.items():
        records.sort()
        previous_end = -1
        for start, end, _ in records:
            if start < previous_end:
                raise ValueError(
                    f"Overlapping DMRs found on {chrom}; unique assignment would be ambiguous"
                )
            previous_end = end
        compact[chrom] = (
            [r[0] for r in records],
            [r[1] for r in records],
            [r[2] for r in records],
        )
    return compact


def count_one_cov(
    task: tuple[str, str, int, dict[str, tuple[list[int], list[int], list[int]]]]
) -> tuple[int, np.ndarray, int, int]:
    cell, cov_path_text, n_dmrs, interval_index = task
    counts = np.zeros(n_dmrs, dtype=np.uint32)
    rows_read = 0
    assigned_unique_loci = 0
    previous_locus: tuple[str, int] | None = None
    try:
        with gzip.open(cov_path_text, "rt") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t", 3)
                if len(fields) < 3:
                    continue
                rows_read += 1
                chrom = fields[0]
                try:
                    position = int(fields[1])
                except ValueError:
                    continue
                locus = (chrom, position)
                if locus == previous_locus:
                    continue
                previous_locus = locus
                chrom_index = interval_index.get(chrom)
                if chrom_index is None:
                    continue
                starts, ends, indices = chrom_index
                local = bisect.bisect_right(starts, position) - 1
                if local >= 0 and position < ends[local]:
                    counts[indices[local]] += 1
                    assigned_unique_loci += 1
    except Exception as exc:
        raise RuntimeError(f"Failed reading {cell}: {cov_path_text}: {exc}") from exc
    return os.getpid(), counts, rows_read, assigned_unique_loci


def atomic_write_matrix(
    output: Path,
    cells: list[str],
    dmrs: list[tuple[str, int, int, str]],
    matrix: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with gzip.open(temp_path, "wt") as handle:
            handle.write("\t".join(["chrom", "start", "end", "dmr_id", *cells]) + "\n")
            for row_index, (chrom, start, end, dmr_id) in enumerate(dmrs):
                values = "\t".join(str(int(x)) for x in matrix[row_index])
                handle.write(f"{chrom}\t{start}\t{end}\t{dmr_id}\t{values}\n")
        temp_path.replace(output)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def process_sample(sample: str, cell_jobs: int, force: bool) -> dict[str, object]:
    ratio_path = ratio_matrix_path(sample)
    output = cpg_count_matrix_path(sample)
    summary_path = cpg_count_dir(sample) / "count_summary.json"
    if output.is_file() and summary_path.is_file() and not force:
        with summary_path.open() as handle:
            summary = json.load(handle)
        summary["status"] = "reused"
        print(f"[{sample_short(sample)} REUSE] {output}", flush=True)
        return summary
    if not ratio_path.is_file():
        raise FileNotFoundError(ratio_path)

    cells, dmrs = read_ratio_header_and_dmrs(ratio_path)
    interval_index = build_interval_index(dmrs)
    cov_root = cov_dir(sample)
    missing_cov = [cell for cell in cells if not (cov_root / f"{cell}.cov.gz").is_file()]
    if missing_cov:
        raise FileNotFoundError(
            f"{sample}: {len(missing_cov)} matrix cells lack cov files; first={missing_cov[0]}"
        )

    matrix = np.zeros((len(dmrs), len(cells)), dtype=np.uint32)
    total_rows = 0
    total_assigned = 0
    worker_pids: set[int] = set()
    tasks = [
        (cell, str(cov_root / f"{cell}.cov.gz"), len(dmrs), interval_index)
        for cell in cells
    ]
    print(
        f"[{sample_short(sample)} RUN] cells={len(cells)} DMRs={len(dmrs)} workers={cell_jobs}",
        flush=True,
    )
    with ProcessPoolExecutor(max_workers=cell_jobs) as executor:
        future_to_index = {
            executor.submit(count_one_cov, task): index for index, task in enumerate(tasks)
        }
        completed = 0
        for future in as_completed(future_to_index):
            column = future_to_index[future]
            pid, counts, rows_read, assigned = future.result()
            matrix[:, column] = counts
            total_rows += rows_read
            total_assigned += assigned
            worker_pids.add(pid)
            completed += 1
            if completed == len(cells) or completed % max(1, len(cells) // 20) == 0:
                print(
                    f"[{sample_short(sample)} PROGRESS] {completed}/{len(cells)} cells",
                    flush=True,
                )

    atomic_write_matrix(output, cells, dmrs, matrix)
    summary = {
        "status": "complete",
        "sample": sample,
        "ratio_matrix": str(ratio_path),
        "cov_dir": str(cov_root),
        "output": str(output),
        "cells": len(cells),
        "dmrs": len(dmrs),
        "cov_rows_read": total_rows,
        "unique_cpg_assignments": total_assigned,
        "covered_matrix_values": int(np.count_nonzero(matrix)),
        "maximum_unique_cpg_count": int(matrix.max(initial=0)),
        "worker_processes_observed": len(worker_pids),
        "interval_rule": "BED half-open: DMR_start <= cov_start < DMR_end",
        "duplicate_locus_rule": "count each chrom/start once per cell",
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[{sample_short(sample)} OK] {output}", flush=True)
    return summary


def main() -> None:
    args = parse_args()
    if args.cell_jobs < 1:
        raise ValueError("--cell-jobs must be >= 1")
    selected = resolve_samples(args.samples)
    summaries = [process_sample(s, args.cell_jobs, args.force) for s in selected]
    combined_dir = RESULT_ROOT / "04_dmr_unique_cpg_counts"
    combined_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "status",
        "cells",
        "dmrs",
        "cov_rows_read",
        "unique_cpg_assignments",
        "covered_matrix_values",
        "maximum_unique_cpg_count",
    ]
    with (combined_dir / "all_samples_count_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == "__main__":
    main()
