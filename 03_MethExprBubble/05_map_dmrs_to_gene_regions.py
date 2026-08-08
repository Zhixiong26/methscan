#!/usr/bin/env python3
"""Map each sample's Top200 DMR set to marker promoters and gene bodies."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from workflow_config import RESULT_ROOT, SAMPLE_NAMES, ratio_matrix_path, sample_short


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--region-dir", type=Path, default=RESULT_ROOT / "03_gene_regions"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "05_dmr_gene_region_map"
    )
    return parser.parse_args()


def read_regions(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    parsed = []
    for row in rows:
        parsed.append(
            {
                **row,
                "start": int(row["start"]),
                "end": int(row["end"]),
            }
        )
    return parsed


def read_dmrs(sample: str) -> list[dict[str, object]]:
    path = ratio_matrix_path(sample)
    dmrs: list[dict[str, object]] = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:4] != ["chrom", "start", "end", "dmr_id"]:
            raise ValueError(f"Unexpected header: {path}")
        for order, line in enumerate(handle, start=1):
            chrom, start, end, dmr_id = line.rstrip("\n").split("\t", 4)[:4]
            dmrs.append(
                {
                    "sample": sample,
                    "sample_short": sample_short(sample),
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "dmr_id": dmr_id,
                    "dmr_order": order,
                }
            )
    return dmrs


def overlaps_for_sample(
    dmrs: list[dict[str, object]], regions: list[dict[str, object]]
) -> tuple[list[dict[str, object]], int]:
    regions_by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for region in regions:
        regions_by_chrom[str(region["chrom"])].append(region)
    for records in regions_by_chrom.values():
        records.sort(key=lambda x: (int(x["start"]), int(x["end"])))

    candidates: list[dict[str, object]] = []
    for dmr in dmrs:
        d_start, d_end = int(dmr["start"]), int(dmr["end"])
        d_length = d_end - d_start
        for region in regions_by_chrom.get(str(dmr["chrom"]), []):
            r_start, r_end = int(region["start"]), int(region["end"])
            if r_start >= d_end:
                break
            if r_end <= d_start:
                continue
            overlap_start = max(d_start, r_start)
            overlap_end = min(d_end, r_end)
            overlap_bp = overlap_end - overlap_start
            if overlap_bp <= 0:
                continue
            region_length = r_end - r_start
            candidates.append(
                {
                    **dmr,
                    "gene_id": region["gene_id"],
                    "gene_symbol": region["gene_symbol"],
                    "gene_type": region["gene_type"],
                    "strand": region["strand"],
                    "region": region["region"],
                    "region_start": r_start,
                    "region_end": r_end,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_bp": overlap_bp,
                    "dmr_length": d_length,
                    "overlap_fraction": overlap_bp / d_length,
                    "region_length": region_length,
                    "region_overlap_fraction": overlap_bp / region_length,
                }
            )

    best: dict[tuple[str, str, str, str], dict[str, object]] = {}
    duplicate_candidates = 0
    for row in candidates:
        key = (
            str(row["sample"]),
            str(row["dmr_id"]),
            str(row["gene_symbol"]),
            str(row["region"]),
        )
        previous = best.get(key)
        if previous is None or int(row["overlap_bp"]) > int(previous["overlap_bp"]):
            if previous is not None:
                duplicate_candidates += 1
            best[key] = row
        else:
            duplicate_candidates += 1
    output = sorted(
        best.values(),
        key=lambda r: (
            str(r["sample"]),
            int(r["dmr_order"]),
            str(r["region"]),
            str(r["gene_symbol"]),
        ),
    )
    return output, duplicate_candidates


def main() -> None:
    args = parse_args()
    promoter_path = args.region_dir / "marker_gene_promoters.bed"
    body_path = args.region_dir / "marker_gene_bodies.bed"
    for path in (promoter_path, body_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    regions = [*read_regions(promoter_path), *read_regions(body_path)]
    if not regions:
        raise RuntimeError("No marker gene regions available")

    all_rows: list[dict[str, object]] = []
    duplicate_candidates = 0
    sample_dmr_counts: dict[str, int] = {}
    for sample in SAMPLE_NAMES:
        dmrs = read_dmrs(sample)
        sample_dmr_counts[sample] = len(dmrs)
        rows, duplicates = overlaps_for_sample(dmrs, regions)
        all_rows.extend(rows)
        duplicate_candidates += duplicates
        print(
            f"[{sample_short(sample)}] DMRs={len(dmrs)} region_overlaps={len(rows)}",
            flush=True,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "sample_short",
        "chrom",
        "start",
        "end",
        "dmr_id",
        "dmr_order",
        "gene_id",
        "gene_symbol",
        "gene_type",
        "strand",
        "region",
        "region_start",
        "region_end",
        "overlap_start",
        "overlap_end",
        "overlap_bp",
        "dmr_length",
        "overlap_fraction",
        "region_length",
        "region_overlap_fraction",
    ]
    output_path = args.output_dir / "dmr_marker_gene_region_overlaps.tsv.gz"
    with gzip.open(output_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    counts = Counter((str(r["region"]) for r in all_rows))
    summary = {
        "samples": len(SAMPLE_NAMES),
        "dmrs_per_sample": sample_dmr_counts,
        "overlap_rows": len(all_rows),
        "overlap_rows_by_region": dict(counts),
        "marker_genes_with_any_overlap": len({r["gene_symbol"] for r in all_rows}),
        "duplicate_symbol_locus_candidates_discarded": duplicate_candidates,
        "interval_rule": "BED half-open overlap",
        "multi_mapping_rule": "DMR may map to multiple genes and both region types",
        "same_symbol_same_region_rule": "retain the locus with maximum overlap_bp",
    }
    (args.output_dir / "mapping_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(f"[OK] {output_path}")


if __name__ == "__main__":
    main()
