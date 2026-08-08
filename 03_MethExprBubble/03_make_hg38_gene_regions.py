#!/usr/bin/env python3
"""Create hg38 promoter/gene-body BED files and marker-only subsets."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from workflow_config import (
    GENCODE_GTF,
    PROMOTER_DOWNSTREAM_BP,
    PROMOTER_UPSTREAM_BP,
    RESULT_ROOT,
    strip_ensembl_version,
)


PRIMARY_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, default=GENCODE_GTF)
    parser.add_argument(
        "--markers",
        type=Path,
        default=RESULT_ROOT / "02_markers" / "marker_genes.tsv",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "03_gene_regions"
    )
    parser.add_argument("--promoter-upstream", type=int, default=PROMOTER_UPSTREAM_BP)
    parser.add_argument(
        "--promoter-downstream", type=int, default=PROMOTER_DOWNSTREAM_BP
    )
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition(" ")
        result[key] = value.strip().strip('"')
    return result


def read_genes(gtf: Path) -> list[dict[str, object]]:
    genes: list[dict[str, object]] = []
    with open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in PRIMARY_CHROMS:
                continue
            attrs = parse_attributes(fields[8])
            symbol = attrs.get("gene_name", "").strip()
            gene_id = strip_ensembl_version(attrs.get("gene_id", ""))
            if not symbol or not gene_id:
                continue
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            genes.append(
                {
                    "chrom": fields[0],
                    "start": start0,
                    "end": end0,
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    "strand": fields[6],
                }
            )
    return genes


def make_regions(
    genes: list[dict[str, object]], upstream: int, downstream: int
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    promoters: list[dict[str, object]] = []
    bodies: list[dict[str, object]] = []
    for gene in genes:
        body = dict(gene)
        body["region"] = "gene_body"
        bodies.append(body)

        strand = str(gene["strand"])
        if strand == "+":
            tss = int(gene["start"])
            p_start = max(0, tss - upstream)
            p_end = tss + downstream
        elif strand == "-":
            tss = int(gene["end"]) - 1
            p_start = max(0, tss - downstream)
            p_end = tss + upstream
        else:
            continue
        promoter = dict(gene)
        promoter.update(start=p_start, end=p_end, region="promoter", tss=tss)
        promoters.append(promoter)
    return promoters, bodies


def write_bed(path: Path, rows: list[dict[str, object]]) -> None:
    header = [
        "chrom",
        "start",
        "end",
        "gene_id",
        "gene_symbol",
        "gene_type",
        "strand",
        "region",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.promoter_upstream < 0 or args.promoter_downstream < 0:
        raise ValueError("Promoter distances must be non-negative")
    if not args.gtf.is_file():
        raise FileNotFoundError(args.gtf)
    if not args.markers.is_file():
        raise FileNotFoundError(args.markers)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    genes = read_genes(args.gtf)
    if not genes:
        raise RuntimeError(f"No primary-chromosome gene records read from {args.gtf}")
    promoters, bodies = make_regions(
        genes, args.promoter_upstream, args.promoter_downstream
    )

    with args.markers.open() as handle:
        marker_rows = list(csv.DictReader(handle, delimiter="\t"))
    marker_symbols = {row["gene_symbol"] for row in marker_rows}

    symbol_loci: dict[str, set[tuple[str, int, int, str]]] = defaultdict(set)
    for gene in genes:
        symbol_loci[str(gene["gene_symbol"])].add(
            (str(gene["chrom"]), int(gene["start"]), int(gene["end"]), str(gene["strand"]))
        )

    marker_promoters = [r for r in promoters if r["gene_symbol"] in marker_symbols]
    marker_bodies = [r for r in bodies if r["gene_symbol"] in marker_symbols]
    sort_key = lambda row: (str(row["chrom"]), int(row["start"]), int(row["end"]), str(row["gene_symbol"]))
    promoters.sort(key=sort_key)
    bodies.sort(key=sort_key)
    marker_promoters.sort(key=sort_key)
    marker_bodies.sort(key=sort_key)

    write_bed(args.output_dir / "gencode_v44_primary_promoters.bed", promoters)
    write_bed(args.output_dir / "gencode_v44_primary_gene_bodies.bed", bodies)
    write_bed(args.output_dir / "marker_gene_promoters.bed", marker_promoters)
    write_bed(args.output_dir / "marker_gene_bodies.bed", marker_bodies)

    audit_path = args.output_dir / "marker_gene_region_audit.tsv"
    with audit_path.open("w", newline="") as handle:
        fields = ["gene_symbol", "marker_cell_type", "matched_loci", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in marker_rows:
            symbol = row["gene_symbol"]
            n_loci = len(symbol_loci.get(symbol, set()))
            writer.writerow(
                {
                    "gene_symbol": symbol,
                    "marker_cell_type": row["marker_cell_type"],
                    "matched_loci": n_loci,
                    "status": "unmatched" if n_loci == 0 else "unique" if n_loci == 1 else "ambiguous",
                }
            )

    status_counts = Counter(
        "unmatched" if not symbol_loci.get(s) else "unique" if len(symbol_loci[s]) == 1 else "ambiguous"
        for s in marker_symbols
    )
    summary = {
        "assembly": "GRCh38/hg38 (from configured GENCODE v44 GTF)",
        "gtf": str(args.gtf),
        "primary_gene_records": len(genes),
        "marker_genes": len(marker_symbols),
        "marker_status": dict(status_counts),
        "promoter_upstream_bp": args.promoter_upstream,
        "promoter_downstream_bp": args.promoter_downstream,
        "promoter_coordinate_rule": "strand-aware BED half-open interval",
    }
    (args.output_dir / "region_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
