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
    PROMOTER_BED,
    PROMOTER_DEFINITION,
    RESULT_ROOT,
    strip_ensembl_version,
)


PRIMARY_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, default=GENCODE_GTF)
    parser.add_argument("--promoter-bed", type=Path, default=PROMOTER_BED)
    parser.add_argument(
        "--markers",
        type=Path,
        default=RESULT_ROOT / "02_markers" / "marker_genes.tsv",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "03_gene_regions"
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


def make_gene_bodies(genes: list[dict[str, object]]) -> list[dict[str, object]]:
    bodies: list[dict[str, object]] = []
    for gene in genes:
        body = dict(gene)
        body["region"] = "gene_body"
        bodies.append(body)
    return bodies


def read_external_promoters(
    path: Path, genes: list[dict[str, object]]
) -> tuple[list[dict[str, object]], dict[str, int]]:
    gene_type_by_id = {str(gene["gene_id"]): str(gene["gene_type"]) for gene in genes}
    promoters: list[dict[str, object]] = []
    seen: set[tuple[str, int, int, str, str, str]] = set()
    audit = {
        "input_rows": 0,
        "primary_rows": 0,
        "nonprimary_rows_excluded": 0,
        "duplicate_rows_excluded": 0,
    }
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            audit["input_rows"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected 6 BED columns, found {len(fields)}"
                )
            chrom, start_text, end_text, gene_id_text, gene_symbol, strand = fields
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: non-integer BED coordinate") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            if strand not in {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid strand: {strand}")
            if not gene_symbol:
                raise ValueError(f"{path}:{line_number}: empty gene symbol")
            if chrom not in PRIMARY_CHROMS:
                audit["nonprimary_rows_excluded"] += 1
                continue
            gene_id = strip_ensembl_version(gene_id_text)
            key = (chrom, start, end, gene_id, gene_symbol, strand)
            if key in seen:
                audit["duplicate_rows_excluded"] += 1
                continue
            seen.add(key)
            promoters.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "gene_type": gene_type_by_id.get(gene_id, ""),
                    "strand": strand,
                    "region": "promoter",
                }
            )
            audit["primary_rows"] += 1
    return promoters, audit


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
    if not args.gtf.is_file():
        raise FileNotFoundError(args.gtf)
    if not args.promoter_bed.is_file():
        raise FileNotFoundError(args.promoter_bed)
    if not args.markers.is_file():
        raise FileNotFoundError(args.markers)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    genes = read_genes(args.gtf)
    if not genes:
        raise RuntimeError(f"No primary-chromosome gene records read from {args.gtf}")
    bodies = make_gene_bodies(genes)
    promoters, promoter_audit = read_external_promoters(args.promoter_bed, genes)
    if not promoters:
        raise RuntimeError(f"No primary-chromosome promoters read from {args.promoter_bed}")

    with args.markers.open() as handle:
        marker_rows = list(csv.DictReader(handle, delimiter="\t"))
    marker_symbols = {row["gene_symbol"] for row in marker_rows}

    body_symbol_loci: dict[str, set[tuple[str, int, int, str]]] = defaultdict(set)
    for gene in genes:
        body_symbol_loci[str(gene["gene_symbol"])].add(
            (str(gene["chrom"]), int(gene["start"]), int(gene["end"]), str(gene["strand"]))
        )
    promoter_symbol_loci: dict[str, set[tuple[str, int, int, str]]] = defaultdict(set)
    for promoter in promoters:
        promoter_symbol_loci[str(promoter["gene_symbol"])].add(
            (
                str(promoter["chrom"]),
                int(promoter["start"]),
                int(promoter["end"]),
                str(promoter["strand"]),
            )
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
        fields = [
            "gene_symbol",
            "marker_cell_type",
            "promoter_loci",
            "gene_body_loci",
            "status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in marker_rows:
            symbol = row["gene_symbol"]
            promoter_loci = len(promoter_symbol_loci.get(symbol, set()))
            body_loci = len(body_symbol_loci.get(symbol, set()))
            if promoter_loci and body_loci:
                status = "complete"
            elif promoter_loci:
                status = "promoter_only"
            elif body_loci:
                status = "gene_body_only"
            else:
                status = "unmatched"
            if promoter_loci > 1 or body_loci > 1:
                status += "_ambiguous"
            writer.writerow(
                {
                    "gene_symbol": symbol,
                    "marker_cell_type": row["marker_cell_type"],
                    "promoter_loci": promoter_loci,
                    "gene_body_loci": body_loci,
                    "status": status,
                }
            )

    def marker_status(symbol: str) -> str:
        promoter_loci = len(promoter_symbol_loci.get(symbol, set()))
        body_loci = len(body_symbol_loci.get(symbol, set()))
        if promoter_loci and body_loci:
            status = "complete"
        elif promoter_loci:
            status = "promoter_only"
        elif body_loci:
            status = "gene_body_only"
        else:
            status = "unmatched"
        return status + "_ambiguous" if promoter_loci > 1 or body_loci > 1 else status

    status_counts = Counter(marker_status(symbol) for symbol in marker_symbols)
    summary = {
        "assembly": "GRCh38/hg38",
        "gtf": str(args.gtf),
        "promoter_bed": str(args.promoter_bed),
        "primary_gene_records": len(genes),
        "primary_promoter_records": len(promoters),
        "promoter_input_audit": promoter_audit,
        "marker_genes": len(marker_symbols),
        "marker_status": dict(status_counts),
        "promoter_definition": PROMOTER_DEFINITION,
        "promoter_coordinate_rule": "use external BED coordinates as BED half-open intervals",
    }
    (args.output_dir / "region_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
