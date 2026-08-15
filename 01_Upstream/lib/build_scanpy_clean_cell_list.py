#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


TRUE_VALUES = {"true", "1", "yes", "y"}
FALSE_VALUES = {"false", "0", "no", "n", ""}


def short_sample(value: str) -> str:
    text = str(value).strip()
    match = re.search(r"(?:^|_)(IR|NR)(\d{2})(?:_|$)", text)
    if match is None:
        match = re.match(r"^(IR|NR)(\d{2})", text)
    if match is None:
        raise ValueError(f"Cannot derive sample from {value!r}")
    return f"{match.group(1)}{match.group(2)}"


def normalize_cell(value: str, sample_name: str, sample_short: str) -> str:
    text = str(value).strip().rsplit("/", 1)[-1]
    for suffix in (".cov.gz", ".cov", ".allc.gz"):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    prefixes = (
        f"{sample_name}__",
        f"{sample_name}_",
        f"{sample_short}__",
        f"{sample_short}_",
    )
    while True:
        for prefix in prefixes:
            if text.startswith(prefix):
                text = text[len(prefix) :]
                break
        else:
            return text


def annotation_dialect(path: Path) -> csv.Dialect:
    with path.open("r", newline="") as handle:
        sample = handle.read(8192)
    return csv.Sniffer().sniff(sample, delimiters=",\t")


def read_clean_annotation_cells(
    path: Path,
    sample_name: str,
    sample_short: str,
) -> set[str]:
    clean_cells: set[str] = set()
    dialect = annotation_dialect(path)
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        fields = set(reader.fieldnames or [])
        cell_column = next(
            (column for column in ("cell_id", "cell") if column in fields),
            None,
        )
        if cell_column is None:
            raise ValueError("Annotation lacks cell_id/cell column")
        sample_column = next(
            (column for column in ("sample", "sample_id") if column in fields),
            None,
        )
        exclude_column = next(
            (
                column
                for column in ("exclude_from_main_analysis", "exclude")
                if column in fields
            ),
            None,
        )

        for row in reader:
            cell_value = (row.get(cell_column) or "").strip()
            if not cell_value:
                raise ValueError("Annotation contains an empty cell identifier")
            row_sample = (
                short_sample(row.get(sample_column, ""))
                if sample_column is not None
                else short_sample(cell_value)
            )
            if row_sample != sample_short:
                continue

            if exclude_column is not None:
                exclude_value = (row.get(exclude_column) or "").strip().lower()
                if exclude_value not in TRUE_VALUES | FALSE_VALUES:
                    raise ValueError(
                        f"Unexpected exclusion value for {cell_value}: {exclude_value!r}"
                    )
                if exclude_value in TRUE_VALUES:
                    continue

            normalized = normalize_cell(cell_value, sample_name, sample_short)
            if normalized in clean_cells:
                raise ValueError(
                    f"Duplicate normalized annotation cell for {sample_short}: {normalized}"
                )
            clean_cells.add(normalized)

    if not clean_cells:
        raise ValueError(f"No clean annotation cells found for {sample_short}")
    return clean_cells


def read_methscan_cells(
    path: Path,
    sample_name: str,
    sample_short: str,
) -> tuple[list[str], dict[str, str]]:
    original_cells = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not original_cells:
        raise ValueError(f"MethSCAn header is empty: {path}")

    normalized_to_original: dict[str, str] = {}
    for original in original_cells:
        normalized = normalize_cell(original, sample_name, sample_short)
        if normalized in normalized_to_original:
            raise ValueError(
                f"Duplicate normalized MethSCAn cell for {sample_short}: {normalized}"
            )
        normalized_to_original[normalized] = original
    return original_cells, normalized_to_original


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an exact MethSCAn --cell-names keep list by intersecting a "
            "coverage-filtered column_header with Scanpy clean-cell annotations."
        )
    )
    parser.add_argument("--methscan-header", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--sample-short", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expected_short = short_sample(args.sample_name)
    if args.sample_short != expected_short:
        raise ValueError(
            f"sample-short={args.sample_short} does not match sample-name={args.sample_name}"
        )
    if not args.methscan_header.is_file():
        raise FileNotFoundError(args.methscan_header)
    if not args.annotation.is_file():
        raise FileNotFoundError(args.annotation)

    clean_cells = read_clean_annotation_cells(
        args.annotation,
        args.sample_name,
        args.sample_short,
    )
    original_cells, normalized_to_original = read_methscan_cells(
        args.methscan_header,
        args.sample_name,
        args.sample_short,
    )
    kept_cells = [
        original
        for original in original_cells
        if normalize_cell(original, args.sample_name, args.sample_short) in clean_cells
    ]
    if not kept_cells:
        raise ValueError(
            f"No coverage-filtered MethSCAn cells match Scanpy clean cells for {args.sample_short}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f"{args.output.name}.tmp.{os.getpid()}")
    temporary.write_text("".join(f"{cell}\n" for cell in kept_cells))
    temporary.replace(args.output)

    clean_not_coverage = len(clean_cells - set(normalized_to_original))
    print(
        "\t".join(
            (
                f"sample={args.sample_short}",
                f"coverage_cells={len(original_cells)}",
                f"scanpy_clean_cells={len(clean_cells)}",
                f"kept_cells={len(kept_cells)}",
                f"removed_after_coverage={len(original_cells) - len(kept_cells)}",
                f"clean_cells_below_or_above_coverage={clean_not_coverage}",
            )
        )
    )


if __name__ == "__main__":
    main()
