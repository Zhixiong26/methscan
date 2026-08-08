#!/usr/bin/env python3
"""Audit RNA, DMR-matrix, cov, annotation, and hg38 inputs before integration."""

from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from workflow_config import (
    EXCLUDED_CELL_TYPES,
    GENCODE_GTF,
    RESULT_ROOT,
    RNA_H5AD,
    SAMPLE_NAMES,
    cov_dir,
    dmr_annotation_path,
    ordered_cell_types,
    ratio_matrix_dir,
    ratio_matrix_path,
    sample_short,
)


PRIMARY_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def fail(message: str) -> None:
    raise ValueError(message)


def bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def read_ratio_header(path: Path) -> tuple[list[str], int, set[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        if header[:4] != ["chrom", "start", "end", "dmr_id"]:
            fail(f"Unexpected ratio-matrix header: {path}")
        cells = header[4:]
        if not cells or len(cells) != len(set(cells)):
            fail(f"Missing or duplicated matrix cells: {path}")
        dmrs = 0
        chroms: set[str] = set()
        previous: tuple[str, int] | None = None
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n\r").split("\t", 4)
            if len(fields) != 5:
                fail(f"Malformed matrix row: {path}:{line_number}")
            chrom, start_text, end_text, _dmr_id, _values = fields
            start = int(start_text)
            end = int(end_text)
            if start < 0 or end <= start:
                fail(f"Invalid DMR interval: {path}:{line_number}")
            key = (chrom, start)
            if previous is not None and key < previous:
                fail(f"DMR rows are not ordered: {path}:{line_number}")
            previous = key
            chroms.add(chrom)
            dmrs += 1
    return cells, dmrs, chroms


def main() -> int:
    output_dir = RESULT_ROOT / "01_audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    required = [RNA_H5AD, GENCODE_GTF]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing global inputs: {missing}")

    adata = ad.read_h5ad(RNA_H5AD, backed="r")
    if adata.raw is None:
        fail("RNA h5ad has no adata.raw; full log-normalized genes are required")
    required_obs = {"sample", "cell_type_integrated", "exclude_from_main_analysis"}
    missing_obs = required_obs.difference(adata.obs.columns)
    if missing_obs:
        fail(f"RNA obs lacks columns: {sorted(missing_obs)}")
    rna_obs = adata.obs[list(required_obs)].copy()
    raw_var_names = pd.Index(adata.raw.var_names.astype(str))
    if raw_var_names.has_duplicates:
        fail("adata.raw.var_names contains duplicate gene identifiers")
    rna_keep = ~bool_series(rna_obs["exclude_from_main_analysis"])
    rna_types = set(
        rna_obs.loc[rna_keep, "cell_type_integrated"].astype(str)
    ).difference(EXCLUDED_CELL_TYPES)

    sample_rows: list[dict[str, object]] = []
    meth_types: set[str] = set()
    for sample_name in SAMPLE_NAMES:
        short = sample_short(sample_name)
        matrix_path = ratio_matrix_path(sample_name)
        annotation_path = ratio_matrix_dir(sample_name) / "cell_annotations.tsv"
        summary_path = ratio_matrix_dir(sample_name) / "matrix_summary.tsv"
        dmr_path = dmr_annotation_path(sample_name)
        sample_cov_dir = cov_dir(sample_name)
        required_sample = [matrix_path, annotation_path, summary_path, dmr_path]
        missing_sample = [str(path) for path in required_sample if not path.is_file()]
        if missing_sample or not sample_cov_dir.is_dir():
            fail(
                f"{sample_name}: missing inputs: {missing_sample}; cov_dir={sample_cov_dir}"
            )

        matrix_cells, dmr_count, chroms = read_ratio_header(matrix_path)
        nonprimary = sorted(chroms.difference(PRIMARY_CHROMS))
        if nonprimary:
            fail(f"{sample_name}: non-primary DMR chromosomes: {nonprimary[:10]}")

        cell_annotations = pd.read_csv(annotation_path, sep="\t", dtype=str)
        needed = {"sample", "cell", "cell_type"}
        missing_columns = needed.difference(cell_annotations.columns)
        if missing_columns:
            fail(f"{annotation_path} lacks columns: {sorted(missing_columns)}")
        observed_cells = cell_annotations["cell"].tolist()
        if matrix_cells != observed_cells:
            fail(f"{sample_name}: ratio header and cell_annotations order differ")
        sample_types = set(cell_annotations["cell_type"].dropna().astype(str))
        meth_types.update(sample_types)
        missing_cov = sum(
            not (sample_cov_dir / f"{cell}.cov.gz").is_file()
            for cell in matrix_cells
        )
        if missing_cov:
            fail(f"{sample_name}: {missing_cov} selected cells lack deduplicated cov")
        cov_files = sum(1 for _ in sample_cov_dir.glob("*.cov.gz"))
        sample_rows.append(
            {
                "sample": short,
                "matrix_cells": len(matrix_cells),
                "dmrs": dmr_count,
                "matrix_chromosomes": len(chroms),
                "selected_cell_types": len(sample_types),
                "cov_files": cov_files,
                "missing_selected_cov": missing_cov,
            }
        )

    common_types = rna_types.intersection(meth_types)
    type_order = ordered_cell_types(common_types)
    if not type_order:
        fail("RNA and methylation have no common cell types")

    pd.DataFrame(sample_rows).to_csv(
        output_dir / "sample_input_summary.tsv", sep="\t", index=False
    )
    all_types = sorted(rna_types.union(meth_types), key=str.casefold)
    pd.DataFrame(
        [
            {
                "cell_type": cell_type,
                "in_rna": cell_type in rna_types,
                "in_methylation": cell_type in meth_types,
                "use_in_joint_plot": cell_type in common_types,
                "plot_order": (
                    type_order.index(cell_type) + 1 if cell_type in common_types else ""
                ),
            }
            for cell_type in all_types
        ]
    ).to_csv(output_dir / "cell_type_crosswalk.tsv", sep="\t", index=False)

    summary = {
        "rna_h5ad": str(RNA_H5AD),
        "rna_cells": int(adata.n_obs),
        "rna_hvg_matrix_genes": int(adata.n_vars),
        "rna_raw_genes": int(adata.raw.n_vars),
        "rna_kept_cells": int(rna_keep.sum()),
        "rna_cell_types": len(rna_types),
        "methylation_cell_types": len(meth_types),
        "common_cell_types": type_order,
        "gencode_gtf": str(GENCODE_GTF),
        "gencode_release": "v44",
        "reference_assembly": "GRCh38/hg38",
        "dmr_chromosome_rule": "chr1-chr22,chrX,chrY",
        "status": "PASS",
    }
    (output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Completed: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
