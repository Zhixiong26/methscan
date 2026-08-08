#!/usr/bin/env python3
"""Shared configuration for the methylation-expression bubble workflow."""

from __future__ import annotations

import os
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(
    os.environ.get("METH_EXPR_PROJECT_ROOT", "/share/home/rzli/METHSCAN/03_MethExprBubble")
)
RESULT_ROOT = Path(os.environ.get("METH_EXPR_RESULT_ROOT", PROJECT_ROOT / "results"))

ALLCOOLS_ROOT = Path(
    os.environ.get("ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools")
)
RNA_H5AD = Path(
    os.environ.get(
        "RNA_H5AD",
        "/share/home/rzli/SCANPY/20260714/result/annotation/02_annotated_final.h5ad",
    )
)
GENCODE_GTF = Path(
    os.environ.get(
        "GENCODE_GTF", "/share/LCZX_Data/ref/gencode.v44.basic.annotation.gtf"
    )
)
PROMOTER_BED = Path(
    os.environ.get(
        "PROMOTER_BED", PROJECT_ROOT / "gencode.promoter_2kb.symbol.bed"
    )
)

QC_TAG = "qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob"
THRESHOLD = "300k"
TOP200_ANALYSIS = "heatmap_top200_rawp0p01_diff0p25"
RATIO_MATRIX_SUBDIR = "single_cell_DMR_mean_of_unique_CpG_ratios_top200"
DMR_ANNOTATION_SUBDIR = "sample_merged_hypo_DMRs_diff0p25_top200"
COV_SUBDIR = "cov_dedup_probability"

SAMPLE_NAMES = tuple(
    f"25110891_{group}{index:02d}_Met"
    for group in ("IR", "NR")
    for index in range(1, 6)
)

PREFERRED_CELL_TYPE_ORDER = (
    "B_cells",
    "B_cells_unresolved",
    "CD14_Monocytes",
    "CD16_Monocytes",
    "CD4_T_cells",
    "CD8_T_cells",
    "Cycling_cells",
    "HLAII_high_APCs",
    "MAIT_cells",
    "NK_cells",
    "Plasma_cells",
    "Treg_cells",
    "cDCs",
    "pDCs",
)
EXCLUDED_CELL_TYPES = frozenset({"Platelet_erythroid_contamination"})

PROMOTER_DEFINITION = "external GENCODE promoter BED; TSS +/- 2000 bp"
MARKER_ADJUSTED_P_MAX = 0.05
MARKER_LOG2FC_MIN = 0.25
MARKER_EXPRESSION_FRACTION_MIN = 0.10
MARKERS_PER_CELL_TYPE = 5
METHYLATION_MIN_VALID_CELLS = 10
METHYLATION_MIN_VALID_FRACTION = 0.20


def sample_short(sample_name: str) -> str:
    prefix = "25110891_"
    suffix = "_Met"
    if not sample_name.startswith(prefix) or not sample_name.endswith(suffix):
        raise ValueError(f"Unsupported sample name: {sample_name}")
    short = sample_name[len(prefix) : -len(suffix)]
    if len(short) != 4 or short[:2] not in {"IR", "NR"} or not short[2:].isdigit():
        raise ValueError(f"Unsupported sample name: {sample_name}")
    return short


def sample_root(sample_name: str) -> Path:
    return ALLCOOLS_ROOT / sample_name


def dmr_root(sample_name: str) -> Path:
    return sample_root(sample_name) / QC_TAG / f"methdiff_celltype_{THRESHOLD}"


def top200_root(sample_name: str) -> Path:
    return dmr_root(sample_name) / TOP200_ANALYSIS


def ratio_matrix_dir(sample_name: str) -> Path:
    return top200_root(sample_name) / RATIO_MATRIX_SUBDIR


def ratio_matrix_path(sample_name: str) -> Path:
    short = sample_short(sample_name)
    return ratio_matrix_dir(sample_name) / f"{short}__single_cell_DMR_mean_CpG_ratio.tsv.gz"


def dmr_annotation_path(sample_name: str) -> Path:
    short = sample_short(sample_name)
    return (
        top200_root(sample_name)
        / DMR_ANNOTATION_SUBDIR
        / f"{short}__merged_DMRs_annotation.tsv"
    )


def cov_dir(sample_name: str) -> Path:
    return sample_root(sample_name) / COV_SUBDIR


def cpg_count_dir(sample_name: str) -> Path:
    return RESULT_ROOT / "04_dmr_unique_cpg_counts" / sample_short(sample_name)


def cpg_count_matrix_path(sample_name: str) -> Path:
    short = sample_short(sample_name)
    return cpg_count_dir(sample_name) / f"{short}__single_cell_DMR_unique_CpG_count.tsv.gz"


def ordered_cell_types(values: set[str] | list[str] | tuple[str, ...]) -> list[str]:
    observed = {str(value) for value in values if str(value)}
    preferred = [value for value in PREFERRED_CELL_TYPE_ORDER if value in observed]
    remaining = sorted(observed.difference(preferred), key=str.casefold)
    return [*preferred, *remaining]


def strip_ensembl_version(value: str) -> str:
    text = str(value).strip()
    if text.startswith("ENSG") and "." in text:
        return text.split(".", 1)[0]
    return text
