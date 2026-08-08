#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCANPY_PYTHON="${SCANPY_PYTHON:-/share/home/rzli/miniconda3/envs/scanpy310/bin/python}"
ACTION="${1:-all}"
CELL_JOBS="${2:-32}"

if [[ ! "$CELL_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: CELL_JOBS must be a positive integer: $CELL_JOBS" >&2
    exit 2
fi
if [[ ! -x "$SCANPY_PYTHON" ]]; then
    echo "ERROR: Python is not executable: $SCANPY_PYTHON" >&2
    exit 1
fi

run_audit() {
    echo "[1/8 RUN] audit RNA, DMR matrices, cell annotations and cov files"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/01_audit_joint_inputs.py"
}

run_markers() {
    echo "[2/8 RUN] one-vs-rest RNA marker selection and global Top5 assignment"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/02_select_top5_markers.py"
}

run_regions() {
    echo "[3/8 RUN] hg38 promoter and gene-body annotation"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/03_make_hg38_gene_regions.py"
}

run_counts() {
    echo "[4/8 RUN] unique CpG counts per cell and DMR; rolling workers=$CELL_JOBS"
    count_args=(all --cell-jobs "$CELL_JOBS")
    if [[ "${FORCE_COUNTS:-0}" == "1" ]]; then
        count_args+=(--force)
    fi
    "$SCANPY_PYTHON" "$SCRIPT_DIR/04_compute_dmr_unique_cpg_counts.py" "${count_args[@]}"
}

run_mapping() {
    echo "[5/8 RUN] DMR overlap with marker promoters and gene bodies"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/05_map_dmrs_to_gene_regions.py"
}

run_methylation() {
    echo "[6/8 RUN] cell-level and cell-type regional methylation"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/06_compute_gene_region_methylation.py"
}

run_expression() {
    echo "[7/8 RUN] marker expression summary from adata.raw.X"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/07_summarize_marker_expression.py"
}

run_plot() {
    echo "[8/8 RUN] aligned promoter/gene-body/expression bubble panels"
    "$SCANPY_PYTHON" "$SCRIPT_DIR/08_merge_plot_joint_bubbles.py"
}

case "$ACTION" in
    audit) run_audit ;;
    markers) run_markers ;;
    regions) run_regions ;;
    counts) run_counts ;;
    map) run_mapping ;;
    methylation) run_methylation ;;
    expression) run_expression ;;
    plot) run_plot ;;
    all)
        run_audit
        run_markers
        run_regions
        run_counts
        run_mapping
        run_methylation
        run_expression
        run_plot
        echo "[ALL OK] methylation-expression bubble workflow completed"
        ;;
    *)
        echo "Usage: $0 {audit|markers|regions|counts|map|methylation|expression|plot|all} [CELL_JOBS]" >&2
        exit 2
        ;;
esac
