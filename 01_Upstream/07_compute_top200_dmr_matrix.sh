#!/usr/bin/env bash

# Step 07: compute a separate single-cell x Top200-DMR matrix for every sample.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_SCRIPT_DIR="${RESULT_SCRIPT_DIR:-${SCRIPT_DIR}/lib/methdiff/python}"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
THRESHOLD="${THRESHOLD:-300k}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites1200000_covdedupprob}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"
SAMPLE_JOBS="${SAMPLE_JOBS:-1}"
CELL_JOBS="${CELL_JOBS:-64}"
CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

sample_short() {
    local sample_name="$1"
    [[ "$sample_name" =~ ^25110891_((IR|NR)[0-9]{2})_Met$ ]] || return 1
    printf '%s\n' "${BASH_REMATCH[1]}"
}

process_sample() {
    local sample_dir="$1"
    local sample_name="${sample_dir##*/}"
    local short dmr_root analysis_root merged_dmr_dir matrix_dir metadata cov_dir
    short="$(sample_short "$sample_name")" || return 1
    dmr_root="$sample_dir/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}"
    analysis_root="$dmr_root/heatmap_top200_rawp0p01_diff0p25"
    merged_dmr_dir="$analysis_root/sample_merged_hypo_DMRs_diff0p25_top200"
    matrix_dir="$analysis_root/single_cell_DMR_mean_of_unique_CpG_ratios_top200"
    metadata="$dmr_root/metadata/cell_metadata.tsv"
    cov_dir="$sample_dir/cov_dedup_probability"

    [[ -s "$merged_dmr_dir/merge_summary.tsv" ]] || {
        echo "ERROR: $short merged Top200 DMRs missing: $merged_dmr_dir" >&2
        return 1
    }
    [[ -s "$metadata" ]] || {
        echo "ERROR: $short cell metadata missing: $metadata" >&2
        return 1
    }
    [[ -d "$cov_dir" ]] || {
        echo "ERROR: $short deduplicated cov directory missing: $cov_dir" >&2
        return 1
    }

    if [[ -s "$matrix_dir/matrix_summary.tsv" && -s "$matrix_dir/parameters.tsv" ]]; then
        echo "[$short REUSE] single-cell DMR matrix"
        return 0
    fi
    [[ ! -e "$matrix_dir" ]] || {
        echo "ERROR: $short partial matrix output exists: $matrix_dir" >&2
        return 1
    }

    echo "[$short RUN] mean of unique CpG ratios; cell_workers=$CELL_JOBS"
    python "$RESULT_SCRIPT_DIR/06_compute_dmr_mean_of_cpg_ratios.py" \
        --cov-base-dir "$BASE_DIR" \
        --cov-dir "$cov_dir" \
        --metadata "$metadata" \
        --dmr-dir "$merged_dmr_dir" \
        --output-dir "$matrix_dir" \
        --jobs 1 \
        --cell-jobs "$CELL_JOBS" || return 1

    echo "[$short OK] DMR matrix: $matrix_dir"
}

[[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
# shellcheck disable=SC1090
source "$CONDA_INIT" || die "failed to initialize Conda"
conda activate "$CONDA_ENV" || die "failed to activate Conda env: $CONDA_ENV"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$EXPECTED_SAMPLES" || die "EXPECTED_SAMPLES must be positive"
is_positive_integer "$SAMPLE_JOBS" || die "SAMPLE_JOBS must be positive"
is_positive_integer "$CELL_JOBS" || die "CELL_JOBS must be positive"

SAMPLE_DIRS=()
while IFS= read -r sample_dir; do
    SAMPLE_DIRS+=("$sample_dir")
done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
[[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
    die "found ${#SAMPLE_DIRS[@]} samples; expected $EXPECTED_SAMPLES"

echo "Matrix concurrency: samples=$SAMPLE_JOBS cells_per_sample=$CELL_JOBS max_workers=$((SAMPLE_JOBS * CELL_JOBS))"
failures=0
for ((offset = 0; offset < ${#SAMPLE_DIRS[@]}; offset += SAMPLE_JOBS)); do
    pids=()
    names=()
    for ((i = offset; i < offset + SAMPLE_JOBS && i < ${#SAMPLE_DIRS[@]}; i++)); do
        process_sample "${SAMPLE_DIRS[$i]}" &
        pids+=("$!")
        names+=("${SAMPLE_DIRS[$i]##*/}")
    done
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            echo "[SAMPLE OK] ${names[$i]}"
        else
            echo "[SAMPLE FAIL] ${names[$i]}" >&2
            failures=$((failures + 1))
        fi
    done
done

[[ "$failures" -eq 0 ]] || die "$failures sample(s) failed matrix calculation"
echo "[ALL SAMPLES OK] independent Top200 DMR matrices complete"
