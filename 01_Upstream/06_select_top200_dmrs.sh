#!/usr/bin/env bash

# Step 06: independently select Top200 hypo-DMRs for every sample.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_SCRIPT_DIR="${RESULT_SCRIPT_DIR:-${SCRIPT_DIR}/lib/methdiff/python}"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
THRESHOLD="${THRESHOLD:-300k}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites1200000_covdedupprob}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"
SAMPLE_JOBS="${SAMPLE_JOBS:-10}"
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
    local short dmr_root analysis_root source_dmr_dir fallback_dmr_dir hypo_dir merged_dmr_dir
    local expected_fallback existing_fallback input_mode
    local -a input_args=()
    short="$(sample_short "$sample_name")" || return 1
    dmr_root="$sample_dir/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}"
    analysis_root="$dmr_root/heatmap_top200_rawp0p01_diff0p25"
    source_dmr_dir="$dmr_root/results"
    fallback_dmr_dir="$dmr_root/rawp_fallback_no_null_dmrs/results"
    hypo_dir="$analysis_root/celltype_hypo_DMRs_diff0p25_top200"
    merged_dmr_dir="$analysis_root/sample_merged_hypo_DMRs_diff0p25_top200"
    input_args=(--result-dir "$source_dmr_dir")
    expected_fallback=""
    if [[ -d "$fallback_dmr_dir" ]]; then
        expected_fallback="$fallback_dmr_dir"
        input_args+=(--fallback-result-dir "$fallback_dmr_dir")
        input_mode="standard+rawp-fallback"
    else
        input_mode="standard-only"
    fi

    [[ -d "$source_dmr_dir" ]] || {
        echo "ERROR: $short standard DMR results directory missing: $source_dmr_dir" >&2
        return 1
    }
    if [[ ! -s "$dmr_root/summary_celltype_pairwise.tsv" &&
        ! -s "$dmr_root/rawp_fallback_no_null_dmrs/fallback_status.tsv" ]]; then
        echo "ERROR: $short has neither a completed standard DMR summary nor a raw-p fallback status" >&2
        return 1
    fi

    if [[ -s "$hypo_dir/overall_summary.tsv" && -s "$hypo_dir/parameters.tsv" ]]; then
        existing_fallback="$(awk -F $'\t' '$1 == "rawp_fallback_result_dirs" {print $2; exit}' "$hypo_dir/parameters.tsv")"
        if [[ "$existing_fallback" != "$expected_fallback" ]]; then
            echo "ERROR: $short existing Top200 selection has different DMR inputs" >&2
            echo "       existing fallback: ${existing_fallback:-none}" >&2
            echo "       required fallback: ${expected_fallback:-none}" >&2
            echo "       Archive both $hypo_dir and $merged_dmr_dir, then rerun." >&2
            return 1
        fi
        echo "[$short 1/2 REUSE] Top200 hypo-DMR selection ($input_mode)"
    else
        [[ ! -e "$hypo_dir" ]] || {
            echo "ERROR: $short partial Top200 selection exists: $hypo_dir" >&2
            return 1
        }
        echo "[$short 1/2 RUN] Top200 hypo-DMR selection"
        python "$RESULT_SCRIPT_DIR/05_extract_celltype_hypo_dmrs_top1500.py" \
            "${input_args[@]}" \
            --output-dir "$hypo_dir" \
            --raw-p 0.01 \
            --min-abs-diff 0.25 \
            --top-dmrs-per-cell 200 \
            --sample "$short" \
            --jobs 1 || return 1
    fi

    if [[ -s "$merged_dmr_dir/merge_summary.tsv" ]]; then
        echo "[$short 2/2 REUSE] merged Top200 DMR intervals"
    else
        [[ ! -e "$merged_dmr_dir" ]] || {
            echo "ERROR: $short partial merged DMR output exists: $merged_dmr_dir" >&2
            return 1
        }
        echo "[$short 2/2 RUN] merge duplicate/overlapping DMR intervals"
        python "$RESULT_SCRIPT_DIR/02_merge_sample_dmrs.py" \
            --input-dir "$hypo_dir" \
            --output-dir "$merged_dmr_dir" \
            --jobs 1 || return 1
    fi

    echo "[$short OK] Top200 DMRs: $merged_dmr_dir"
}

[[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
# shellcheck disable=SC1090
source "$CONDA_INIT" || die "failed to initialize Conda"
conda activate "$CONDA_ENV" || die "failed to activate Conda env: $CONDA_ENV"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$EXPECTED_SAMPLES" || die "EXPECTED_SAMPLES must be positive"
is_positive_integer "$SAMPLE_JOBS" || die "SAMPLE_JOBS must be positive"

SAMPLE_DIRS=()
while IFS= read -r sample_dir; do
    SAMPLE_DIRS+=("$sample_dir")
done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
[[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
    die "found ${#SAMPLE_DIRS[@]} samples; expected $EXPECTED_SAMPLES"

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

[[ "$failures" -eq 0 ]] || die "$failures sample(s) failed Top200 processing"
echo "[ALL SAMPLES OK] independent Top200 DMR selection complete"
