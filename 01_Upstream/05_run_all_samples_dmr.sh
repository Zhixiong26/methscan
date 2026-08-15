#!/usr/bin/env bash

# Step 05: run cell-type pairwise DMRs independently within each sample.
# Each sample uses its own filtered_data_single_<threshold>/smoothed background.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
DMR_SCRIPT="${DMR_SCRIPT:-${SCRIPT_DIR}/lib/methdiff/run_single_sample_dmr.sh}"
ANNOTATION_CSV="${ANNOTATION_CSV:-/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv}"
THRESHOLD="${THRESHOLD:-300k}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"
DEFAULT_PREPARE_JOBS="${DEFAULT_PREPARE_JOBS:-2}"
DEFAULT_SAMPLE_JOBS="${DEFAULT_SAMPLE_JOBS:-2}"
DEFAULT_COMPARISON_JOBS="${DEFAULT_COMPARISON_JOBS:-2}"
DEFAULT_THREADS="${DEFAULT_THREADS:-24}"

usage() {
    cat <<'EOF'
Usage:
  bash 05_run_all_samples_dmr.sh prepare [sample_jobs]
  bash 05_run_all_samples_dmr.sh run [sample_jobs] [comparison_jobs_per_sample] [threads_per_comparison]
  bash 05_run_all_samples_dmr.sh status
  bash 05_run_all_samples_dmr.sh summarize

Examples:
  bash 05_run_all_samples_dmr.sh prepare 2
  bash 05_run_all_samples_dmr.sh run 2 2 24
  bash 05_run_all_samples_dmr.sh status

The run example permits at most 2 x 2 x 24 = 96 MethSCAn diff threads.
All comparisons remain within one sample; no merged MethSCAn input is used.
The default annotation is the Scanpy 20260810 Result0810 annotation. Override it
with ANNOTATION_CSV=/path/to/02_cell_annotation_all_cells.csv when needed.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

collect_samples() {
    SAMPLE_DIRS=()
    while IFS= read -r sample_dir; do
        SAMPLE_DIRS+=("$sample_dir")
    done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
    [[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
        die "found ${#SAMPLE_DIRS[@]} sample directories; expected $EXPECTED_SAMPLES"
}

sample_short() {
    local sample_name="$1"
    if [[ "$sample_name" =~ ^25110891_((IR|NR)[0-9]{2})_Met$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        return 1
    fi
}

run_one_sample() {
    local sample_dir="$1"
    local action="$2"
    shift 2
    local sample_name="${sample_dir##*/}"
    local short
    short="$(sample_short "$sample_name")" || {
        echo "ERROR: invalid sample directory: $sample_dir" >&2
        return 1
    }

    echo ">>> $short $action"
    env SAMPLE_NAME="$sample_name" SAMPLE_SHORT="$short" THRESHOLD="$THRESHOLD" \
        ANNOTATION_CSV="$ANNOTATION_CSV" \
        bash "$DMR_SCRIPT" "$action" "$@"
}

run_parallel() {
    local action="$1"
    local max_jobs="$2"
    shift 2
    local sample_dir i failures=0
    local -a pids=()
    local -a names=()

    wait_batch() {
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[SAMPLE OK] ${names[$i]} $action"
            else
                echo "[SAMPLE FAIL] ${names[$i]} $action" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    }

    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        run_one_sample "$sample_dir" "$action" "$@" &
        pids+=("$!")
        names+=("${sample_dir##*/}")
        if [[ "${#pids[@]}" -ge "$max_jobs" ]]; then
            wait_batch
        fi
    done
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch

    if [[ "$failures" -gt 0 ]]; then
        echo "[ALL SAMPLES FAIL] $failures sample(s) failed during $action" >&2
        return 1
    fi
    echo "[ALL SAMPLES OK] $action complete"
}

[[ -s "$DMR_SCRIPT" ]] || die "single-sample DMR implementation missing: $DMR_SCRIPT"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$EXPECTED_SAMPLES" || die "EXPECTED_SAMPLES must be positive"
action="${1:-}"
if [[ "$action" == -h || "$action" == --help || "$action" == help ]]; then
    usage
    exit 0
fi
collect_samples

case "$action" in
    prepare)
        sample_jobs="${2:-$DEFAULT_PREPARE_JOBS}"
        is_positive_integer "$sample_jobs" || die "sample_jobs must be positive"
        run_parallel prepare "$sample_jobs"
        ;;
    run)
        sample_jobs="${2:-$DEFAULT_SAMPLE_JOBS}"
        comparison_jobs="${3:-$DEFAULT_COMPARISON_JOBS}"
        threads="${4:-$DEFAULT_THREADS}"
        is_positive_integer "$sample_jobs" || die "sample_jobs must be positive"
        is_positive_integer "$comparison_jobs" || die "comparison_jobs must be positive"
        is_positive_integer "$threads" || die "threads must be positive"
        echo "Independent DMR concurrency: samples=$sample_jobs comparisons_per_sample=$comparison_jobs threads_per_comparison=$threads"
        echo "Maximum requested MethSCAn threads: $((sample_jobs * comparison_jobs * threads))"
        run_parallel run "$sample_jobs" "$comparison_jobs" "$threads"
        ;;
    status)
        for sample_dir in "${SAMPLE_DIRS[@]}"; do
            run_one_sample "$sample_dir" status || exit 1
        done
        ;;
    summarize)
        run_parallel summarize 1
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
