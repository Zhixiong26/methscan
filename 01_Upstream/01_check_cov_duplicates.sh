#!/usr/bin/env bash

# Step 01: audit duplicate CpG coordinates in every sample's original cov files.
#
# Usage:
#   bash 01_check_cov_duplicates.sh all [sample_jobs] [file_jobs_per_sample]
#   bash 01_check_cov_duplicates.sh one <cov_dir> <output_dir> [file_jobs]
#
# Backward-compatible single-sample usage:
#   bash 01_check_cov_duplicates.sh <cov_dir> <output_dir> [file_jobs]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMPLEMENTATION="$SCRIPT_DIR/lib/check_cov_duplicates_one_sample.sh"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

run_all_samples() {
    local sample_jobs="${1:-2}"
    local file_jobs="${2:-48}"
    local sample_dir sample failures=0 i
    local -a sample_dirs=() pids=() names=()

    is_positive_integer "$sample_jobs" || die "sample_jobs must be positive"
    is_positive_integer "$file_jobs" || die "file_jobs_per_sample must be positive"
    mapfile -t sample_dirs < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
    [[ "${#sample_dirs[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
        die "found ${#sample_dirs[@]} sample directories; expected $EXPECTED_SAMPLES"

    wait_batch() {
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[SAMPLE OK] ${names[$i]}"
            else
                echo "[SAMPLE FAIL] ${names[$i]}" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    }

    echo "Auditing ${#sample_dirs[@]} samples: sample_jobs=$sample_jobs file_jobs_per_sample=$file_jobs"
    for sample_dir in "${sample_dirs[@]}"; do
        sample="${sample_dir##*/}"
        bash "$IMPLEMENTATION" \
            "$sample_dir/cov" \
            "$sample_dir/cov_duplicate_qc" \
            "$file_jobs" &
        pids+=("$!")
        names+=("$sample")
        if [[ "${#pids[@]}" -ge "$sample_jobs" ]]; then
            wait_batch
        fi
    done
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch
    [[ "$failures" -eq 0 ]] || die "$failures sample(s) failed duplicate audit"
    echo "[ALL SAMPLES OK] duplicate audit complete"
}

[[ -s "$IMPLEMENTATION" ]] || die "implementation missing: $IMPLEMENTATION"
action="${1:-all}"
case "$action" in
    all)
        run_all_samples "${2:-2}" "${3:-48}"
        ;;
    one)
        shift
        exec bash "$IMPLEMENTATION" "$@"
        ;;
    -h|--help|help)
        sed -n '3,10p' "$0"
        ;;
    *)
        exec bash "$IMPLEMENTATION" "$@"
        ;;
esac

