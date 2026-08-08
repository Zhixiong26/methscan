#!/usr/bin/env bash

# Step 02: deduplicate cov files for all samples using the probability rule.
#
# Usage:
#   bash 02_deduplicate_cov_by_probability.sh all [sample_jobs] [file_jobs_per_sample]
#   bash 02_deduplicate_cov_by_probability.sh one <cov_dir> <output_dir> [file_jobs]
#
# Backward-compatible single-sample usage:
#   bash 02_deduplicate_cov_by_probability.sh <cov_dir> <output_dir> [file_jobs]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMPLEMENTATION="$SCRIPT_DIR/lib/deduplicate_cov_by_probability_one_sample.sh"
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

    echo "Deduplicating ${#sample_dirs[@]} samples: sample_jobs=$sample_jobs file_jobs_per_sample=$file_jobs"
    for sample_dir in "${sample_dirs[@]}"; do
        sample="${sample_dir##*/}"
        bash "$IMPLEMENTATION" \
            "$sample_dir/cov" \
            "$sample_dir/cov_dedup_probability" \
            "$file_jobs" &
        pids+=("$!")
        names+=("$sample")
        if [[ "${#pids[@]}" -ge "$sample_jobs" ]]; then
            wait_batch
        fi
    done
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch
    [[ "$failures" -eq 0 ]] || die "$failures sample(s) failed cov deduplication"
    echo "[ALL SAMPLES OK] cov deduplication complete"
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

