#!/usr/bin/env bash

# ==============================================================================
# Methscan 全样本上游流程：状态检查、阈值分析与断点续跑
#
# 运行逻辑序号：
#   [1/8] Preflight：环境、参数、QC 配置与既有完成状态检查
#   [2/8] Prepare：cov -> compact（阈值无关，可复用）
#   [3/8] Profile：TSS QC（阈值无关，按 TSS SHA-256 复用）
#   [4/8] Filter：按 min/max sites 和 min-meth 过滤细胞
#   [5/8] Smooth：对 filtered 数据进行伪 bulk 平滑
#   [6/8] Scan：发现 VMR
#   [7/8] Matrix：生成细胞 × VMR 矩阵
#   [8/8] Summary：汇总每个样本成功/失败状态
# 说明：prepare/profile 与过滤阈值无关，优先复用已有合格产物。
#
# 用法：
#   bash 03_run_upstream_pipeline.sh status [10k|20k|30k|50k|300k]
#   bash 03_run_upstream_pipeline.sh run <threshold> [max_jobs] [threads] [sample|all]
#   bash 03_run_upstream_pipeline.sh run-to-compact <threshold> [max_jobs] [threads] [sample|all]
#   bash 03_run_upstream_pipeline.sh run-to-smooth <threshold> [max_jobs] [threads] [sample|all]
# ==============================================================================

set -uo pipefail

# ==============================================================================
# 1. 全局配置
#    所有配置均可在命令前通过同名环境变量覆盖。
# ==============================================================================

BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
DATA_TAG="${DATA_TAG:-covdedupprob}"
COV_SUBDIR="${COV_SUBDIR:-cov_dedup_probability}"
COMPACT_SUBDIR="${COMPACT_SUBDIR:-compact_data_dedup_probability}"
PROFILE_BASENAME="${PROFILE_BASENAME:-TSS_profile_dedup_probability}"
TSS_BED="${TSS_BED:-/share/LCZX_Data/ref/human_hg38_TSS.bed}"
CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"
DEFAULT_MAX_JOBS="${DEFAULT_MAX_JOBS:-1}"
DEFAULT_THREADS="${DEFAULT_THREADS:-20}"
FILTER_MIN_METH="${FILTER_MIN_METH:-55}"
FILTER_MAX_METH="${FILTER_MAX_METH:-}"
FILTER_MAX_SITES="${FILTER_MAX_SITES:-10000000}"
VALID_THRESHOLDS=(10k 20k 30k 50k 300k)
TSS_SHA256=""
SCRIPT_SHA256=""
STOP_AFTER_SMOOTH="${STOP_AFTER_SMOOTH:-0}"
STOP_AFTER_PREPARE="${STOP_AFTER_PREPARE:-0}"

FILTER_MAX_METH_LABEL="${FILTER_MAX_METH:-none}"
QC_TAG="minmeth${FILTER_MIN_METH}_maxmeth${FILTER_MAX_METH_LABEL}_maxsites${FILTER_MAX_SITES}"
QC_TAG="${QC_TAG//./p}"

# ==============================================================================
# 2. 命令行帮助与基础参数校验
# ==============================================================================

# 打印命令格式和常用示例。
usage() {
    cat <<'EOF'
Usage:
  bash 03_run_upstream_pipeline.sh status [10k|20k|30k|50k|300k]
  bash 03_run_upstream_pipeline.sh run <threshold> [max_jobs] [threads] [sample|all]
  bash 03_run_upstream_pipeline.sh run-to-compact <threshold> [max_jobs] [threads] [sample|all]
  bash 03_run_upstream_pipeline.sh run-to-smooth <threshold> [max_jobs] [threads] [sample|all]

Examples:
  bash 03_run_upstream_pipeline.sh status
  bash 03_run_upstream_pipeline.sh status 30k
  bash 03_run_upstream_pipeline.sh run-to-compact 300k 10 1 all
  bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
EOF
}

# 打印错误并终止主进程。
die() {
    echo "ERROR: $*" >&2
    exit 1
}

# 判断阈值是否属于允许的五档阈值。
is_threshold() {
    local value="$1"
    local item
    for item in "${VALID_THRESHOLDS[@]}"; do
        [[ "$value" == "$item" ]] && return 0
    done
    return 1
}

# 并发数和线程数必须是正整数。
is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

# 甲基化阈值使用 0–100 的百分数；空字符串仅允许表示“不设上限”。
is_percentage() {
    local value="$1"
    [[ "$value" =~ ^[0-9]+([.][0-9]+)?$ ]] &&
        awk -v value="$value" 'BEGIN { exit(value >= 0 && value <= 100 ? 0 : 1) }'
}

# 校验过滤配置，防止异常环境变量生成错误路径或传给 MethSCAn。
validate_filter_config() {
    is_percentage "$FILTER_MIN_METH" ||
        die "FILTER_MIN_METH must be a percentage between 0 and 100"
    if [[ -n "$FILTER_MAX_METH" ]]; then
        is_percentage "$FILTER_MAX_METH" ||
            die "FILTER_MAX_METH must be empty or a percentage between 0 and 100"
        awk -v min="$FILTER_MIN_METH" -v max="$FILTER_MAX_METH" \
            'BEGIN { exit(max >= min ? 0 : 1) }' ||
            die "FILTER_MAX_METH must be greater than or equal to FILTER_MIN_METH"
    fi
    is_positive_integer "$FILTER_MAX_SITES" ||
        die "FILTER_MAX_SITES must be a positive integer"
    [[ "$FILTER_MAX_SITES" -ge 300000 ]] ||
        die "FILTER_MAX_SITES must be at least the largest min-sites threshold (300000)"
    [[ "$QC_TAG" =~ ^[A-Za-z0-9._-]+$ ]] || die "invalid derived QC tag: $QC_TAG"
}

# ==============================================================================
# 3. 产物计数与完整性判定
#    不以“目录或日志存在”作为成功标准，检查关键非空文件。
# ==============================================================================

# 统计目录第一层中匹配给定模式的普通文件。
count_files() {
    local dir="$1"
    local pattern="${2:-*}"
    find "$dir" -maxdepth 1 -type f -name "$pattern" 2>/dev/null | wc -l
}

# 统计数据目录中的真实细胞数；每行对应一个细胞。
count_cells() {
    local dir="$1"
    local header="$dir/column_header.txt"

    if [[ ! -s "$header" ]]; then
        printf '0\n'
        return 0
    fi
    awk 'NF { n += 1 } END { print n + 0 }' "$header"
}

# compact 必须包含细胞名、细胞统计信息和至少一个染色体 NPZ。
valid_compact() {
    local dir="$1"
    [[ -s "$dir/column_header.txt" ]] &&
        [[ -s "$dir/cell_stats.csv" ]] &&
        find "$dir" -maxdepth 1 -type f -name '*.npz' -print -quit 2>/dev/null |
        grep -q .
}

# 检查 filtered 数据是否由当前过滤参数生成。
valid_filter_provenance() {
    local dir="$1"
    local min_sites="$2"
    local metadata="$dir/filter_provenance.tsv"

    [[ -s "$metadata" ]] &&
        awk -F '\t' \
            -v min_sites="$min_sites" \
            -v max_sites="$FILTER_MAX_SITES" \
            -v min_meth="$FILTER_MIN_METH" \
            -v max_meth="$FILTER_MAX_METH_LABEL" '
            $1 == "min_sites" && $2 == min_sites { a = 1 }
            $1 == "max_sites" && $2 == max_sites { b = 1 }
            $1 == "min_meth" && $2 == min_meth { c = 1 }
            $1 == "max_meth" && $2 == max_meth { d = 1 }
            END { exit(a && b && c && d ? 0 : 1) }
        ' "$metadata"
}

# filtered 目录结构与 compact 相同；NPZ 是染色体矩阵，不是单个细胞。
valid_filtered() {
    local dir="$1"
    local min_sites="$2"
    valid_compact "$dir" &&
        [[ "$(count_cells "$dir")" -gt 0 ]] &&
        valid_filter_provenance "$dir" "$min_sites"
}

# smooth 会在 filtered 目录下创建非空 smoothed 子目录。
valid_smooth() {
    local dir="$1"
    [[ -d "$dir/smoothed" ]] &&
        [[ -n "$(find "$dir/smoothed" -mindepth 1 -print -quit 2>/dev/null)" ]]
}

# scan 必须产生非空 VMRs.bed。
valid_scan() {
    [[ -s "$1/VMRs.bed" ]]
}

# matrix 必须包含非空 total_sites.csv.gz，且目录至少有 4 个文件。
valid_matrix() {
    local dir="$1"
    [[ -s "$dir/total_sites.csv.gz" ]] &&
        [[ "$(count_files "$dir")" -ge 4 ]]
}

# ==============================================================================
# 4. 样本发现与阈值无关产物复用
# ==============================================================================

# 收集并排序所有 *_Met 样本；样本数异常时立即停止。
collect_samples() {
    mapfile -t SAMPLE_DIRS < <(
        find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort
    )
    [[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
        die "found ${#SAMPLE_DIRS[@]} sample directories; expected ${EXPECTED_SAMPLES}"
}

# 所有样本的概率去重数据使用独立 QC 根目录。
qc_root_for_sample() {
    local sample_dir="$1"
    printf '%s/qc_%s_%s\n' "$sample_dir" "$QC_TAG" "$DATA_TAG"
}

# 只允许复用由概率去重 cov 生成的 compact，避免混入原始 cov 产物。
choose_compact() {
    local sample_dir="$1"
    local candidate="$sample_dir/$COMPACT_SUBDIR"
    if valid_compact "$candidate"; then
        printf '%s\n' "$candidate"
        return 0
    fi
    return 1
}

# 判断 common profile 是否由当前 TSS_BED 生成。
valid_profile() {
    local sample_dir="$1"
    local compact_dir="$2"
    local profile_file="$sample_dir/${PROFILE_BASENAME}.csv"
    local metadata_file="$sample_dir/${PROFILE_BASENAME}.meta.tsv"

    [[ -s "$profile_file" ]] &&
        [[ -s "$metadata_file" ]] &&
        awk -F '\t' -v expected="$TSS_SHA256" -v compact="$compact_dir" '
            $1 == "tss_sha256" && $2 == expected { tss_matched = 1 }
            $1 == "input_compact" && $2 == compact { compact_matched = 1 }
            END { exit(tss_matched && compact_matched ? 0 : 1) }
        ' "$metadata_file"
}

# 仅复用带有匹配 SHA-256 元数据的 common profile。
# 不再自动采用来源不明的 TSS_profile_single_*.csv。
choose_profile() {
    local sample_dir="$1"
    local compact_dir="$2"
    if valid_profile "$sample_dir" "$compact_dir"; then
        printf '%s/%s.csv\n' "$sample_dir" "$PROFILE_BASENAME"
        return 0
    fi
    return 1
}

# 判断 QC 根目录是否与本次过滤配置一致。
valid_qc_config() {
    local qc_root="$1"
    local metadata="$qc_root/pipeline_config.tsv"

    [[ -s "$metadata" ]] &&
        awk -F '\t' \
            -v qc_tag="$QC_TAG" \
            -v max_sites="$FILTER_MAX_SITES" \
            -v min_meth="$FILTER_MIN_METH" \
            -v max_meth="$FILTER_MAX_METH_LABEL" '
            $1 == "qc_tag" && $2 == qc_tag { a = 1 }
            $1 == "max_sites" && $2 == max_sites { b = 1 }
            $1 == "min_meth" && $2 == min_meth { c = 1 }
            $1 == "max_meth" && $2 == max_meth { d = 1 }
            END { exit(a && b && c && d ? 0 : 1) }
        ' "$metadata"
}

# 首次运行时创建 QC 配置；已有内容但配置不匹配时拒绝继续。
ensure_qc_config() {
    local qc_root="$1"
    local metadata="$qc_root/pipeline_config.tsv"
    local tmp="${metadata}.tmp.$$"

    if valid_qc_config "$qc_root"; then
        return 0
    fi
    if [[ -d "$qc_root" ]] &&
        [[ -n "$(find "$qc_root" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "ERROR: QC directory exists without matching configuration: $qc_root" >&2
        echo "       Archive it or use the parameters recorded in pipeline_config.tsv." >&2
        return 1
    fi

    mkdir -p "$qc_root"
    if ! {
        printf 'qc_tag\t%s\n' "$QC_TAG"
        printf 'max_sites\t%s\n' "$FILTER_MAX_SITES"
        printf 'min_meth\t%s\n' "$FILTER_MIN_METH"
        printf 'max_meth\t%s\n' "$FILTER_MAX_METH_LABEL"
        printf 'script_sha256\t%s\n' "$SCRIPT_SHA256"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"$tmp"; then
        rm -f "$tmp"
        echo "ERROR: failed to write QC configuration: $metadata" >&2
        return 1
    fi
    mv "$tmp" "$metadata"
}

# 在 filtered 目录中记录阈值、输入目录和过滤前后细胞数。
write_filter_provenance() {
    local filtered_dir="$1"
    local compact_dir="$2"
    local min_sites="$3"
    local metadata="$filtered_dir/filter_provenance.tsv"
    local tmp="${metadata}.tmp.$$"

    if ! {
        printf 'qc_tag\t%s\n' "$QC_TAG"
        printf 'min_sites\t%s\n' "$min_sites"
        printf 'max_sites\t%s\n' "$FILTER_MAX_SITES"
        printf 'min_meth\t%s\n' "$FILTER_MIN_METH"
        printf 'max_meth\t%s\n' "$FILTER_MAX_METH_LABEL"
        printf 'input_compact\t%s\n' "$compact_dir"
        printf 'cells_before\t%s\n' "$(count_cells "$compact_dir")"
        printf 'cells_after\t%s\n' "$(count_cells "$filtered_dir")"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"$tmp"; then
        rm -f "$tmp"
        echo "ERROR: failed to write filter provenance: $metadata" >&2
        return 1
    fi
    mv "$tmp" "$metadata"
}

# ==============================================================================
# 5. 日志管理与状态盘点
# ==============================================================================

# 新一轮运行前保留旧日志，并追加时间戳，避免覆盖诊断证据。
rotate_log() {
    local log="$1"
    if [[ -e "$log" ]]; then
        mv "$log" "${log}.previous.$(date +%Y%m%d_%H%M%S)"
    fi
}

# 执行单个 Methscan 步骤：记录日志、检查退出码、成功后写入 .ok 标记。
stage_number() {
    case "$1" in
        prepare) printf '2' ;;
        profile) printf '3' ;;
        filter) printf '4' ;;
        smooth) printf '5' ;;
        scan) printf '6' ;;
        matrix) printf '7' ;;
        *) printf '?' ;;
    esac
}

run_logged() {
    local step="$1"
    local log="$2"
    local ok_file="$3"
    local sequence
    shift 3

    sequence="$(stage_number "$step")"

    rotate_log "$log"
    rm -f "$ok_file"
    echo "    [$sequence/8 RUN] $step"

    if "$@" >"$log" 2>&1; then
        date -Is >"$ok_file"
        echo "    [$sequence/8 OK]  $step"
        return 0
    else
        local rc=$?
        echo "    [$sequence/8 FAIL] $step (exit $rc); see $log" >&2
        return "$rc"
    fi
}

# 输出一个样本在一个阈值下的产物数量。
status_one() {
    local sample_dir="$1"
    local threshold="$2"
    local min_sites="${threshold%k}000"
    local sample="${sample_dir##*/}"
    local compact=0 filtered=0 vmrs=0 matrix logs
    local compact_dir
    local qc_root
    qc_root="$(qc_root_for_sample "$sample_dir")"
    local filtered_dir="$qc_root/filtered_data_single_${threshold}"
    local scan_dir="$qc_root/scan_results_single_${threshold}"
    local matrix_dir="$qc_root/VMR_matrix_single_${threshold}"
    local log_dir="$qc_root/logs_single_${threshold}"

    if compact_dir="$(choose_compact "$sample_dir" "$threshold")"; then
        compact="$(count_cells "$compact_dir")"
    fi
    if valid_qc_config "$qc_root"; then
        if valid_filtered "$filtered_dir" "$min_sites"; then
            filtered="$(count_cells "$filtered_dir")"
        fi
        if valid_scan "$scan_dir"; then
            vmrs="$(wc -l <"$scan_dir/VMRs.bed")"
        fi
        matrix="$(count_files "$matrix_dir")"
        logs="$(count_files "$log_dir")"
    else
        matrix=0
        logs=0
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$threshold" "$compact" "$filtered" "$vmrs" "$matrix" "$logs"
}

# 输出指定阈值或全部阈值的状态表。
show_status() {
    local requested="${1:-}"
    local sample_dir threshold

    if [[ -n "$requested" ]]; then
        is_threshold "$requested" || die "invalid threshold: $requested"
    fi

    collect_samples
    printf '# qc_tag=%s min_meth=%s max_meth=%s max_sites=%s\n' \
        "$QC_TAG" "$FILTER_MIN_METH" "$FILTER_MAX_METH_LABEL" "$FILTER_MAX_SITES"
    printf 'sample\tthreshold\tcompact_cells\tfiltered_cells\tVMRs\tmatrix_files\tlogs\n'
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        if [[ -n "$requested" ]]; then
            status_one "$sample_dir" "$requested"
        else
            for threshold in "${VALID_THRESHOLDS[@]}"; do
                status_one "$sample_dir" "$threshold"
            done
        fi
    done
}

# ==============================================================================
# 6. 残缺产物保护
# ==============================================================================

# 此函数仅在产物未通过完整性校验时调用。
# 无论 .ok 是否存在，只要目录非空就拒绝覆盖并要求先归档。
refuse_untrusted_partial() {
    local dir="$1"
    local ok_file="$2"
    local label="$3"

    if [[ -d "$dir" ]] &&
        [[ -n "$(find "$dir" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "ERROR: invalid or unverified $label output exists: $dir" >&2
        if [[ -s "$ok_file" ]]; then
            echo "       A success marker exists, but the output failed validation." >&2
        fi
        echo "       Archive that directory, then rerun this sample." >&2
        return 1
    fi
}

# ==============================================================================
# 7. 单样本 Methscan 主流程
#    已存在完整 matrix 时整条跳过；否则逐步运行并验证。
# ==============================================================================

run_one_sample() {
    local sample_dir="$1"
    local threshold="$2"
    local threads="$3"
    local min_sites="${threshold%k}000"
    local sample="${sample_dir##*/}"
    local qc_root
    qc_root="$(qc_root_for_sample "$sample_dir")"
    local log_dir="$qc_root/logs_single_${threshold}"
    local filtered_dir="$qc_root/filtered_data_single_${threshold}"
    local scan_dir="$qc_root/scan_results_single_${threshold}"
    local matrix_dir="$qc_root/VMR_matrix_single_${threshold}"
    local compact_dir profile_file profile_metadata cov_dir
    local -a cov_files filter_args

    # 全部样本固定使用各自的概率去重 cov 目录。
    cov_dir="$sample_dir/$COV_SUBDIR"

    echo ">>> $sample $threshold"

    # ---------- [1/8] Preflight：QC 配置与整体完成检查 ----------
    echo "    [1/8 CHECK] QC configuration and existing outputs"
    ensure_qc_config "$qc_root" || return 1

    if compact_dir="$(choose_compact "$sample_dir" "$threshold")" &&
        valid_profile "$sample_dir" "$compact_dir" &&
        valid_filtered "$filtered_dir" "$min_sites" &&
        [[ -s "$log_dir/smooth.ok" ]] && valid_smooth "$filtered_dir" &&
        valid_scan "$scan_dir" &&
        valid_matrix "$matrix_dir"; then
        echo "    [2/8 SKIP] prepare/compact already exists"
        echo "    [3/8 SKIP] profile already exists"
        echo "    [4/8 SKIP] filter already exists"
        echo "    [5/8 SKIP] smooth already completed"
        echo "    [6/8 SKIP] scan already exists"
        echo "    [7/8 SKIP] matrix already exists"
        echo "<<< $sample $threshold complete"
        return 0
    fi

    mkdir -p "$log_dir"

    # ---------- [2/8] Prepare：优先复用 compact，无合格产物时才生成 ----------
    if ! compact_dir="$(choose_compact "$sample_dir" "$threshold")"; then
        compact_dir="$sample_dir/$COMPACT_SUBDIR"
        refuse_untrusted_partial "$compact_dir" "$log_dir/prepare.ok" prepare || return 1

        shopt -s nullglob
        cov_files=("$cov_dir"/*.cov.gz)
        shopt -u nullglob
        [[ "${#cov_files[@]}" -gt 0 ]] || {
            echo "ERROR: no *.cov.gz files for $sample in $cov_dir" >&2
            return 1
        }

        echo "    [2/8 INPUT] cov directory: $cov_dir (${#cov_files[@]} files)"
        mkdir -p "$compact_dir"
        run_logged prepare "$log_dir/prepare.log" "$log_dir/prepare.ok" \
            methscan prepare "${cov_files[@]}" "$compact_dir" || return 1
        valid_compact "$compact_dir" || {
            rm -f "$log_dir/prepare.ok"
            echo "ERROR: prepare exited successfully but compact output is invalid" >&2
            return 1
        }
    else
        echo "    [2/8 REUSE] prepare/compact: $compact_dir"
    fi

    if [[ "$STOP_AFTER_PREPARE" == 1 ]]; then
        echo "    [3/8 STOP] profile not requested"
        echo "    [4/8 STOP] filter not requested"
        echo "    [5/8 STOP] smooth not requested"
        echo "<<< $sample compact ready"
        return 0
    fi

    # ---------- [3/8] Profile：优先复用已有 TSS QC ----------
    if profile_file="$(choose_profile "$sample_dir" "$compact_dir")"; then
        echo "    [3/8 REUSE] profile: $profile_file"
    else
        profile_file="$sample_dir/${PROFILE_BASENAME}.csv"
        profile_metadata="$sample_dir/${PROFILE_BASENAME}.meta.tsv"

        # 旧 common profile 未提供参考文件校验信息，先改名保留再重建。
        if [[ -e "$profile_file" ]]; then
            mv "$profile_file" "${profile_file}.previous.$(date +%Y%m%d_%H%M%S)"
        fi
        if [[ -e "$profile_metadata" ]]; then
            mv "$profile_metadata" "${profile_metadata}.previous.$(date +%Y%m%d_%H%M%S)"
        fi

        run_logged profile "$log_dir/profile.log" "$log_dir/profile.ok" \
            methscan profile --strand-column 6 "$TSS_BED" \
            "$compact_dir" "$profile_file" || return 1
        [[ -s "$profile_file" ]] || {
            rm -f "$log_dir/profile.ok"
            echo "ERROR: profile exited successfully but CSV is empty" >&2
            return 1
        }

        if ! {
            printf 'tss_bed\t%s\n' "$TSS_BED"
            printf 'tss_sha256\t%s\n' "$TSS_SHA256"
            printf 'input_compact\t%s\n' "$compact_dir"
            printf 'created_at\t%s\n' "$(date -Is)"
        } >"${profile_metadata}.tmp.$$"; then
            rm -f "$log_dir/profile.ok" "${profile_metadata}.tmp.$$"
            echo "ERROR: failed to write profile provenance metadata" >&2
            return 1
        fi
        if ! mv "${profile_metadata}.tmp.$$" "$profile_metadata"; then
            rm -f "$log_dir/profile.ok" "${profile_metadata}.tmp.$$"
            echo "ERROR: failed to finalize profile provenance metadata" >&2
            return 1
        fi
    fi

    # ---------- [4/8] Filter：应用覆盖度和最低甲基化阈值 ----------
    if [[ -s "$log_dir/filter.ok" ]] && valid_filtered "$filtered_dir" "$min_sites"; then
        echo "    [4/8 SKIP] filter"
    else
        refuse_untrusted_partial "$filtered_dir" "$log_dir/filter.ok" filter || return 1
        mkdir -p "$filtered_dir"

        filter_args=(
            --min-sites "$min_sites"
            --max-sites "$FILTER_MAX_SITES"
            --min-meth "$FILTER_MIN_METH"
        )
        if [[ -n "$FILTER_MAX_METH" ]]; then
            filter_args+=(--max-meth "$FILTER_MAX_METH")
        fi

        run_logged filter "$log_dir/filter.log" "$log_dir/filter.ok" \
            methscan filter "${filter_args[@]}" \
            "$compact_dir" "$filtered_dir" || return 1
        write_filter_provenance "$filtered_dir" "$compact_dir" "$min_sites" || {
            rm -f "$log_dir/filter.ok"
            return 1
        }
        valid_filtered "$filtered_dir" "$min_sites" || {
            rm -f "$log_dir/filter.ok"
            echo "ERROR: filtered output or provenance is invalid for $sample" >&2
            return 1
        }
    fi

    # ---------- [5/8] Smooth：在 filtered 数据上原位平滑 ----------
    if [[ -s "$log_dir/smooth.ok" ]] && valid_smooth "$filtered_dir"; then
        echo "    [5/8 SKIP] smooth"
    else
        refuse_untrusted_partial "$filtered_dir/smoothed" "$log_dir/smooth.ok" smooth ||
            return 1
        run_logged smooth "$log_dir/smooth.log" "$log_dir/smooth.ok" \
            methscan smooth "$filtered_dir" || return 1
        valid_smooth "$filtered_dir" || {
            rm -f "$log_dir/smooth.ok"
            echo "ERROR: smooth output failed validation" >&2
            return 1
        }
    fi

    # MethSCAn diff/DMR 只需要 filtered + smoothed 数据。
    if [[ "$STOP_AFTER_SMOOTH" == 1 ]]; then
        echo "    [6/8 SKIP] scan not requested (DMR input mode)"
        echo "    [7/8 SKIP] matrix not requested (DMR input mode)"
        echo "<<< $sample $threshold DMR input ready"
        return 0
    fi

    # ---------- [6/8] Scan：发现 VMR ----------
    if [[ -s "$log_dir/scan.ok" ]] && valid_scan "$scan_dir"; then
        echo "    [6/8 SKIP] scan"
    else
        refuse_untrusted_partial "$scan_dir" "$log_dir/scan.ok" scan || return 1
        mkdir -p "$scan_dir"
        run_logged scan "$log_dir/scan.log" "$log_dir/scan.ok" \
            methscan scan --threads "$threads" "$filtered_dir" "$scan_dir/VMRs.bed" || return 1
        valid_scan "$scan_dir" || {
            rm -f "$log_dir/scan.ok"
            echo "ERROR: scan exited successfully but VMRs.bed is empty" >&2
            return 1
        }
    fi

    # ---------- [7/8] Matrix：生成每细胞 × VMR 矩阵 ----------
    if [[ -s "$log_dir/matrix.ok" ]] && valid_matrix "$matrix_dir"; then
        echo "    [7/8 SKIP] matrix"
    else
        refuse_untrusted_partial "$matrix_dir" "$log_dir/matrix.ok" matrix || return 1
        mkdir -p "$matrix_dir"
        run_logged matrix "$log_dir/matrix.log" "$log_dir/matrix.ok" \
            methscan matrix --threads "$threads" "$scan_dir/VMRs.bed" \
            "$filtered_dir" "$matrix_dir" || return 1
        valid_matrix "$matrix_dir" || {
            rm -f "$log_dir/matrix.ok"
            echo "ERROR: matrix exited successfully but expected files are missing" >&2
            return 1
        }
    fi

    echo "<<< $sample $threshold complete"
}

# ==============================================================================
# 8. 环境初始化、受控并发与失败汇总
# ==============================================================================

run_samples() {
    local threshold="$1"
    local max_jobs="$2"
    local threads="$3"
    local selected_sample="${4:-all}"
    local sample_dir
    local failures=0
    local i
    local -a pids=()
    local -a names=()

    # 仅 run 模式需要加载计算环境；status 模式保持只读且不依赖 Conda。
    echo "[1/8 CHECK] initialize environment and validate global inputs"
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization script not found: $CONDA_INIT"
    source "$CONDA_INIT" || die "failed to initialize Conda from $CONDA_INIT"
    conda activate "$CONDA_ENV" || die "failed to activate Conda environment: $CONDA_ENV"
    command -v methscan >/dev/null 2>&1 || die "methscan is not available in $CONDA_ENV"
    [[ -s "$TSS_BED" ]] || die "TSS BED not found or empty: $TSS_BED"
    command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required for TSS provenance checks"
    TSS_SHA256="$(sha256sum "$TSS_BED" | awk '{print $1}')"
    [[ "$TSS_SHA256" =~ ^[0-9a-f]{64}$ ]] || die "failed to calculate TSS BED SHA-256"
    SCRIPT_SHA256="$(sha256sum "${BASH_SOURCE[0]}" | awk '{print $1}')"
    [[ "$SCRIPT_SHA256" =~ ^[0-9a-f]{64}$ ]] || die "failed to calculate pipeline SHA-256"
    collect_samples

    # 指定样本时只运行该目录；all 则保留完整样本列表。
    if [[ "$selected_sample" != all ]]; then
        local selected_dir="$BASE_DIR/$selected_sample"
        [[ -d "$selected_dir" ]] || die "sample directory not found: $selected_dir"
        [[ "$selected_sample" == *_Met ]] || die "invalid sample name: $selected_sample"
        SAMPLE_DIRS=("$selected_dir")
    fi

    # 等待当前批次全部后台任务，并逐个汇总退出状态。
    wait_batch() {
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[8/8 SAMPLE OK] ${names[$i]}"
            else
                echo "[8/8 SAMPLE FAIL] ${names[$i]}" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    }

    # 按 max_jobs 分批启动，避免旧脚本一次性并发全部样本。
    echo "=== threshold=$threshold max_jobs=$max_jobs threads_per_job=$threads sample=$selected_sample ==="
    echo "=== qc_tag=$QC_TAG min_meth=$FILTER_MIN_METH max_meth=$FILTER_MAX_METH_LABEL max_sites=$FILTER_MAX_SITES ==="
    echo "=== tss_bed=$TSS_BED sha256=$TSS_SHA256 ==="
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        run_one_sample "$sample_dir" "$threshold" "$threads" &
        pids+=("$!")
        names+=("${sample_dir##*/}")

        if [[ "${#pids[@]}" -ge "$max_jobs" ]]; then
            wait_batch
        fi
    done
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch

    if [[ "$failures" -gt 0 ]]; then
        echo "[8/8 FAIL] $failures sample(s) failed" >&2
        return 1
    fi
    echo "[8/8 OK] ALL SAMPLES COMPLETE"
}

# ==============================================================================
# 9. 命令行入口
# ==============================================================================

main() {
    local action="${1:-}"
    local threshold="${2:-}"
    local max_jobs="${3:-$DEFAULT_MAX_JOBS}"
    local threads="${4:-$DEFAULT_THREADS}"
    local selected_sample="${5:-all}"

    validate_filter_config

    case "$action" in
        status)
            show_status "$threshold"
            ;;
        run|run-to-compact|run-to-smooth)
            [[ -n "$threshold" ]] || die "run requires a threshold"
            is_threshold "$threshold" || die "invalid threshold: $threshold"
            is_positive_integer "$max_jobs" || die "max_jobs must be a positive integer"
            is_positive_integer "$threads" || die "threads must be a positive integer"
            if [[ "$action" == run-to-compact ]]; then
                STOP_AFTER_PREPARE=1
            elif [[ "$action" == run-to-smooth ]]; then
                STOP_AFTER_SMOOTH=1
            fi
            run_samples "$threshold" "$max_jobs" "$threads" "$selected_sample"
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            usage >&2
            exit 1
            ;;
    esac
}

# 直接执行时进入主程序；被测试脚本 source 时只加载函数。
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
