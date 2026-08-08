# METHSCAN 10样本独立分析流程

本目录顶层`01–08`为当前正式流程。分析单位始终是单个样本：10个样本可并行运行，但不合并compact、不构建联合pseudobulk，也不执行跨样本DMR。

完整命令和参数见[`00_ALL_SAMPLES_WORKFLOW.md`](00_ALL_SAMPLES_WORKFLOW.md)，cov去重统计见[`Report.md`](Report.md)。

## 数据流

```text
<sample>/cov
  → <sample>/cov_dedup_probability
  → <sample>/compact_data_dedup_probability
  → <sample>/qc_..._covdedupprob/filtered_data_single_300k
  → <sample>/qc_..._covdedupprob/methdiff_celltype_300k
  → 每个样本自己的Top200、matrix和可选热图
```

## 必需的远程文件

```text
/share/home/rzli/METHSCAN/01_Upstream/
├── 00_ALL_SAMPLES_WORKFLOW.md
├── 01_check_cov_duplicates.sh
├── 02_deduplicate_cov_by_probability.sh
├── 03_run_upstream_pipeline.sh
├── 04_run_all_samples_to_smooth.sh
├── 05_run_all_samples_dmr.sh
├── 06_select_top200_dmrs.sh
├── 07_compute_top200_dmr_matrix.sh
├── 08_plot_top200_dmr_heatmap.sh
└── lib/
    ├── check_cov_duplicates_one_sample.sh
    └── deduplicate_cov_by_probability_one_sample.sh

/share/home/rzli/METHSCAN/02_Methdiff/
├── run_single_sample_dmr.sh
└── Result/
    ├── 02_merge_sample_dmrs.py
    ├── 04_plot_single_cell_dmr_heatmaps.py
    ├── 05_extract_celltype_hypo_dmrs_top1500.py
    └── 06_compute_dmr_mean_of_cpg_ratios.py
```

旧的`04_merge_all_samples_to_smooth.sh`和`run_ir01_single_sample_dmr.sh`已被替换，不应与新脚本同时保留在正式流程中。

## 上线前检查

```bash
cd /share/home/rzli/METHSCAN/01_Upstream

chmod 750 0[1-8]_*.sh lib/*.sh \
  /share/home/rzli/METHSCAN/02_Methdiff/run_single_sample_dmr.sh

for script in 0[1-8]_*.sh lib/*.sh \
  /share/home/rzli/METHSCAN/02_Methdiff/run_single_sample_dmr.sh
do
  echo "CHECK $script"
  bash -n "$script" || exit 1
done

python -m py_compile \
  /share/home/rzli/METHSCAN/02_Methdiff/Result/02_merge_sample_dmrs.py \
  /share/home/rzli/METHSCAN/02_Methdiff/Result/04_plot_single_cell_dmr_heatmaps.py \
  /share/home/rzli/METHSCAN/02_Methdiff/Result/05_extract_celltype_hypo_dmrs_top1500.py \
  /share/home/rzli/METHSCAN/02_Methdiff/Result/06_compute_dmr_mean_of_cpg_ratios.py
```

## 已废弃的联合目录

当前流程不读取下列目录；若存在，只作为历史记录归档：

```text
/share/LCZX_Data/data/allcools/merged_10samples_covdedupprob*
```
