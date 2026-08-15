# METHSCAN 10样本独立分析流程

本目录只保留已运行成功的`01–09`单样本流程，并在`lib/methdiff/`内集成原`02_Methdiff`的必要实现。分析单位始终是单个样本；10个样本可并行，但不合并compact、不构建联合pseudobulk、不执行跨样本DMR。

主要数据流：

```text
cov
  → cov_dedup_probability
  → compact_data_dedup_probability
  → filtered_coverage_single_300k (30万–120万 coverage)
  → filtered_data_single_300k (Scanpy clean singlets)
  → smoothed
  → methdiff_celltype_300k
  → Top200 → single-cell DMR matrix → 8类热图
```

注释与 smooth 前 clean-cell 名单分别固定为：

```text
/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_all_cells.csv
/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_clean_cells.csv
```

统一热图入口：

```bash
bash 08_plot_all_top200_heatmaps.sh all
bash 08_plot_all_top200_heatmaps.sh status
bash 08_plot_all_top200_heatmaps.sh links
```

其中`DMRwise_Zscore_ColorClip1`保留未截断Z-score，仅把颜色范围设为`[-1,1]`。所有热图都排除“在该样本中没有自身特异DMR”的cell type细胞行。

完整参数与命令见[`00_ALL_SAMPLES_WORKFLOW.md`](00_ALL_SAMPLES_WORKFLOW.md)，cov去重统计见[`Report.md`](Report.md)。

未经本轮完整验证的联合流程、旧包装脚本和旧结果已转移到项目根目录的`archive_local/not_in_recent_success_20260815/`。
