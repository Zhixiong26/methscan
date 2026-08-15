# 内置 Methdiff 实现

该目录由原`02_Methdiff`合并而来，是`01_Upstream`已验证单样本流程的内部实现：

```text
run_single_sample_dmr.sh
python/
├── 02_merge_sample_dmrs.py
├── 04_plot_single_cell_dmr_heatmaps.py
├── 05_extract_celltype_hypo_dmrs_top1500.py
└── 06_compute_dmr_mean_of_cpg_ratios.py
```

正常情况下不直接运行这些文件，而是使用顶层`05–08`入口。合并只改变代码位置和调用路径，不改变算法、参数、数据路径或输出目录。
