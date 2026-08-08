# MethSCAn 单样本内细胞类型DMR

当前正式入口是`run_single_sample_dmr.sh`。它读取一个样本自己的`filtered_data_single_300k`及其smooth结果，只在该样本内部进行细胞类型两两比较。

通常不直接逐样本调用，而由：

```text
../01_Upstream/05_run_all_samples_dmr.sh
```

统一调度10个样本。

## 单样本直接调用

```bash
SAMPLE_NAME=25110891_IR02_Met THRESHOLD=300k \
bash run_single_sample_dmr.sh prepare

SAMPLE_NAME=25110891_IR02_Met THRESHOLD=300k \
bash run_single_sample_dmr.sh run 2 24

SAMPLE_NAME=25110891_IR02_Met THRESHOLD=300k \
bash run_single_sample_dmr.sh status
```

## 固定规则

- `min_cells=10`
- 只比较同一样本内不同细胞类型
- 只向MethSCAn diff暴露`chr1–chr22, chrX, chrY`
- 空主染色体和所有非主染色体contig均排除
- 每个比较输出12列DMR BED、完成标记和显著性汇总
- 默认注释：`/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv`

## 旧联合流程

旧的10样本合并DMR脚本仅作历史复核，位于：

```text
archive_merged_workflow/run_methdiff_pipeline_merged_legacy.sh
```

当前正式流程不调用它，也不读取`merged_10samples_covdedupprob`。
