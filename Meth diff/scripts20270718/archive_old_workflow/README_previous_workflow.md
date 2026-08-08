# MethScan DMR与threshold VMR运行结果

服务器目录：

```text
cd /share/home/rzli/METHSCAN/Meth_diff/20260716
```

`individual-effect mask`：由组内样本DMR组成的一组BED区间；任何与这些区间
重叠的VMR都会被整体排除。

主流程：

```text
输入
├── MethSCAn filtered data（16,241个细胞）
└── ALL_annotation_200k.csv（匹配15,500个细胞）
    ↓
03_meth_diff_celltype_sample_pairwise_200k.sh
├── 调用04_generate_celltype_sample_pairwise_groups.R
│   生成“同细胞类型、同response组内”的样本两两分组
└── 运行每种细胞类型内IR组和NR组的样本两两DMR
    ↓
08_merge_celltype_sample_pairwise_dmrs.sh
按cell type × response筛选并合并q<0.05 DMR
├── IR：14个*_sample_pairwise_union_q005.bed
└── NR：14个*_sample_pairwise_union_q005.bed
    ↓
10_prepare_individual_effect_mask.sh
├── 查找第8步IR/NR共28个并集BED
├── 记录source_files.tsv（28个来源，26个非空）
├── 从非空BED提取chr/start/end前三列
├── 合并得到76,193个原始区域
├── 按chr/start/end排序
└── bedtools merge合并重叠或相邻区域
    输出individual_effect_union_q005.bed（67,167个区域）
    ↓
11_run_threshold_vmrs_remove_individual.sh
├── threshold005：重新methscan scan（var-threshold=0.05）
├── threshold002：复用现有All VMR（var-threshold=0.02）
└── threshold001：重新methscan scan（var-threshold=0.01）
    ↓
bedtools将每套All VMR分别减去同一个individual-effect mask
Clean VMR = All VMR − 与individual_effect_union_q005.bed任意重叠的VMR
├── bedtools intersect -v：保留不重叠VMR → clean_VMRs.bed
└── bedtools intersect -u：记录被删除VMR → removed_individual_effect_VMRs.bed
    ├── threshold005：159,457 − 26,327 = 133,130 Clean VMRs
    ├── threshold002： 83,245 − 13,135 =  70,110 Clean VMRs
    └── threshold001： 45,119 −  7,252 =  37,867 Clean VMRs
    ↓
methscan matrix
为每套Clean VMR生成16,241个细胞的甲基化矩阵
    ↓
validate_threshold_matrix.py
验证matrix细胞数、VMR数量及Clean VMR坐标集合
    ↓
13_run_threshold_clean_vmr_reclustering.sh
└── 13_recluster_threshold_clean_vmrs.R
    ├── 匹配15,500个Scanpy注释细胞
    ├── 删除全NA/零方差VMR
    ├── iterative PCA（20 PCs）
    ├── UMAP（30 neighbors，min_dist=0.05，seed=2）
    ├── Leiden（resolution=0.001）
    └── 输出cell type、sample、response、Leiden图及比较指标
        ↓
14_collect_threshold_metrics.sh
合并三套comparison_metrics.tsv
输出comparison_metrics_all.tsv并比较最终threshold
```

`01/02`、`05/06`是补充DMR分析；`07/09/12`属于旧流程，不进入上述主链。

第10–14步的矩阵和聚类统计接在本文第5–9节。

## 1. 输入与细胞统计

| 项目 | 统计 |
|---|---:|
| MethSCAn filtered data细胞 | 16,241 |
| 匹配Scanpy注释 | 15,500 |
| 未匹配 | 741 |
| 正式DMR有效细胞 | 15,466 |
| IR有效细胞 | 8,227 |
| NR有效细胞 | 7,239 |
| 排除污染细胞 | 34 |
| `MIN_CELLS` | 6 |

样本有效细胞数：

| 样本 | 细胞数 | 样本 | 细胞数 |
|---|---:|---|---:|
| IR01 | 1,523 | NR01 | 917 |
| IR02 | 2,012 | NR02 | 1,981 |
| IR03 | 2,173 | NR03 | 1,556 |
| IR04 | 1,132 | NR04 | 1,785 |
| IR05 | 1,387 | NR05 | 1,000 |

## 2. 脚本运行结果

| 脚本 | 运行结果 | 状态 |
|---|---|---|
| `01_meth_diff_pairwise_200k.sh` | 15种细胞类型IR vs NR；Job 162782 | 完成 |
| `02_generate_cell_groups.R` | 生成15套IR/NR cell-group | 完成 |
| `03_meth_diff_celltype_sample_pairwise_200k.sh` | 每种细胞类型内IR、NR样本两两DMR，理论最多300组 | 完成 |
| `04_generate_celltype_sample_pairwise_groups.R` | 生成cell type × response × sample pair分组 | 完成 |
| `05_meth_diff_sample_pairwise_200k.sh` | IR 10组、NR 10组；Job 162785 | 完成 |
| `06_generate_sample_pairwise_groups.R` | 生成20套全细胞样本两两分组 | 完成 |
| `07_merge_celltype_ir_vs_nr_dmrs.sh` | 合并15种细胞类型IR vs NR DMR；新主流程不使用 | 完成 |
| `08_merge_celltype_sample_pairwise_dmrs.sh` | IR 14个并集、NR 14个并集 | 完成 |
| `10_prepare_individual_effect_mask.sh` | 28个来源文件合并为67,167个mask regions | 完成 |
| `11_run_threshold_vmrs_remove_individual.sh` | 0.05、0.02、0.01三套Clean VMR矩阵 | 完成 |
| `13_run_threshold_clean_vmr_reclustering.sh` | 0.01、0.02完成；0.05运行中 | 进行中 |
| `14_collect_threshold_metrics.sh` | 等0.05完成后汇总 | 待运行 |

## 3. 上游DMR结果

### 3.1 细胞类型内IR vs NR

```text
Job ID：162782
比较数：15
结果目录：result/DMR_results_200k/3_same_cell_type_IR_vs_NR
```

### 3.2 细胞类型内组内样本比较

第8步合并结果：

| response | 并集BED数 | 非空BED数 |
|---|---:|---:|
| IR | 14 | 14 |
| NR | 14 | 12 |
| 合计 | 28 | 26 |

结果目录：

```text
result/celltype_sample_pairwise/merged_DMRs_200k/q005/IR
result/celltype_sample_pairwise/merged_DMRs_200k/q005/NR
```

这28个BED是`individual-effect mask`的正式输入。

### 3.3 不区分细胞类型的组内样本比较

```text
Job ID：162785
IR比较数：10
NR比较数：10
总比较数：20
```

| response | DMR总数 | 单个比较最少 | 单个比较最多 |
|---|---:|---:|---:|
| IR | 801,366 | 70,773 | 87,058 |
| NR | 828,608 | 75,520 | 88,157 |

结果目录：

```text
result/sample_pairwise/DMR_results_200k/5_IR_sample_pairwise
result/sample_pairwise/DMR_results_200k/5_NR_sample_pairwise
```

该结果不进入当前`individual-effect mask`。

## 4. 上下游交接结果

`individual_effect_union_q005.bed`的生成过程：

```text
第8步IR 14个并集BED + NR 14个并集BED
→ 排除2个空文件
→ 从26个非空BED提取前三列
→ 拼接为76,193个区域
→ sort -k1,1V -k2,2n -k3,3n
→ bedtools merge
→ 67,167个非重叠mask regions
```

```text
第8步来源BED：28
非空来源BED：26
来源区域：76,193
合并后mask regions：67,167
```

第10步输出：

```text
result/individual_effect_mask/individual_effect_union_q005.bed
result/individual_effect_mask/source_files.tsv
result/individual_effect_mask/mask_summary.tsv
```

## 5. Threshold矩阵统计

三套均使用同一个shared `individual-effect mask`，并执行：

```text
Clean VMRs = All VMRs − removed individual-effect VMRs
```

| variant | threshold | All VMRs | Clean VMRs | 删除VMRs | 删除比例 | matrix cells |
|---|---:|---:|---:|---:|---:|---:|
| `threshold005` | 0.05 | 159,457 | 133,130 | 26,327 | 16.51% | 16,241 |
| `threshold002` | 0.02 | 83,245 | 70,110 | 13,135 | 15.78% | 16,241 |
| `threshold001` | 0.01 | 45,119 | 37,867 | 7,252 | 16.07% | 16,241 |

矩阵验证：

| variant | matrix VMRs | missing VMRs | extra VMRs | BED顺序一致 | 状态 |
|---|---:|---:|---:|---|---|
| `threshold005` | 133,130 | 0 | 0 | 否 | 通过 |
| `threshold002` | 70,110 | 0 | 0 | 否 | 通过 |
| `threshold001` | 37,867 | 0 | 0 | 否 | 通过 |

`BED顺序一致=否`只表示MethSCAn重新排列矩阵列；VMR集合完全一致。

矩阵Job记录：

| Job ID | variant | 结果 |
|---:|---|---|
| 162808 | threshold005 | scan和matrix完成；旧validator因列顺序退出 |
| 162810 | threshold005 | 更新validator后验证成功 |
| 162811 | threshold002 | matrix成功，峰值内存约39.6G |
| 162813 | threshold001 | 已有完整结果，验证后跳过 |

## 6. 重聚类统计

统一分析参数：

```text
矩阵细胞：16,241
匹配Scanpy注释：15,500
PCA：20 PCs
UMAP：30 neighbors，min_dist 0.05，seed 2
Leiden resolution：0.001
```

运行记录：

| variant | Job ID | OpenBLAS线程 | 运行时间 | VMRs | Leiden簇 | 状态 |
|---|---:|---:|---:|---:|---:|---|
| `threshold001` | 162815 | 32 | 41分25秒 | 37,867 | 59 | COMPLETE |
| `threshold002` | 162818 | 64 | 1时14分45秒 | 70,110 | 42 | COMPLETE |
| `threshold005` | 162819 | 32 | 待完成 | 133,130 | 待完成 | RUNNING |

特征检查：

| variant | all-NA VMRs | zero-variance VMRs | 缺失注释细胞 |
|---|---:|---:|---:|
| `threshold001` | 0 | 0 | 0 |
| `threshold002` | 0 | 0 | 0 |
| `threshold005` | 待完成 | 待完成 | 待完成 |

## 7. 聚类指标

原始200k参考：

| cell-type purity | sample purity | sample mixing entropy |
|---:|---:|---:|
| 0.7203 | 0.5154 | 0.5688 |

Threshold结果：

| variant | cell-type purity | sample purity | sample mixing entropy | response purity | response entropy | ARI vs 200k |
|---|---:|---:|---:|---:|---:|---:|
| `threshold001` | 0.7126 | 0.3005 | 0.8443 | 0.6094 | 0.9167 | 0.4033 |
| `threshold002` | 0.7135 | 0.3369 | 0.8090 | 0.6312 | 0.8910 | 0.4932 |
| `threshold005` | 待完成 | 待完成 | 待完成 | 待完成 | 待完成 | 待完成 |

相对参考变化：

| variant | Δ cell-type purity | Δ sample purity | Δ sample entropy |
|---|---:|---:|---:|
| `threshold001` | -0.0077 | -0.2149 | +0.2754 |
| `threshold002` | -0.0068 | -0.1785 | +0.2402 |
| `threshold005` | 待完成 | 待完成 | 待完成 |

当前结果：

```text
threshold001：样本效应去除最强
threshold002：更接近原始200k聚类
threshold005：等待Job 162819
```

## 8. 输出文件

每套矩阵：

```text
result/threshold_VMR_remove_individual/<variant>/
├── all_VMRs.bed
├── clean_VMRs.bed
├── removed_individual_effect_VMRs.bed
├── VMR_matrix/mean_shrunken_residuals.csv.gz
├── matrix_validation.tsv
├── run_metadata.tsv
└── .complete
```

每套聚类：

```text
result/threshold_VMR_remove_individual_reclustering/<variant>/
├── feature_qc.tsv
├── feature_qc_summary.tsv
├── comparison_metrics.tsv
├── plots/
└── .complete
```

旧的threshold005单线程部分结果：

```text
result/threshold_VMR_remove_individual_reclustering/
  threshold005.singlethread_partial_162814/
```

## 9. 最终汇总

Job 162819完成后运行：

```bash
bash 14_collect_threshold_metrics.sh
```

汇总文件：

```text
result/threshold_VMR_remove_individual_reclustering/comparison_metrics_all.tsv
```
