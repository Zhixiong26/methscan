# VMR 聚类去除个体效应流程（LEGACY：二次截取 top 5%/2%/1%）

## 前置结果

第 1–6 步已经生成：

- 每种细胞类型的 IR vs NR DMR；
- 每种细胞类型内 IR、NR 样本间两两比较的 DMR。

本 README 从第 7 步开始，不再重复 DMR 计算过程。

## DMR 比较对象

本流程中的 DMR 均在**同一细胞类型内部**计算。以 `CD14_Monocytes` 为例：

### A：细胞类型内 IR vs NR DMR

```text
group_A = IR01–IR05 中所有有效的 CD14_Monocytes
group_B = NR01–NR05 中所有有效的 CD14_Monocytes
```

每种细胞类型产生 1 个 IR vs NR 比较，15 种细胞类型共 15 个比较。

### B：IR 组内同细胞类型的样本间 DMR

```text
IR01的CD14_Monocytes vs IR02的CD14_Monocytes
IR01的CD14_Monocytes vs IR03的CD14_Monocytes
...
IR04的CD14_Monocytes vs IR05的CD14_Monocytes
```

每种细胞类型理论上有 `choose(5,2)=10` 个 IR 组内样本比较。

### C：NR 组内同细胞类型的样本间 DMR

```text
NR01的CD14_Monocytes vs NR02的CD14_Monocytes
NR01的CD14_Monocytes vs NR03的CD14_Monocytes
...
NR04的CD14_Monocytes vs NR05的CD14_Monocytes
```

每种细胞类型理论上有 `choose(5,2)=10` 个 NR 组内样本比较。

任一比较组少于 6 个细胞时，该样本组合不运行。未匹配注释的细胞和 `Platelet_erythroid_contamination` 不参与比较。

本流程**不使用**忽略细胞类型的全细胞样本比较结果。

最终集合关系：

```text
Clean DMR = A − union(B) − union(C)
```

## 第 7 步：合并每种细胞类型的 IR vs NR DMR

脚本：

```text
07_merge_celltype_ir_vs_nr_dmrs.sh
```

输入来自第 1、2 步。每个文件表示一种细胞类型中，合并后的 IR01–IR05 细胞与合并后的 NR01–NR05 细胞之间的 DMR：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/DMR_results_200k/3_same_cell_type_IR_vs_NR/<cell_type>_IR_vs_NR_DMRs.bed
```

第 7 步对每种细胞类型的 DMR 区域分别执行 `bedtools merge`。这一步：

- 不筛选 q 值；
- 保留该比较中的全部 DMR；
- 只合并同一个细胞类型文件内发生重叠或相邻的 DMR；
- 不把不同细胞类型合并在一起。

输出：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/merged_celltype_IR_vs_NR_200k/all_dmr_union/<cell_type>_IR_vs_NR_union_all.bed
```

每种细胞类型合并前后的 DMR 数记录在：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/merged_celltype_IR_vs_NR_200k/merge_summary.tsv
```

## 第 8 步：合并 IR、NR 组内样本间 q<0.05 DMR

脚本：

```text
08_merge_celltype_sample_pairwise_dmrs.sh
```

输入来自第 3、4 步：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/DMR_results_200k/6_IR_within_celltype_sample_pairwise/
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/DMR_results_200k/6_NR_within_celltype_sample_pairwise/
```

理论上共有 300 个样本组合，实际完成 214 个比较结果 BED：

| 分组 | 可运行比较 | 非空 DMR BED | 0 DMR BED |
|---|---:|---:|---:|
| IR 组内同细胞类型样本比较 | 116 | 110 | 6 |
| NR 组内同细胞类型样本比较 | 98 | 96 | 2 |
| 合计 | 214 | 206 | 8 |

另外 86 个组合因任一侧少于 6 个细胞而未运行。214 是比较结果 BED 数量，不是 DMR 区域数。

对每种细胞类型分别执行：

```text
IR并集 = union(该细胞类型所有IR样本对的q<0.05 DMR)
NR并集 = union(该细胞类型所有NR样本对的q<0.05 DMR)
```

例如，`CD14_Monocytes` 的 `IR01_vs_IR02`、`IR01_vs_IR03`，直到 `IR04_vs_IR05` 的 q<0.05 DMR 合并为一个 IR 并集；NR 单独生成 NR 并集。IR 与 NR 不混合，不同细胞类型也不混合。

IR 输出：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/merged_DMRs_200k/q005/IR/<cell_type>_IR_sample_pairwise_union_q005.bed
```

NR 输出：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/merged_DMRs_200k/q005/NR/<cell_type>_NR_sample_pairwise_union_q005.bed
```

合并统计：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/merged_DMRs_200k/merge_summary.tsv
```

## 第 9 步：相减得到 Clean DMR

脚本：

```text
09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh
```

对每种细胞类型执行：

```text
Clean DMR
= 第7步得到的该细胞类型全部 IR vs NR DMR 并集
- 第8步得到的该细胞类型 IR 组内样本间 q<0.05 DMR 并集
- 第8步得到的该细胞类型 NR 组内样本间 q<0.05 DMR 并集
```

第 9 步使用 `bedtools intersect -v`：第 7 步的一个 IR vs NR DMR 只要与 IR 或 NR 扣除并集发生重叠，就从结果中移除。

Clean DMR 输出：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/clean_celltype_IR_vs_NR/clean/<cell_type>_IR_vs_NR_clean.bed
```

扣除前、扣除后及移除数量记录在：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/clean_celltype_IR_vs_NR/subtract_summary.tsv
```

### 实际运行结果

作业 `162802` 成功完成，第 7–9 步退出状态为 0。15 种细胞类型共得到：

```text
第7步 IR vs NR 合并区域：820,175
第9步移除区域：            13,912
最终 Clean DMR：           806,263
总体移除比例：               1.696%
```

各细胞类型结果：

| 细胞类型 | 第7步 IR vs NR 合并区域 | 第8步 IR+NR 扣除掩膜区域 | Clean DMR | 移除 DMR |
|---|---:|---:|---:|---:|
| B_cells | 92,458 | 244 | 92,413 | 45 |
| B_cells_unresolved | 9,008 | 2 | 9,007 | 1 |
| CD14_Monocytes | 86,270 | 21,424 | 79,998 | 6,272 |
| CD16_Monocytes | 90,049 | 38 | 90,035 | 14 |
| CD4_T_cells | 88,056 | 51,271 | 80,577 | 7,479 |
| CD8_T_cells | 92,252 | 215 | 92,186 | 66 |
| cDCs | 17 | 0 | 17 | 0 |
| Cycling_cells | 28,918 | 10 | 28,917 | 1 |
| Gamma_delta_T_cells | 10,854 | 3 | 10,853 | 1 |
| HLAII_high_APCs | 56,176 | 24 | 56,172 | 4 |
| MAIT_cells | 63,761 | 4 | 63,759 | 2 |
| NK_cells | 86,365 | 102 | 86,348 | 17 |
| pDCs | 6,940 | 15 | 6,937 | 3 |
| Plasma_cells | 56,418 | 13 | 56,415 | 3 |
| Treg_cells | 52,633 | 39 | 52,629 | 4 |

状态说明：

- `B_cells_unresolved` 和 `Cycling_cells` 的 NR q<0.05 并集为 0，IR 并集分别为 2 和 10；因此仍然使用有效的 IR 掩膜完成扣除。
- `cDCs` 没有满足细胞数要求的组内样本比较，因此没有 IR 或 NR 扣除掩膜。其 17 个 IR vs NR 区域全部保留，但这表示“无法评估组内个体差异”，不能解释为“确认不存在个体差异”。
- MethSCAn 非 debug DMR BED 的第 12 列为 `adjusted_p`，第 8 步使用该列执行 q<0.05 筛选。

本流程不使用忽略细胞类型的全细胞样本比较结果。

## 第 10 步：映射回 All VMR

脚本：

```text
10_map_clean_dmrs_to_all_vmrs.sh
```

将第 9 步每种细胞类型的 Clean DMR 映射到统一的 All VMR 区域，得到可从原始 VMR 矩阵中提取的 VMR ID 和坐标列表：

```text
Clean VMR
= All VMR 中与 Clean DMR 重叠的 VMR
```

Clean DMR 输入：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/clean_celltype_IR_vs_NR/clean/<cell_type>_IR_vs_NR_clean.bed
```

All VMR 坐标及矩阵 ID：

```text
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.bed
```

脚本先为每种细胞类型分别映射，再对所有细胞类型映射到的 VMR 去重。不同细胞类型可以映射到同一个 All VMR，最终并集中只保留一次。

每种细胞类型的映射结果：

```text
result/clean_celltype_IR_vs_NR/clean_VMRs/by_cell_type/<cell_type>_clean_VMRs.bed
```

所有细胞类型的 Clean VMR 去重并集：

```text
result/clean_celltype_IR_vs_NR/clean_VMRs/all_celltypes_clean_VMRs.bed
```

用于提取原始 VMR 矩阵列的 VMR ID：

```text
result/clean_celltype_IR_vs_NR/clean_VMRs/all_celltypes_clean_VMR_IDs.txt
```

统计文件：

```text
result/clean_celltype_IR_vs_NR/clean_VMRs/map_summary.tsv
result/clean_celltype_IR_vs_NR/clean_VMRs/union_summary.tsv
```

作业 `162803` 的实际结果：

```text
All VMR：          82,257
Clean VMR 并集：  60,284
Clean VMR ID：    60,284
```

Clean VMR 占标准染色体 All VMR 的约 73.29%。Clean VMR BED 行数与去重后的 VMR ID 数完全一致。

第 9 步的 806,263 个 Clean DMR 最终映射到 60,284 个不同的 All VMR。Clean DMR 数量大于 Clean VMR 数量，是因为不同细胞类型、不同坐标的多个 Clean DMR 可以重叠到同一个 All VMR；最终按 VMR ID 合并去重。

## 第 11 步：选择 top 5%、2%、1% VMR

脚本：

```text
11_select_top_clean_vmrs.sh
```

对第 10 步得到的 Clean VMR 去重并集进行排序。这里的 MethScan 差异强度定义为 All 数据运行 `methscan scan` 时输出的 VMR methylation-variance score，即 `VMRs.bed` 第 4 列：

```text
/share/LCZX_Data/data/All/scan_results/VMRs.bed
```

排序范围是第 10 步得到的 Clean VMR 集合，不是全部 All VMR。按照 score 从高到低分别保留：

```text
top 5%
top 2%
top 1%
```

数量按向上取整计算，因此：

```text
top 1% ⊆ top 2% ⊆ top 5% ⊆ 全部 Clean VMR
```

带 score 的完整排序表：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMRs/all_clean_VMRs_ranked.tsv
```

top BED：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top5pct.bed
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top2pct.bed
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top1pct.bed
```

用于提取矩阵列的 VMR ID：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top5pct_IDs.txt
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top2pct_IDs.txt
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top1pct_IDs.txt
```

如果任一 Clean VMR 无法与 `VMRs.bed` 中的坐标精确匹配，脚本会报错并停止，避免在排序时静默丢失区域。

作业 `162803` 的实际结果：

| 集合 | 比例 | VMR 数量 |
|---|---:|---:|
| 全部 Clean VMR | 100% | 60,284 |
| top 5% | 5% | 3,015 |
| top 2% | 2% | 1,206 |
| top 1% | 1% | 603 |

三个 top 集合均从同一份降序表的开头提取，因此严格满足：

```text
top 1% ⊆ top 2% ⊆ top 5%
```

第 10、11 步的完整数量关系：

```text
82,257 个 All VMR
        ↓ 与全部细胞类型的 Clean DMR 重叠并按 VMR ID 去重
60,284 个 Clean VMR
        ↓ 按 MethScan scan VMR score 从高到低排序
├── top 5%：ceil(60,284 × 0.05) = 3,015 个 VMR
├── top 2%：ceil(60,284 × 0.02) = 1,206 个 VMR
└── top 1%：ceil(60,284 × 0.01) =   603 个 VMR
```

因此：

```text
top 1%（603）
⊆ top 2%（1,206）
⊆ top 5%（3,015）
⊆ 全部 Clean VMR（60,284）
⊆ All VMR（82,257）
```

## 第 12 步：提取 VMR 矩阵

脚本：

```text
12_subset_top_clean_vmr_matrices.py
12_run_subset_top_clean_vmr_matrices.sh
```

原始矩阵：

```text
/share/LCZX_Data/data/All/VMR_matrix/mean_shrunken_residuals.csv.gz
```

分别使用第 11 步生成的 top 5%、2%、1% VMR ID，从原始矩阵提取对应列。Python 脚本只读取一次大型原始矩阵，并同时写出三套 gzip 压缩子矩阵，保留原始细胞顺序和原始 VMR 列顺序。

输入 VMR ID：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top5pct_IDs.txt
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top2pct_IDs.txt
result/clean_celltype_IR_vs_NR/top_clean_VMRs/clean_VMRs_top1pct_IDs.txt
```

输出：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMR_matrices/top5/mean_shrunken_residuals.csv.gz
result/clean_celltype_IR_vs_NR/top_clean_VMR_matrices/top2/mean_shrunken_residuals.csv.gz
result/clean_celltype_IR_vs_NR/top_clean_VMR_matrices/top1/mean_shrunken_residuals.csv.gz
```

预期矩阵维度：

```text
原始矩阵全部细胞 × 3,015 个 top 5% Clean VMR
原始矩阵全部细胞 × 1,206 个 top 2% Clean VMR
原始矩阵全部细胞 ×   603 个 top 1% Clean VMR
```

矩阵提取统计：

```text
result/clean_celltype_IR_vs_NR/top_clean_VMR_matrices/subset_matrix_summary.tsv
```

脚本会检查：

- top 1% ID 必须是 top 2% 的子集；
- top 2% ID 必须是 top 5% 的子集；
- 所有 VMR ID 必须存在于原始矩阵；
- 每一行的列数必须与矩阵表头一致。

任一检查失败时，脚本立即停止。

作业 `162804` 的实际矩阵提取结果：

| 矩阵 | 原始矩阵细胞数 | VMR 数 |
|---|---:|---:|
| top 5% | 16,241 | 3,015 |
| top 2% | 16,241 | 1,206 |
| top 1% | 16,241 | 603 |

三个 VMR ID 集合均完整匹配原始矩阵，没有缺失特征。

## 第 13 步：重新聚类、注释和比较

脚本：

```text
13_recluster_top_clean_vmrs.R
13_run_top_clean_vmr_reclustering.sh
```

顺序运行第 12、13 步：

```text
run_steps_12_13.sh
```

聚类输入是第 12 步的三套 `mean_shrunken_residuals.csv.gz`，不是 DMR BED。

每套矩阵分别独立运行：

```text
top 5% VMR矩阵 → PCA → UMAP → Leiden → 注释与比较
top 2% VMR矩阵 → PCA → UMAP → Leiden → 注释与比较
top 1% VMR矩阵 → PCA → UMAP → Leiden → 注释与比较
```

### 与原始 200k 分析保持一致的口径

第 13 步复用以下脚本的分析方法和参数：

```text
/share/home/rzli/METHSCAN/Annotation/20260716/02_All_200k_analysis.R
```

使用的 Scanpy 细胞类型注释：

```text
/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv
```

用于比较的原始 200k 聚类结果：

```text
/share/home/rzli/METHSCAN/Annotation/20260716/result/ALL_annotation_200k.csv
```

三个 top 比例与原始分析使用相同的：

```text
Scanpy匹配细胞：约15,500
PCA：20 PCs
UMAP n_neighbors：30
UMAP min_dist：0.05
UMAP seed：2
UMAP SGD：1线程
Leiden resolution：0.001
```

为了公平比较，三个 top 比例使用相同随机种子和相同参数。原始矩阵中没有 Scanpy 注释的细胞在 PCA 前排除；`Platelet_erythroid_contamination` 与原始 200k 注释流程一致，仍保留在聚类中。

这里的“注释”是将已有的 `cell_type_integrated` 标签映射到新 UMAP 和新 Leiden cluster，并根据每个新 cluster 的多数细胞类型生成 cluster 注释。它不是重新根据 marker 从零定义细胞类型。

### 每个 top 比例的输出

结果目录：

```text
result/clean_celltype_IR_vs_NR/reclustering/top5/
result/clean_celltype_IR_vs_NR/reclustering/top2/
result/clean_celltype_IR_vs_NR/reclustering/top1/
```

每个目录包含：

```text
<top>_PCA.RData
<top>_PCA_coordinates.csv
<top>_UMAP_coordinates.csv
<top>_annotation.csv
cluster_majority_annotation.csv
cluster_cell_type_composition.csv
cluster_sample_composition.csv
cluster_response_composition.csv
cell_cluster_comparison_to_original_200k.csv
comparison_metrics.tsv
run_metadata.tsv
plots/
```

三个比例的汇总比较表：

```text
result/clean_celltype_IR_vs_NR/reclustering/comparison_metrics_all.tsv
```

### 比较指标

`comparison_metrics_all.tsv` 记录：

- `cell_type_cluster_purity`：新 cluster 的细胞类型纯度，越高表示主要细胞类型结构保留越好；
- `sample_cluster_purity`：新 cluster 被单一样本主导的程度，越低越好；
- `sample_mixing_entropy`：样本混合度，越高表示不同样本混合越充分；
- `response_cluster_purity` 和 `response_mixing_entropy`：IR/NR 在新 cluster 中的组成；
- `ARI_vs_original_200k_leiden`：新聚类与原始 200k Leiden 聚类的一致性；
- 新 Leiden cluster 数、输入细胞数和 VMR 数。

最终结合指标和图片检查：

- 主要细胞类型结构是否保留；
- sample 聚集是否减弱；
- 聚类是否不再主要由个体或批次驱动。

不能只追求 sample 混合度：如果细胞类型纯度明显下降，说明 VMR 数量可能过少。top 5%、2%、1% 应在“细胞类型结构保留”和“样本效应减弱”之间综合选择。

### 实际聚类结果

作业 `162804` 成功完成，退出状态为 0。三套分析均使用 15,500 个 Scanpy 匹配细胞：

| 集合 | VMR | Leiden clusters | 细胞类型纯度 | sample纯度 | sample混合熵 | response混合熵 | ARI vs 原始200k |
|---|---:|---:|---:|---:|---:|---:|---:|
| top 5% | 3,015 | 28 | 0.6268 | 0.2875 | 0.8302 | 0.9618 | 0.2918 |
| top 2% | 1,206 | 39 | 0.5754 | 0.2599 | 0.8351 | 0.9884 | 0.2280 |
| top 1% | 603 | 144 | 0.5383 | 0.2446 | 0.8506 | 0.9860 | 0.1593 |

结果解释：

- top 5% 的细胞类型纯度最高、cluster 数最少、与原始 200k Leiden 的 ARI 最高，说明主要生物学结构保留最好；
- top 2% 的 sample 混合略优于 top 5%，但细胞类型纯度和 ARI 下降，cluster 数增加；
- top 1% 的 sample 混合熵最高，但产生 144 个 cluster，细胞类型纯度及 ARI 最低，提示 VMR 数量过少后出现明显聚类碎片化；
- 在这三个比例中，top 5% 是当前最稳妥的主分析选择，top 2% 和 top 1% 可作为比较结果；
- 这里只能比较三个新方案之间的 sample 混合程度。要严格判断相对于未过滤的原始 200k 分析是否改善，还需要用同一公式计算原始 200k 的 sample purity 和 sample mixing entropy 基线。

本次作业资源使用：

```text
Job ID：162804
运行时间：520秒
峰值内存：约5.2 GB
退出状态：0
```

注意：该流程属于个体差异特征过滤，不等同于在统计模型中加入 batch 协变量。
