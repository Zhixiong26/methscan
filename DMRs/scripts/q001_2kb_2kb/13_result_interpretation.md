# q<0.01：肺癌免疫治疗应答相关 DMR、promoter 与 DNA-RNA 整合结果解读

## 1. 分析对象

```text
IR = immunotherapy responder，肺癌免疫治疗应答者
NR = immunotherapy non-responder，肺癌免疫治疗不应答者
```

本分析在同一细胞类型内部比较 IR 与 NR。`group_A=IR`、`group_B=NR`；Meth_diff 第 10 列为低甲基化组，因此：

```text
group_B → NR 低甲基化 → IR_hyper
group_A → IR 低甲基化 → IR_hypo
```

promoter 定义为 GENCODE v44 注释的 TSS 上游 2 kb 至下游 2 kb。本目录使用严格阈值 `q<0.01`，用于获得高置信度候选。

## 2. IR vs NR DMR 结果

| cell type | q<0.01 DMR | IR-hyper | IR-hypo |
|---|---:|---:|---:|
| B cells | 0 | 0 | 0 |
| CD4 T cells | 2 | 2 | 0 |
| CD8 T cells | 4 | 3 | 1 |
| CD14 Monocytes | 1,016 | 427 | 589 |
| CD16 Monocytes | 3 | 2 | 1 |
| NK cells | 8 | 2 | 6 |
| pDCs | 0 | 0 | 0 |
| Plasma cells | 1 | 0 | 1 |

主要结论：在严格 q<0.01 阈值下，应答相关 DNA methylation 差异几乎完全由 **CD14 Monocytes** 驱动。NK cells 从 q<0.05 的 921 个 DMR 降至 8 个，表明其多数 q<0.05 信号未通过严格筛选。

## 3. Clean cell-type DMR / VMR 分支

该分支用于构建较少受 IR/NR 应答状态影响的细胞类型甲基化特征，不用于发现疾病/应答 DMR。

```text
clean cell-type DMR 映射到 All VMR matrix 的 VMR features：25,706
```

q<0.01 比 q<0.05 更严格，因此相应 clean VMR feature 数少于 q<0.05 的 29,863。

## 4. Promoter DMR 结果

| cell type | direction | all DMR | unique DMR with promoter overlap |
|---|---|---:|---:|
| CD4 T cells | IR-hyper | 2 | 2 |
| CD4 T cells | IR-hypo | 0 | 0 |
| CD14 Monocytes | IR-hyper | 427 | 83 |
| CD14 Monocytes | IR-hypo | 589 | 118 |
| NK cells | IR-hyper | 2 | 2 |
| NK cells | IR-hypo | 6 | 4 |

严格阈值下，promoter DMR 仍集中在 CD14 Monocytes：共 201 个唯一 DMR 与 promoter 重叠。NK cells 只保留 6 个 promoter DMR，CD4 T cells 保留 2 个。

## 5. DMR-to-gene 与 RNA pseudobulk 整合

```text
promoter DMR-gene overlap rows：313
protein-coding DMR-gene rows：85
去重后 protein-coding candidate genes：84
RNA raw matrix 中存在的候选基因：80
缺失于 RNA raw matrix 的候选基因：1
具有 RNA expression delta 的 DMR-gene rows：83
promoter methylation 与 RNA expression 反向方向候选：36
```

反向方向定义：

```text
IR-hyper promoter + IR expression down
IR-hypo promoter  + IR expression up
```

36 个反向方向候选的分布：

| cell type | IR-hyper + expression down | IR-hypo + expression up |
|---|---:|---:|
| CD14 Monocytes | 17 | 17 |
| NK cells | 1 | 1 |

在 q<0.01 下，CD4 T cells 的 promoter DMR 虽然存在，但没有同时满足 RNA 反向方向的候选；严格筛选后的 methylation-expression 候选几乎全部集中在 CD14 Monocytes。

## 6. Sample-level DNA methylation–RNA correlation

DNA pseudobulk 定义为同一 sample × cell type 内：

```text
sum(methylated sites) / sum(total sites)
```

DNA 与 RNA 使用相同的 10 个样本：

```text
IR01–IR05
NR01–NR05
```

相关性分析结果：

```text
candidate DMR-gene rows：84
不与任何 VMR feature 重叠的候选：40
实际可进行 DNA-RNA 配对的候选：约 43
selected VMR features：41
DNA-RNA paired records：430
DNA sample × cell type records excluded：0
```

每个 cell type 内仅有 5 个 IR 与 5 个 NR 样本。当前设置 `min_paired_samples=6`，因此：

```text
pooled IR+NR：n=10，可计算相关性
IR 内部：n=5，不计算
NR 内部：n=5，不计算
```

`pooled_descriptive` 结果可用于描述性筛选，但可能受 IR/NR 组间差异驱动，不能替代组内相关性或因果证据。

### 多重检验结果

```text
BH-FDR < 0.05：0 个候选
```

因此，严格 q<0.01 候选中没有 DNA methylation–RNA correlation 在多重检验校正后显著。

### 原始 p 值层面的探索性候选

要求 Pearson 与 Spearman 均为负相关，且两种原始 p value 均 <0.05，仅得到 1 个候选：

| cell type | gene | DMR direction | Pearson r / p | Spearman rho / p |
|---|---|---|---|---|
| CD14 Monocytes | ZFP57 | IR-hyper | -0.796 / 0.0059 | -0.711 / 0.021 |

ZFP57 在 q<0.05 和 q<0.01 结果中均保留，DMR 坐标均为：

```text
chr6:29678956-29682956
```

ZFP57 编码 KRAB 锌指蛋白，参与 DNA methylation imprinting/表观遗传调控；其在 CD14 Monocytes 中具有 IR-hyper promoter DMR，并显示一致的 DNA methylation–RNA expression 负相关趋势。因此，ZFP57 是当前最高优先级的探索性候选。

但其相关性 FDR 未显著：

```text
Pearson FDR = 0.213
Spearman FDR = 0.591
```

因此应定义为：

> High-priority, cross-threshold, nominal inverse methylation-expression candidate.

不能定义为 FDR-significant correlation gene。

## 7. 推荐汇报结论

> 在严格 q<0.01 阈值下，肺癌免疫治疗应答相关 DNA methylation 差异主要稳定地集中在 CD14 Monocytes。Promoter DMR 与 RNA pseudobulk 整合后，36 个反向方向候选中绝大多数仍来自 CD14 Monocytes。匹配样本的 DNA-RNA correlation 分析未发现经 FDR 校正后显著的相关性；但 ZFP57 在 CD14 Monocytes 中同时通过 q<0.05 和 q<0.01 的 DMR 筛选，并具有一致的 promoter 高甲基化、RNA 表达降低及原始 p 值层面负相关趋势，因此是最值得后续独立队列和功能实验验证的探索性候选。

## 8. 输出文件

```text
DMR summary:
DMR_disease_200k_q001/same_cell_type_IR_vs_NR_DMR_q001_summary_direction_fixed.tsv

Promoter summary:
DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_q001_summary_2kb_2kb.tsv

Protein-coding candidates:
DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/protein_coding_promoter_DMR_candidate_genes_q001_2kb_2kb.tsv

RNA integration:
DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_promoter_methylation_expression_candidate_genes.tsv

DNA-RNA correlation:
DMR_disease_200k_q001/expression_integration_2kb_2kb/DNA_RNA_correlation/DNA_methylation_RNA_expression_correlation.tsv
```
