# q<0.05：肺癌免疫治疗应答相关 DMR、promoter 与 DNA-RNA 整合结果解读

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

promoter 定义为 GENCODE v44 注释的 TSS 上游 2 kb 至下游 2 kb。

## 2. IR vs NR DMR 结果

| cell type | q<0.05 DMR | IR-hyper | IR-hypo |
|---|---:|---:|---:|
| B cells | 0 | 0 | 0 |
| CD4 T cells | 53 | 18 | 35 |
| CD8 T cells | 4 | 3 | 1 |
| CD14 Monocytes | 3,540 | 1,457 | 2,083 |
| CD16 Monocytes | 3 | 2 | 1 |
| NK cells | 921 | 101 | 820 |
| pDCs | 0 | 0 | 0 |
| Plasma cells | 1 | 0 | 1 |

主要结论：IR/NR 甲基化差异高度集中在 **CD14 Monocytes**，其次为 **NK cells**；CD4 T cells 仅有少量差异，其他细胞类型信号有限。

## 3. Clean cell-type DMR / VMR 分支

该分支用于构建较少受 IR/NR 应答状态影响的细胞类型甲基化特征，不用于疾病应答 DMR 发现。

```text
cell-type pairwise DMRs q<0.05：149,867
IR-vs-NR sample/response component：4,353
clean cell-type DMRs：147,974
映射到 All VMR matrix 的 clean VMR features：29,863
```

## 4. Promoter DMR 结果

| cell type | direction | all DMR | unique DMR with promoter overlap |
|---|---|---:|---:|
| CD4 T cells | IR-hyper | 18 | 9 |
| CD4 T cells | IR-hypo | 35 | 9 |
| CD14 Monocytes | IR-hyper | 1,457 | 237 |
| CD14 Monocytes | IR-hypo | 2,083 | 370 |
| NK cells | IR-hyper | 101 | 23 |
| NK cells | IR-hypo | 820 | 140 |

promoter-level 结果与全局 DMR 一致：主要信号集中在 CD14 Monocytes 和 NK cells，且两者都以 IR-hypo DMR 较多。

## 5. DMR-to-gene 与 RNA pseudobulk 整合

```text
promoter DMR-gene overlap rows：1,065
protein-coding DMR-gene rows：294
去重后 protein-coding candidate genes：289
RNA raw matrix 中存在的候选基因：271
缺失于 RNA raw matrix 的候选基因：8
具有 RNA expression delta 的 DMR-gene rows：281
promoter methylation 与 RNA expression 反向方向候选：115
```

反向方向定义：

```text
IR-hyper promoter + IR expression down
IR-hypo promoter  + IR expression up
```

115 个反向方向候选的分布：

| cell type | IR-hyper + expression down | IR-hypo + expression up |
|---|---:|---:|
| CD4 T cells | 1 | 2 |
| CD14 Monocytes | 47 | 47 |
| NK cells | 6 | 12 |

因此，CD14 Monocytes 是最主要的免疫治疗应答相关 promoter methylation-expression 候选细胞类型。

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
candidate DMR-gene rows：289
不与任何 VMR feature 重叠的候选：137
实际可进行 DNA-RNA 配对的候选：约 152
DNA-RNA paired records：1,492
```

每个 cell type 内仅有 5 个 IR 与 5 个 NR 样本。当前设置 `min_paired_samples=6`，因此：

```text
pooled IR+NR：n=10，可计算相关性
IR 内部：n=5，不计算
NR 内部：n=5，不计算
```

`pooled_descriptive` 结果可用于描述性筛选，但可能受 IR/NR 组间差异驱动，不能替代组内相关性或因果证据。

### 多重检验结果

在所有 pooled correlation 中：

```text
BH-FDR < 0.05：0 个候选
```

因此，没有 DNA methylation–RNA correlation 在多重检验校正后显著。

### 原始 p 值层面的探索性候选

要求 Pearson 与 Spearman 均为负相关，且两种原始 p value 均 <0.05，共得到 2 个候选：

| cell type | gene | DMR direction | Pearson r / p | Spearman rho / p |
|---|---|---|---|---|
| CD4 T cells | PTPN7 | IR-hyper | -0.711 / 0.021 | -0.745 / 0.013 |
| CD14 Monocytes | ZFP57 | IR-hyper | -0.796 / 0.0059 | -0.711 / 0.021 |

两者均体现：应答者相关的 promoter 高甲基化，且在 pooled 样本中甲基化升高与表达降低趋势一致。ZFP57 的原始相关强度更高；PTPN7 与 T 细胞 TCR/MAPK 信号调节相关，ZFP57 与 DNA methylation imprinting/表观遗传调控相关。

但二者的相关性 FDR 均未达到显著，因此只能定义为：

> Exploratory nominal inverse methylation-expression candidates.

不能定义为 FDR-significant correlation genes。

## 7. 推荐汇报结论

> 在肺癌免疫治疗应答者与不应答者的同细胞类型比较中，应答相关 DNA methylation 差异主要集中在 CD14 Monocytes，其次为 NK cells。Promoter DMR 与 RNA pseudobulk 整合后，CD14 Monocytes 仍表现出最强的候选 methylation-expression coupling。匹配样本的 DNA-RNA correlation 分析发现 PTPN7 和 ZFP57 具有原始 p 值层面的负相关趋势，但在 10 个样本和多重检验校正后没有相关性达到 FDR 显著。因此，应将二者作为优先验证的探索性候选，而非确定性调控基因。

## 8. 输出文件

```text
DMR summary:
DMR_disease_200k_q005/same_cell_type_IR_vs_NR_DMR_q005_summary_direction_fixed.tsv

Promoter summary:
DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_q005_summary_2kb_2kb.tsv

Protein-coding candidates:
DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/protein_coding_promoter_DMR_candidate_genes_q005_2kb_2kb.tsv

RNA integration:
DMR_disease_200k_q005/expression_integration_2kb_2kb/negative_direction_promoter_methylation_expression_candidate_genes.tsv

DNA-RNA correlation:
DMR_disease_200k_q005/expression_integration_2kb_2kb/DNA_RNA_correlation/DNA_methylation_RNA_expression_correlation.tsv
```
