# Meth_diff IR vs NR：q<0.01 与 q<0.05 结果总结

## 1. 分析目标

本次分析基于 Meth_diff 的同一细胞类型内 IR vs NR 差异甲基化结果，分别使用两个显著性阈值：

- **q < 0.05**：相对宽松，用于发现更多候选 DMR / 基因。
- **q < 0.01**：更严格，用于筛选更高置信度候选 DMR / 基因。

后续分析统一使用：

- **promoter 定义：TSS ±2 kb**
- **DMR 方向定义：**
  - `IR_hyper`：IR 组甲基化水平高于 NR 组。
  - `IR_hypo`：IR 组甲基化水平低于 NR 组。
- **候选基因筛选逻辑：**
  1. 同一 cell type 内做 IR vs NR DMR。
  2. 筛选 promoter 区域 DMR。
  3. 注释到 gene。
  4. 保留 protein-coding genes。
  5. 与 RNA pseudobulk expression 做整合。
  6. 保留 promoter methylation 与 RNA expression 呈负相关方向的候选基因。

---

## 2. q<0.05 主结果

### 2.1 q<0.05 DMR summary

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/same_cell_type_IR_vs_NR_DMR_q005_summary_direction_fixed.tsv
```

| cell_type | total_DMR | q005_DMR | IR_hyper_DMR | IR_hypo_DMR |
|---|---:|---:|---:|---:|
| B_cells | 94548 | 0 | 0 | 0 |
| CD4_T_cells | 90585 | 53 | 18 | 35 |
| CD8_T_cells | 35183 | 4 | 3 | 1 |
| Monocytes_CD14 | 86070 | 3540 | 1457 | 2083 |
| Monocytes_CD16 | 91041 | 3 | 2 | 1 |
| NK_cells | 89008 | 921 | 101 | 820 |
| pDCs | 5088 | 0 | 0 | 0 |
| Plasma_cells | 79728 | 1 | 0 | 1 |

### 2.2 q<0.05 promoter DMR overlap，TSS ±2 kb

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_q005_summary_2kb_2kb.tsv
```

| cell_type | direction | all_DMR | promoter_DMR_overlap |
|---|---|---:|---:|
| CD4_T_cells | IR_hyper | 18 | 9 |
| CD4_T_cells | IR_hypo | 35 | 9 |
| CD8_T_cells | IR_hyper | 3 | 0 |
| CD8_T_cells | IR_hypo | 1 | 0 |
| Monocytes_CD14 | IR_hyper | 1457 | 237 |
| Monocytes_CD14 | IR_hypo | 2083 | 370 |
| Monocytes_CD16 | IR_hyper | 2 | 0 |
| Monocytes_CD16 | IR_hypo | 1 | 0 |
| NK_cells | IR_hyper | 101 | 23 |
| NK_cells | IR_hypo | 820 | 140 |
| Plasma_cells | IR_hypo | 1 | 0 |

主要结论：

- q<0.05 下，promoter DMR 主要集中在 **Monocytes_CD14** 和 **NK_cells**。
- Monocytes_CD14 的 promoter DMR 数量最多：
  - IR_hyper promoter DMR：237
  - IR_hypo promoter DMR：370
- NK_cells 也有较明显信号：
  - IR_hyper promoter DMR：23
  - IR_hypo promoter DMR：140

### 2.3 q<0.05 promoter DMR-to-gene / protein-coding gene

主要文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_to_gene_q005_2kb_2kb.tsv

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_to_protein_coding_gene_q005_2kb_2kb.tsv

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/protein_coding_promoter_DMR_candidate_genes_q005_2kb_2kb.tsv
```

### 2.4 q<0.05 RNA integration 结果

RNA integration 结果目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration_2kb_2kb
```

主要结果：

| category | count |
|---|---:|
| All merged promoter protein-coding candidates | 289 |
| Genes with RNA expression delta | 281 |
| Missing expression or methylation | 8 |
| Negative-direction candidates | 115 |

负相关方向定义：

- `IR_hyper_promoter_and_IR_expression_down`
- `IR_hypo_promoter_and_IR_expression_up`

q<0.05 负相关方向分布：

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| CD4_T_cells | 1 | 2 |
| Monocytes_CD14 | 47 | 47 |
| NK_cells | 6 | 12 |

### 2.5 q<0.05 最终过滤结果

过滤条件：

```text
abs(delta_methylation_IR_minus_NR) >= 0.03
abs(delta_expression_IR_minus_NR) >= 0.005
```

最终过滤后：

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| CD4_T_cells | 1 | 1 |
| Monocytes_CD14 | 20 | 11 |
| NK_cells | 5 | 2 |

总数：**40 个候选基因**。

核心文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration_2kb_2kb/negative_direction_candidate_genes_clean_table_2kb_2kb.tsv

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_absMeth003_absExpr0005_2kb_2kb.tsv
```

---

## 3. q<0.01 主结果

### 3.1 q<0.01 DMR summary

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/same_cell_type_IR_vs_NR_DMR_q001_summary_direction_fixed.tsv
```

| cell_type | total_DMR | q001_DMR | IR_hyper_DMR | IR_hypo_DMR |
|---|---:|---:|---:|---:|
| B_cells | 94548 | 0 | 0 | 0 |
| CD4_T_cells | 90585 | 2 | 2 | 0 |
| CD8_T_cells | 35183 | 4 | 3 | 1 |
| Monocytes_CD14 | 86070 | 1016 | 427 | 589 |
| Monocytes_CD16 | 91041 | 3 | 2 | 1 |
| NK_cells | 89008 | 8 | 2 | 6 |
| pDCs | 5088 | 0 | 0 | 0 |
| Plasma_cells | 79728 | 1 | 0 | 1 |

主要结论：

- q<0.01 后，DMR 数量明显减少。
- 主要信号集中在 **Monocytes_CD14**。
- NK_cells 从 q<0.05 的 921 个 DMR 降到 q<0.01 的 8 个 DMR，说明 NK_cells 中较多信号在严格阈值下不稳定。

### 3.2 q<0.01 promoter DMR overlap，TSS ±2 kb

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_q001_summary_2kb_2kb.tsv
```

| cell_type | direction | all_DMR | promoter_DMR_overlap |
|---|---|---:|---:|
| CD4_T_cells | IR_hyper | 2 | 2 |
| CD4_T_cells | IR_hypo | 0 | 0 |
| CD8_T_cells | IR_hyper | 3 | 0 |
| CD8_T_cells | IR_hypo | 1 | 0 |
| Monocytes_CD14 | IR_hyper | 427 | 83 |
| Monocytes_CD14 | IR_hypo | 589 | 118 |
| Monocytes_CD16 | IR_hyper | 2 | 0 |
| Monocytes_CD16 | IR_hypo | 1 | 0 |
| NK_cells | IR_hyper | 2 | 2 |
| NK_cells | IR_hypo | 6 | 4 |
| Plasma_cells | IR_hypo | 1 | 0 |

主要结论：

- q<0.01 下 promoter DMR 仍主要集中在 **Monocytes_CD14**。
- Monocytes_CD14：
  - IR_hyper promoter DMR：83
  - IR_hypo promoter DMR：118
- NK_cells 只剩少量 promoter DMR：
  - IR_hyper：2
  - IR_hypo：4

### 3.3 q<0.01 promoter DMR-to-gene / protein-coding gene

第 8 步结果：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_to_gene_q001_2kb_2kb.tsv
```

结果规模：

```text
All promoter DMR-gene rows: 313
```

gene type 主要分布：

| gene_type | count |
|---|---:|
| lncRNA | 90 |
| protein_coding | 85 |
| processed_pseudogene | 53 |
| IG_D_gene | 16 |
| unprocessed_pseudogene | 14 |

第 9 步 protein-coding 结果：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/promoter_DMR_to_protein_coding_gene_q001_2kb_2kb.tsv

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/protein_coding_promoter_DMR_candidate_genes_q001_2kb_2kb.tsv
```

结果规模：

| category | count |
|---|---:|
| All promoter DMR-gene rows | 313 |
| Protein-coding DMR-gene rows | 85 |
| Deduplicated candidate genes | 84 |

protein-coding DMR-gene overlaps：

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| CD4_T_cells | 1 | 0 |
| Monocytes_CD14 | 24 | 58 |
| NK_cells | 1 | 1 |

unique protein-coding candidate genes：

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| CD4_T_cells | 1 | 0 |
| Monocytes_CD14 | 24 | 57 |
| NK_cells | 1 | 1 |

### 3.4 q<0.01 RNA integration 结果

RNA integration 结果目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb
```

主要输出文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/RNA_pseudobulk_metadata.tsv
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/RNA_pseudobulk_mean_expression_candidate_genes_by_sample_response_celltype.tsv
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/RNA_expression_delta_IR_vs_NR_candidate_genes.tsv
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/promoter_DMR_candidate_genes_with_RNA_expression_delta.tsv
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_promoter_methylation_expression_candidate_genes.tsv
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/candidate_genes_missing_in_RNA_raw.tsv
```

RNA integration 结果规模：

| category | count |
|---|---:|
| Candidate table | 84 |
| Unique candidate gene symbols | 81 |
| Genes present in RNA raw | 80 |
| Genes missing in RNA raw | 1 |
| Cells used for RNA pseudobulk | 37187 |
| Pseudobulk groups | 30 |
| All merged | 84 |
| Expression delta not null | 83 |
| Missing expression or methylation | 1 |
| Negative-direction candidates | 36 |

RNA pseudobulk 使用细胞数：

| methdiff_cell_type | IR | NR |
|---|---:|---:|
| CD4_T_cells | 9373 | 5590 |
| Monocytes_CD14 | 8515 | 5233 |
| NK_cells | 5099 | 3377 |

expected negative pattern：

| expected_negative_pattern | count |
|---|---:|
| other | 47 |
| IR_hyper_promoter_and_IR_expression_down | 18 |
| IR_hypo_promoter_and_IR_expression_up | 18 |
| missing_expression_or_methylation | 1 |

negative-direction candidates：

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| Monocytes_CD14 | 17 | 17 |
| NK_cells | 1 | 1 |

### 3.5 q<0.01 最终过滤结果

所有 methylation-expression 反向候选基因：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_candidate_genes_clean_table_q001_2kb_2kb.tsv
```

总数：**36 个**。

过滤条件 1：

```text
abs(delta_methylation_IR_minus_NR) >= 0.03
abs(delta_expression_IR_minus_NR) >= 0.005
```

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_q001_absMeth003_absExpr0005_2kb_2kb.tsv
```

结果：**13 个候选基因**。

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| Monocytes_CD14 | 7 | 5 |
| NK_cells | 0 | 1 |

过滤条件 2：

```text
abs(delta_methylation_IR_minus_NR) >= 0.05
abs(delta_expression_IR_minus_NR) >= 0.005
```

结果文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_q001_absMeth005_absExpr0005_2kb_2kb.tsv
```

结果：**10 个候选基因**。

| cell_type | IR_hyper | IR_hypo |
|---|---:|---:|
| Monocytes_CD14 | 4 | 5 |
| NK_cells | 0 | 1 |

---

## 4. q<0.05 与 q<0.01 对比

### 4.1 DMR 数量对比

| cell_type | q<0.05 DMR | q<0.01 DMR | 变化 |
|---|---:|---:|---|
| CD4_T_cells | 53 | 2 | 大幅减少 |
| CD8_T_cells | 4 | 4 | 基本不变 |
| Monocytes_CD14 | 3540 | 1016 | 仍保留大量信号 |
| Monocytes_CD16 | 3 | 3 | 基本不变，但总数很少 |
| NK_cells | 921 | 8 | 大幅减少 |
| Plasma_cells | 1 | 1 | 总数很少 |

主要结论：

- **Monocytes_CD14 是最稳定、最主要的 IR vs NR 甲基化差异来源。**
- q<0.01 后，Monocytes_CD14 仍保留 1016 个 DMR，说明该细胞类型内 IR vs NR 差异较强。
- NK_cells 在 q<0.05 下有明显信号，但 q<0.01 后大幅减少，说明 NK_cells 信号对阈值比较敏感。

### 4.2 promoter DMR 对比

| cell_type | direction | q<0.05 promoter DMR | q<0.01 promoter DMR |
|---|---|---:|---:|
| CD4_T_cells | IR_hyper | 9 | 2 |
| CD4_T_cells | IR_hypo | 9 | 0 |
| Monocytes_CD14 | IR_hyper | 237 | 83 |
| Monocytes_CD14 | IR_hypo | 370 | 118 |
| NK_cells | IR_hyper | 23 | 2 |
| NK_cells | IR_hypo | 140 | 4 |

主要结论：

- promoter 层面与整体 DMR 层面一致，q<0.01 后主要保留 **Monocytes_CD14**。
- Monocytes_CD14 中 IR_hypo promoter DMR 数量略多于 IR_hyper promoter DMR。

### 4.3 RNA integration 负相关候选基因对比

| threshold | all negative-direction candidates | filtered absMeth>=0.03 & absExpr>=0.005 | filtered absMeth>=0.05 & absExpr>=0.005 |
|---|---:|---:|---:|
| q<0.05 | 115 | 40 | not generated / not used as main table |
| q<0.01 | 36 | 13 | 10 |

主要结论：

- q<0.05 得到更多候选，适合探索。
- q<0.01 得到更少但更高置信度的候选，适合作为最终重点结果。
- 两个阈值共同指向 **Monocytes_CD14**。

---

## 5. 推荐用于汇报的重点结论

1. 在同一 cell type 内比较 IR vs NR 后，DMR 信号主要集中在 **Monocytes_CD14**。
2. q<0.05 下可以观察到 Monocytes_CD14 和 NK_cells 都有较多 DMR，但 q<0.01 后 NK_cells 信号大幅减少。
3. q<0.01 下，Monocytes_CD14 仍保留大量 DMR 和 promoter DMR，说明 Monocytes_CD14 是更稳定的疾病相关甲基化差异细胞类型。
4. promoter DMR 与 RNA pseudobulk expression 整合后，q<0.01 得到 36 个 methylation-expression 反向候选基因。
5. 进一步按甲基化差异幅度和表达差异幅度过滤后：
   - `absMeth >= 0.03 & absExpr >= 0.005`：13 个候选基因。
   - `absMeth >= 0.05 & absExpr >= 0.005`：10 个高置信候选基因。
6. 最终主结果建议使用 q<0.01 的 13 个候选基因作为主要候选表，q<0.05 的 40 个候选基因作为补充探索结果。

---

## 6. 最终推荐主表

### q<0.05 推荐补充表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_absMeth003_absExpr0005_2kb_2kb.tsv
```

### q<0.01 推荐主表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_q001_absMeth003_absExpr0005_2kb_2kb.tsv
```

### q<0.01 高置信表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/expression_integration_2kb_2kb/negative_direction_candidate_genes_filtered_q001_absMeth005_absExpr0005_2kb_2kb.tsv
```

---

## 7. 一句话总结

q<0.05 用于探索时能发现更多候选基因，但 q<0.01 后仍然稳定保留的信号主要集中在 **Monocytes_CD14**；结合 promoter DMR 与 RNA expression 负相关方向后，q<0.01 最终筛得 13 个主要候选基因，其中 10 个在更严格甲基化差异阈值下仍保留，适合作为后续重点分析对象。
