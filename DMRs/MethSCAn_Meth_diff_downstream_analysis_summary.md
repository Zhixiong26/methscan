# 去 batch 后 MethSCAn / Meth_diff 下游分析流程整理

## 总体分析目标

本轮分析主线是：

```text
1. 基于 Meth_diff 结果，筛选 q value < 0.05 的 DMR
2. 确认 same-cell-type IR vs NR DMR 的方向
3. 提取每个 cell type 内 IR vs NR 的疾病相关 DMR
4. 注释 promoter 区域 DMR
5. 将 promoter DMR 映射到 gene
6. 整合 RNA pseudobulk expression
7. 筛选 promoter methylation 与 gene expression 变化方向相反的候选基因
```

核心逻辑：

```text
IR_hyper promoter + IR expression down
IR_hypo promoter  + IR expression up
```

也就是筛选 promoter 甲基化和表达变化方向相反的候选调控基因。

---

# 一、主要输入数据

## 1. Meth_diff DMR 输入目录

### 每个 cell type 内 IR vs NR DMR

这是疾病组 vs 对照组 DMR 的主输入：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR
```

包括：

```text
B_cells_IR_vs_NR_DMRs.bed
CD4_T_cells_IR_vs_NR_DMRs.bed
CD8_T_cells_IR_vs_NR_DMRs.bed
Monocytes_CD14_IR_vs_NR_DMRs.bed
Monocytes_CD16_IR_vs_NR_DMRs.bed
NK_cells_IR_vs_NR_DMRs.bed
pDCs_IR_vs_NR_DMRs.bed
Plasma_cells_IR_vs_NR_DMRs.bed
```

DMR BED 文件为 12 列，其中关键列为：

```text
第 1-3 列：chr, start, end
第 8 列：group_A mean methylation
第 9 列：group_B mean methylation
第 10 列：hypomethylated group
第 11 列：p value
第 12 列：q value
```

后续筛选使用：

```bash
$12 < 0.05
```

---

## 2. cell group 文件

用于确认 group_A / group_B 对应 IR 还是 NR：

```bash
/share/home/rzli/METHSCAN/Meth_diff/cell_groups_200k/3_same_cell_type_IR_vs_NR
```

示例：

```bash
/share/home/rzli/METHSCAN/Meth_diff/cell_groups_200k/3_same_cell_type_IR_vs_NR/CD4_T_cells_IR_vs_NR_cell_groups.csv
```

已确认：

```text
group_A = IR
group_B = NR
```

因此 Meth_diff 第 10 列方向解释为：

```text
group_A = IR hypomethylated = IR_hypo / NR_hyper
group_B = NR hypomethylated = IR_hyper / NR_hypo
```

也就是说：

```text
第10列 group_B → IR_hyper DMR
第10列 group_A → IR_hypo DMR
```

---

## 3. promoter 注释输入

使用 GENCODE v44 basic annotation：

```bash
/share/LCZX_Data/ref/gencode.v44.basic.annotation.gtf
```

从该 GTF 生成 promoter BED：

```text
promoter = TSS upstream 2 kb 到 downstream 500 bp
```

生成后的 promoter 文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/gencode_v44_basic_promoter_2kb_up_500bp_down.bed
```

promoter 数量：

```text
62663 promoters
```

---

## 4. RNA 表达数据输入

使用 Scanpy 整合后的 h5ad 文件：

```bash
/share/home/rzli/SCANPY/result/ALL_batch_corrected_pbmc.h5ad
```

检查结果：

```text
AnnData: 58534 cells × 2000 HVG genes
obs:
  sample
  batch
  group
  leiden_integrated
  cell_type_integrated

layers:
  counts
  log1p_uncorrected

raw:
  58534 cells × 38606 genes
```

表达整合时使用：

```text
adata.raw.X
```

原因：

```text
adata.X 只有 2000 个 highly variable genes；
adata.raw.X 有 38606 个 gene，更适合查 promoter DMR candidate genes 的表达。
```

---

# 二、主要脚本和作用

## 1. 统计每个 cell type 内 IR vs NR q<0.05 DMR 数量

脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/summarize_same_cell_type_IR_vs_NR_q005_direction_fixed.sh
```

输出：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/same_cell_type_IR_vs_NR_DMR_q005_summary_direction_fixed.tsv
```

核心结果：

```text
cell_type       total_DMR  q005_DMR  IR_hyper_DMR  IR_hypo_DMR
B_cells         94548      0         0             0
CD4_T_cells     90585      53        18            35
CD8_T_cells     35183      4         3             1
Monocytes_CD14  86070      3540      1457          2083
Monocytes_CD16  91041      3         2             1
NK_cells        89008      921       101           820
pDCs            5088       0         0             0
Plasma_cells    79728      1         0             1
```

解释：

```text
疾病相关 DMR 主要集中在 Monocytes_CD14 和 NK_cells。
Monocytes_CD14 的 DMR 数量最多，q<0.05 后有 3540 个。
NK_cells 次之，有 921 个。
```

---

## 2. 生成方向明确的 IR_hyper / IR_hypo DMR 文件

脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/filter_same_cell_type_IR_vs_NR_DMR_q005_direction_fixed.sh
```

输出目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed
```

主要输出文件格式：

```text
{cell_type}_IR_vs_NR_DMRs_q005.full.bed
{cell_type}_IR_vs_NR_DMRs_q005.bed
{cell_type}_IR_hyper_DMRs_q005.full.bed
{cell_type}_IR_hyper_DMRs_q005.bed
{cell_type}_IR_hypo_DMRs_q005.full.bed
{cell_type}_IR_hypo_DMRs_q005.bed
```

重点文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/Monocytes_CD14_IR_hyper_DMRs_q005.bed
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/Monocytes_CD14_IR_hypo_DMRs_q005.bed

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/NK_cells_IR_hyper_DMRs_q005.bed
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/NK_cells_IR_hypo_DMRs_q005.bed

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/CD4_T_cells_IR_hyper_DMRs_q005.bed
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed/CD4_T_cells_IR_hypo_DMRs_q005.bed
```

---

## 3. 生成 GENCODE v44 promoter BED

脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/make_gencode_v44_promoter_bed.py
```

输入：

```bash
/share/LCZX_Data/ref/gencode.v44.basic.annotation.gtf
```

输出：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/gencode_v44_basic_promoter_2kb_up_500bp_down.bed
```

输出格式：

```text
chr
promoter_start
promoter_end
gene_id
gene_name
gene_type
strand
tss
```

结果：

```text
62663 promoter regions
```

---

## 4. DMR 与 promoter overlap

脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/intersect_disease_DMR_with_promoter_q005.sh
```

输入：

```bash
DMR_DIR=/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed

PROMOTER=/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/gencode_v44_basic_promoter_2kb_up_500bp_down.bed
```

输出目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs
```

summary 文件：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_q005_summary.tsv
```

核心结果：

```text
cell_type       direction  all_DMR  promoter_DMR_overlap
CD4_T_cells     IR_hyper   18       10
CD4_T_cells     IR_hypo    35       19
Monocytes_CD14  IR_hyper   1457     240
Monocytes_CD14  IR_hypo    2083     424
NK_cells        IR_hyper   101      22
NK_cells        IR_hypo    820      159
```

解释：

```text
promoter DMR 仍主要集中在 Monocytes_CD14 和 NK_cells。
Monocytes_CD14 中 IR_hypo promoter DMR overlap 最多，有 424 个。
```

---

## 5. 整理 promoter DMR-to-gene 表

脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/format_promoter_DMR_gene_table_q005.py
```

输入目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs
```

输出主表：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_to_gene_q005.tsv
```

该表包含：

```text
cell_type
direction
gene_id
gene_name
gene_type
strand
tss
dmr_chr
dmr_start
dmr_end
promoter_chr
promoter_start
promoter_end
dmr_score
dmr_col5
dmr_col6
dmr_col7
dmr_col8
dmr_col9
hypomethylated_group
pvalue
qvalue
```

其中：

```text
dmr_col8 = IR mean methylation
dmr_col9 = NR mean methylation
```

结果：

```text
Total promoter DMR-gene overlaps: 874
```

gene_type 分布：

```text
lncRNA                                273
protein_coding                        232
processed_pseudogene                  144
unprocessed_pseudogene                 61
miRNA                                  39
misc_RNA                               20
TEC                                    17
```

protein-coding promoter DMR-gene overlap：

```text
direction       IR_hyper  IR_hypo
cell_type
CD4_T_cells            2        4
Monocytes_CD14        51      137
NK_cells               5       33
```

去重后 protein-coding gene：

```text
direction       IR_hyper  IR_hypo
cell_type
CD4_T_cells            2        4
Monocytes_CD14        51      133
NK_cells               5       33
```

---

## 6. 生成 protein-coding promoter DMR candidate gene 表

该步骤通过临时 Python 命令完成。

输入：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_to_gene_q005.tsv
```

输出：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_to_protein_coding_gene_q005.tsv
```

去重 candidate gene 表：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/protein_coding_promoter_DMR_candidate_genes_q005.tsv
```

这个表是后续 RNA 表达整合的输入。

---

## 7. RNA pseudobulk expression 与 promoter DMR 整合

主脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/pseudobulk_expression_for_promoter_DMR_genes.py
```

wrapper 脚本：

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/run_pseudobulk_expr_promoter_DMR.sh
```

提交命令：

```bash
dsub -n pseudobulk_expr_promoter_DMR \
  -R "cpu=64;mem=180G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/pseudobulk_expr_promoter_DMR.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/pseudobulk_expr_promoter_DMR.%J.err \
  bash /share/home/rzli/METHSCAN/Meth_diff/scripts/run_pseudobulk_expr_promoter_DMR.sh
```

JobID：

```text
162629
```

输入：

```bash
RNA h5ad:
/share/home/rzli/SCANPY/result/ALL_batch_corrected_pbmc.h5ad

candidate genes:
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/protein_coding_promoter_DMR_candidate_genes_q005.tsv
```

使用表达矩阵：

```text
adata.raw.X
```

分组单位：

```text
sample + response + cell_type
```

细胞类型映射：

```text
CD14_Monocytes → Monocytes_CD14
CD16_Monocytes → Monocytes_CD16
CD4_T_cells    → CD4_T_cells
CD8_T_cells    → CD8_T_cells
NK_cells       → NK_cells
B_cells        → B_cells
Plasma_cells   → Plasma_cells
pDCs           → pDCs
```

---

# 三、RNA expression integration 输出结果

输出目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration
```

文件列表：

```text
candidate_genes_missing_in_RNA_raw.tsv
negative_direction_promoter_methylation_expression_candidate_genes.tsv
promoter_DMR_candidate_genes_with_RNA_expression_delta.tsv
RNA_expression_delta_IR_vs_NR_candidate_genes.tsv
RNA_pseudobulk_mean_expression_candidate_genes_by_sample_response_celltype.tsv
RNA_pseudobulk_metadata.tsv
```

---

## 1. RNA pseudobulk 表达矩阵

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/RNA_pseudobulk_mean_expression_candidate_genes_by_sample_response_celltype.tsv
```

含义：

```text
行 = sample + response + cell_type pseudobulk
列 = candidate gene
值 = raw expression 的 pseudobulk mean
```

metadata：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/RNA_pseudobulk_metadata.tsv
```

---

## 2. RNA IR vs NR expression delta 表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/RNA_expression_delta_IR_vs_NR_candidate_genes.tsv
```

包含：

```text
cell_type
gene_name
IR_mean_expression
NR_mean_expression
delta_expression_IR_minus_NR
n_IR_pseudobulk
n_NR_pseudobulk
```

---

## 3. promoter methylation + RNA expression 整合表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/promoter_DMR_candidate_genes_with_RNA_expression_delta.tsv
```

结果：

```text
All merged: 228 rows
有 expression delta 的 gene: 224
缺失 expression delta 的 gene: 4
```

重要列：

```text
cell_type
direction
gene_name
delta_methylation_IR_minus_NR
delta_expression_IR_minus_NR
negative_direction
expected_negative_pattern
qvalue
```

---

## 4. 负相关方向候选基因表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/negative_direction_promoter_methylation_expression_candidate_genes.tsv
```

这是当前最重要的结果表。

筛选逻辑：

```text
delta_methylation_IR_minus_NR > 0 且 delta_expression_IR_minus_NR < 0
或者
delta_methylation_IR_minus_NR < 0 且 delta_expression_IR_minus_NR > 0
```

也就是：

```text
IR_hyper promoter + IR expression down
IR_hypo promoter  + IR expression up
```

结果统计：

```text
Negative direction by cell_type and direction:

direction       IR_hyper  IR_hypo
cell_type
CD4_T_cells            1        2
Monocytes_CD14        33       39
NK_cells               4       12
```

总数：

```text
CD4_T_cells:      3
Monocytes_CD14:  72
NK_cells:        16
Total:           91
```

---

# 四、当前结果的主要生物学解读

## 1. 疾病相关 DMR 主要集中在 CD14 Monocytes 和 NK cells

按 q<0.05 筛选后：

```text
Monocytes_CD14: 3540 DMRs
NK_cells:        921 DMRs
CD4_T_cells:      53 DMRs
```

说明 IR vs NR 的甲基化差异不是均匀分布在所有细胞类型，而是主要集中在：

```text
Monocytes_CD14
NK_cells
```

---

## 2. promoter DMR 也主要集中在 CD14 Monocytes 和 NK cells

promoter overlap 后：

```text
Monocytes_CD14:
  IR_hyper promoter DMR overlap = 240
  IR_hypo  promoter DMR overlap = 424

NK_cells:
  IR_hyper promoter DMR overlap = 22
  IR_hypo  promoter DMR overlap = 159

CD4_T_cells:
  IR_hyper promoter DMR overlap = 10
  IR_hypo  promoter DMR overlap = 19
```

说明疾病相关 promoter 甲基化差异也主要集中在 Monocytes_CD14 和 NK_cells。

---

## 3. protein-coding promoter DMR gene 主要集中在 Monocytes_CD14

protein-coding promoter gene 统计：

```text
Monocytes_CD14:
  IR_hyper = 51 genes
  IR_hypo  = 133 genes

NK_cells:
  IR_hyper = 5 genes
  IR_hypo  = 33 genes

CD4_T_cells:
  IR_hyper = 2 genes
  IR_hypo  = 4 genes
```

说明后续功能解释和候选基因筛选应优先放在：

```text
Monocytes_CD14
NK_cells
```

---

## 4. promoter methylation-expression 负相关方向候选基因

整合 RNA pseudobulk expression 后，得到：

```text
91 个 promoter methylation-expression inverse-direction candidate genes
```

其中：

```text
Monocytes_CD14: 72
NK_cells:       16
CD4_T_cells:     3
```

这说明 IR 相关的 promoter methylation-expression coupling 主要发生在：

```text
Monocytes_CD14
NK_cells
```

---

# 五、目前较值得关注的候选基因

## Monocytes_CD14: IR_hyper promoter + IR expression down

代表基因：

```text
ATG16L1
GTF3C6
PRCP
TIPRL
PLSCR2
MRPS22
RRP1B
TAF2
ZNF133
ITPKC
ABCA1
INTS6
TTC21B
```

其中相对更值得优先看：

```text
ATG16L1
ABCA1
INTS6
TAF2
PRCP
ITPKC
ZNF133
```

解释示例：

```text
ATG16L1 在 Monocytes_CD14 中表现为 IR promoter hypermethylation，同时 IR expression down，符合 promoter 甲基化升高抑制表达的方向。
```

---

## Monocytes_CD14: IR_hypo promoter + IR expression up

代表基因：

```text
SPSB3
VSTM1
TMEM14C
```

这些基因表现为：

```text
IR promoter hypomethylation
+
IR expression up
```

---

## CD4_T_cells

代表基因：

```text
PTPN7
CHI3L2
```

其中：

```text
PTPN7: IR_hyper promoter + IR expression down
CHI3L2: IR_hypo promoter + IR expression up
```

---

# 六、当前结果的注意事项

## 1. 现在是“负相关方向候选基因”，不是严格显著负相关基因

当前做法是：

```text
IR mean methylation - NR mean methylation
IR mean expression - NR mean expression
```

然后筛选方向相反的基因。

因此现在结果应称为：

```text
promoter methylation-expression inverse-direction candidate genes
```

中文：

```text
promoter 甲基化与表达变化方向相反的候选基因
```

或者：

```text
负相关方向候选基因
```

暂时不要称为：

```text
显著负相关基因
```

因为还没有做 sample-level Pearson/Spearman correlation 和 p value。

---

## 2. 需要注意表达变化幅度

有些 gene 虽然方向相反，但：

```text
delta_expression_IR_minus_NR 非常接近 0
```

例如：

```text
0.000038
0.000063
0.000138
```

这种基因生物意义可能较弱。

建议后续加一个温和过滤：

```text
abs(delta_methylation_IR_minus_NR) >= 0.03
abs(delta_expression_IR_minus_NR) >= 0.005
```

---

# 七、建议下一步分析

## 1. 生成 clean negative-direction candidate table

建议输出两个表：

```text
1. all negative-direction genes
2. filtered negative-direction genes
```

过滤条件：

```text
abs(delta_methylation) >= 0.03
abs(delta_expression) >= 0.005
```

推荐输出路径：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/negative_direction_candidate_genes_clean_table.tsv

/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/negative_direction_candidate_genes_filtered_absMeth003_absExpr0005.tsv
```

---

## 2. 做 sample-level consistency / correlation

下一步更严格的分析应该是：

```text
每个 sample-cell_type pseudobulk 中：
promoter methylation value
RNA expression value
```

然后对每个 candidate gene 做：

```text
Spearman correlation
Pearson correlation
correlation p value
FDR
```

目的：

```text
排除只由单个 sample 驱动的 gene
筛选真正稳健的 promoter methylation-expression negative correlation genes
```

---

## 3. 做可视化

建议最终图包括：

```text
1. 每个 cell type 的 IR vs NR DMR 数量柱状图
2. 每个 cell type 的 promoter DMR 数量柱状图
3. protein-coding promoter DMR gene 数量图
4. negative-direction candidate gene 数量图
5. Monocytes_CD14 top candidate gene scatter plot
   x = promoter methylation
   y = RNA expression
   color = IR / NR
6. top candidate gene heatmap
```

优先展示：

```text
Monocytes_CD14
NK_cells
```

候选基因可优先展示：

```text
ATG16L1
ABCA1
INTS6
TAF2
PRCP
ITPKC
VSTM1
TMEM14C
PTPN7
CHI3L2
```

---

# 八、目前最重要路径汇总

## DMR 方向修正 summary

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/same_cell_type_IR_vs_NR_DMR_q005_summary_direction_fixed.tsv
```

## 方向明确的 DMR 文件目录

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed
```

## promoter BED

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/gencode_v44_basic_promoter_2kb_up_500bp_down.bed
```

## promoter DMR summary

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_q005_summary.tsv
```

## promoter DMR-to-gene 总表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_to_gene_q005.tsv
```

## protein-coding promoter DMR gene 表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/promoter_DMR_to_protein_coding_gene_q005.tsv
```

## 去重 protein-coding candidate gene 表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/promoter_annotation/per_cell_type_promoter_DMRs/protein_coding_promoter_DMR_candidate_genes_q005.tsv
```

## RNA h5ad

```bash
/share/home/rzli/SCANPY/result/ALL_batch_corrected_pbmc.h5ad
```

## RNA pseudobulk expression 输出目录

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration
```

## promoter methylation + RNA expression 整合表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/promoter_DMR_candidate_genes_with_RNA_expression_delta.tsv
```

## 负相关方向候选基因表

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/expression_integration/negative_direction_promoter_methylation_expression_candidate_genes.tsv
```

---

# 九、可以写进汇报的简洁版本

```text
We performed cell-type-specific IR vs NR DMR analysis using Meth_diff results and filtered DMRs by q value < 0.05. The group direction was confirmed using the cell group files, where group_A corresponds to IR and group_B corresponds to NR. Since the 10th column of the Meth_diff BED indicates the hypomethylated group, group_B DMRs were interpreted as IR-hypermethylated DMRs, while group_A DMRs were interpreted as IR-hypomethylated DMRs.

Disease-associated DMRs were mainly detected in CD14 monocytes and NK cells. After intersecting DMRs with promoter regions defined as TSS -2 kb to +500 bp using GENCODE v44 annotation, CD14 monocytes showed the largest number of promoter DMRs. Protein-coding promoter DMR genes were then integrated with RNA pseudobulk expression from the Scanpy h5ad object using adata.raw.X.

We identified 91 promoter methylation-expression inverse-direction candidate genes, including 72 in CD14 monocytes, 16 in NK cells, and 3 in CD4 T cells. These candidates include genes showing either IR promoter hypermethylation with reduced IR expression or IR promoter hypomethylation with increased IR expression. The results suggest that IR-associated promoter methylation-expression coupling is mainly concentrated in CD14 monocytes and NK cells.
```
