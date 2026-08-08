# 甲基化概率与基因表达联合气泡图

本目录从现有10样本Top200 DMR matrix和Scanpy注释对象生成三个并排面板：

1. Promoter methylation
2. Gene-body methylation
3. Gene expression

默认输出一张合并10个样本的描述性图。中间结果保留`sample`、`cell`和`cell_type`，便于后续改成分样本或IR/NR分组图。

注意：10个样本使用各自筛选出的Top200 DMR集合，因此合并图适合描述“现有样本特异DMR所覆盖区域”的模式，不是使用统一全基因组区域做的无偏样本间比较。

## 固定输入

- RNA：`/share/home/rzli/SCANPY/20260714/result/annotation/02_annotated_final.h5ad`
- RNA表达：`adata.raw.X`（完整基因、log-normalized）
- RNA cell type：`cell_type_integrated`
- hg38：`/share/LCZX_Data/ref/gencode.v44.basic.annotation.gtf`（GENCODE v44）
- 甲基化：每个样本自己的300k Top200 DMR mean-ratio matrix
- cov：每个样本自己的`cov_dedup_probability`

## 脚本顺序

| 步骤 | 脚本 | 作用 |
|---:|---|---|
| 01 | `01_audit_joint_inputs.py` | 检查RNA、DMR、cov、cell type、染色体和基因ID |
| 02 | `02_select_top5_markers.py` | cell type one-vs-rest marker并全局去重，每类最多Top5 |
| 03 | `03_make_hg38_gene_regions.py` | 从GENCODE v44生成promoter和gene-body BED |
| 04 | `04_compute_dmr_unique_cpg_counts.py` | 补算与mean-ratio matrix逐元素对应的unique-CpG-count matrix |
| 05 | `05_map_dmrs_to_gene_regions.py` | DMR与marker promoter/gene body区间重叠 |
| 06 | `06_compute_gene_region_methylation.py` | CpG数×overlap比例加权，汇总单细胞和cell type甲基化 |
| 07 | `07_summarize_marker_expression.py` | 汇总marker表达均值、表达比例和逐基因Z-score |
| 08 | `08_merge_plot_joint_bubbles.py` | 合并长表并绘制三个并排气泡图 |

总控：`run_meth_expr_bubble_workflow.sh`。

Marker默认条件：adjusted P value `<0.05`、log2FC `>0.25`、目标细胞表达比例 `>=10%`，只保留正向marker。一个gene symbol如果出现在多个cell type，只归入log2FC最大的cell type，然后每类最多取Top5。

## 核心定义

DMR-cell mean ratio沿用原流程：DMR内所有唯一CpG ratio的等权算术平均。

区域汇总权重：

```text
effective_weight = unique_CpG_count(cell, DMR) × overlap_bp / DMR_length
region_probability = Σ(DMR_probability × effective_weight) / Σ(effective_weight)
```

Promoter固定为链特异的`TSS -2000 bp / +500 bp`；gene body使用GENCODE gene feature完整区间。DMR可同时属于promoter和gene body，也可对应多个基因。

甲基化汇总要求至少10个有效细胞且有效比例至少20%；否则均值保留为`NA`。`NA`不补0。

最终结果只能解释为“DMR覆盖部分所支持的区域甲基化概率”，不能称为完整promoter或gene-body CpG平均甲基化。

## 运行

```bash
cd /share/home/rzli/METHSCAN/03_MethExprBubble
chmod 750 run_meth_expr_bubble_workflow.sh

bash run_meth_expr_bubble_workflow.sh audit
bash run_meth_expr_bubble_workflow.sh markers
bash run_meth_expr_bubble_workflow.sh regions
bash run_meth_expr_bubble_workflow.sh counts 64
bash run_meth_expr_bubble_workflow.sh map
bash run_meth_expr_bubble_workflow.sh methylation
bash run_meth_expr_bubble_workflow.sh expression
bash run_meth_expr_bubble_workflow.sh plot
```

不要在未检查`01_audit/`和`02_markers/marker_genes.tsv`前直接运行全部步骤。

完整运行可以提交：

```bash
cd /share/home/rzli/METHSCAN/03_MethExprBubble
mkdir -p scheduler_logs

dsub \
  -n meth_expr_bubble_all_64c \
  -R "cpu=64;mem=180G" \
  --cwd /share/home/rzli/METHSCAN/03_MethExprBubble \
  -oo scheduler_logs/meth_expr_bubble_all_64c.out \
  -eo scheduler_logs/meth_expr_bubble_all_64c.err \
  bash /share/home/rzli/METHSCAN/03_MethExprBubble/run_meth_expr_bubble_workflow.sh \
  all 64
```

第04步对每个样本依次处理、样本内部按cell滚动并行；任一cell完成后会立即补入下一个cell，不会等待固定批次全部结束。

## 主要输出

- `results/02_markers/marker_genes.tsv`：最终全局唯一Top5 marker
- `results/03_gene_regions/marker_gene_promoters.bed`
- `results/03_gene_regions/marker_gene_bodies.bed`
- `results/04_dmr_unique_cpg_counts/*/*unique_CpG_count.tsv.gz`
- `results/05_dmr_gene_region_map/dmr_marker_gene_region_overlaps.tsv.gz`
- `results/06_region_methylation/single_cell_gene_region_methylation.tsv.gz`
- `results/06_region_methylation/celltype_gene_region_methylation.tsv`
- `results/07_expression/celltype_marker_gene_expression.tsv`
- `results/08_joint_plot/joint_methylation_expression_bubble_data.tsv.gz`
- `results/08_joint_plot/methylation_expression_top5_marker_bubble_plot.png`
- `results/08_joint_plot/methylation_expression_top5_marker_bubble_plot.pdf`
