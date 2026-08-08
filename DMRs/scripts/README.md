# Meth_diff DMR 分析脚本说明

本目录包含两套使用相同分析流程、但显著性阈值不同的脚本：

```text
q001_2kb_2kb/    q value < 0.01（高置信度、严格阈值）
q005_2kb_2kb/    q value < 0.05（候选发现、相对宽松阈值）
```

两套目录均按 `00` 到 `12` 编号，文件名相同；具体 q 值由所在目录决定。`common/` 中保存不依赖 q 阈值的共享资源脚本。

## 分析目标与分组

```text
IR = immunotherapy responder，肺癌免疫治疗应答者
NR = immunotherapy non-responder，肺癌免疫治疗不应答者
```

分析分为两条不能混用的线：

```text
A. clean cell-type DMR / VMR matrix 分析（01–03）
   用于构建较少受 IR/NR 应答状态影响的细胞类型甲基化特征。

B. IR vs NR 应答相关 DMR、promoter、RNA 与 DNA-RNA 相关性分析（04–12）
   用原始 same-cell-type IR vs NR DMR 发现应答相关候选基因。
```

特别注意：步骤 A 中会从 cell-type DMR 中排除与 IR/NR DMR 重叠的区域。因此，步骤 B 必须直接使用原始的 same-cell-type IR vs NR DMR，不能使用 clean cell-type DMR，否则会删除真正的应答相关差异。

## 总流程

```text
细胞类型两两 DMR ─┐
                   ├─ 01：减去 IR/NR 重叠区域 → clean cell-type DMR
同细胞类型 IR vs NR ┘
                                      ↓
                         02–03：抽取对应 clean VMR matrix

同细胞类型 IR vs NR DMR
        ↓
04：统计 DMR 数量
        ↓
05：按 q 值筛选，并区分 IR-hyper / IR-hypo
        ↓
06：调用共享脚本生成 promoter 区域（TSS ±2 kb）
        ↓
07：DMR 与 promoter 坐标 overlap
        ↓
08：整理 DMR–gene 注释表
        ↓
09：保留 protein-coding genes
        ↓
10–11：与 RNA pseudobulk expression 整合
        ↓
保留 promoter 甲基化与 RNA 表达反向变化的候选基因
        ↓
12：以匹配 sample × cell type 的 DNA/RNA pseudobulk 计算相关性
```

## 每个脚本的作用

| 步骤 | 文件名 | 作用 | 主要输出 |
|---:|---|---|---|
| 00 | `00_make_dna_cell_metadata.py` | 读取每个 cell type 的 Meth_diff group CSV；将 `group_A/group_B` 转换为 `IR/NR`，从 `cell_id` 的 `__` 前缀提取 sample ID。 | DNA cell metadata：`cell_id`、`sample`、`response`、`cell_type`。 |
| 01 | `01_merge_subtract_dmr.sh` | 合并显著 cell-type pairwise DMR 与显著 IR-vs-NR DMR；从前者中去除与后者重叠的区域。 | clean cell-type DMR BED；对应 VMR region list。 |
| 02 | `02_subset_vmr_matrix.py` | 根据 VMR region list，从完整 All VMR matrix 中按列提取目标 VMR。 | 四类子矩阵：methylation fraction、methylated sites、total sites、mean shrunken residuals。 |
| 03 | `03_run_subset_vmr_matrix.sh` | 配置 Conda/Python 环境并调用步骤 02。 | 步骤 02 的矩阵输出目录。 |
| 04 | `04_summarize_ir_vs_nr_dmr.sh` | 统计每种细胞类型的总 DMR、显著 DMR、IR-hyper 与 IR-hypo DMR 数量。 | IR-vs-NR DMR summary TSV。 |
| 05 | `05_filter_ir_vs_nr_dmr.sh` | 按当前目录对应 q 阈值筛选 DMR；根据 Meth_diff direction 列和 group_A=IR、group_B=NR 赋予方向。 | 每个 cell type 的 q 值筛选 DMR、IR-hyper BED、IR-hypo BED。 |
| 06 | `06_make_promoter_bed_2kb_2kb.py` | 调用 `common/make_promoter_bed_2kb_2kb.py`，从 GENCODE 注释生成共享 promoter BED。promoter 定义为 TSS 上游 2 kb 至下游 2 kb。 | 共享 GENCODE promoter BED。 |
| 07 | `07_intersect_dmr_promoter.sh` | 将 IR-hyper/IR-hypo DMR 与 promoter BED 做基因组坐标 overlap。 | 每个 cell type、每个方向的 promoter DMR overlap 文件与数量汇总。 |
| 08 | `08_format_promoter_dmr_gene_table.py` | 将 bedtools overlap 输出整理为结构化 DMR–gene 表，并保留 DMR 统计量、IR/NR 平均甲基化值、基因信息和 promoter 信息。 | `promoter_DMR_to_gene` 主表。 |
| 09 | `09_filter_protein_coding_genes.py` | 从 DMR–gene 主表中筛选 protein-coding genes，并按 cell type、方向、gene 去重。 | protein-coding DMR–gene 表及候选 gene list。 |
| 10 | `10_integrate_pseudobulk_expression.py` | 对候选基因计算同一 cell type 内 IR 与 NR 的 RNA pseudobulk 平均表达，并与 promoter 甲基化差异整合；记录 `adata.raw.X` 的来源和未自动判定的归一化状态。 | 含 methylation delta、expression delta 与方向标签的候选表；RNA run metadata。 |
| 11 | `11_run_pseudobulk_expression.sh` | 配置运行环境并调用步骤 10。 | 步骤 10 的 RNA integration 结果。 |
| 12 | `12_correlate_dna_methylation_rna_expression.py` | 对每个样本、每个 cell type 聚合 DMR 重叠 VMR 的 DNA 甲基化（`sum(methylated_sites)/sum(total_sites)`）；仅分块读取候选 DMR 所需 VMR 列，再与匹配 RNA pseudobulk 均值计算相关性。 | DNA/RNA 配对表、相关性/FDR 表、诊断表和排除记录。 |

## 方向解释

Meth_diff 输出中第 10 列表示低甲基化的分组。当前 cell group 定义为：

```text
group_A = IR
group_B = NR
```

因此：

```text
第 10 列为 group_B
→ NR 低甲基化
→ IR 高甲基化
→ IR_hyper

第 10 列为 group_A
→ IR 低甲基化
→ IR_hypo
```

## RNA 整合后的候选筛选

优先保留 promoter methylation 与表达呈反向方向的候选：

```text
IR_hyper promoter + IR expression down
IR_hypo promoter  + IR expression up
```

进一步可采用效应大小筛选：

```text
|delta_methylation_IR_minus_NR| >= 0.03
|delta_expression_IR_minus_NR| >= 0.005
```

以上阈值用于候选优先级排序，不是通用固定标准；表达阈值应结合 RNA 矩阵的归一化尺度解释。

## 运行与维护注意事项

1. 各脚本中仍可能保留服务器上的绝对路径（如 `/share/home/rzli/...`）。将本地整理后的脚本复制回服务器或直接运行前，请先检查这些路径。
2. q001 与 q005 目录中的 `01`、`03` 脚本分别写入对应 `q001` 或 `q005` 输出路径；不要混用两个目录的中间结果。
3. `02` 使用 VMR region 名称与矩阵列名的精确字符串匹配；region list 与矩阵列名必须一致。
4. promoter overlap 表示位置关联，结合 RNA 反向变化可优先筛选候选，但不能单独证明甲基化对表达的因果调控。
5. 步骤 06 的共享 promoter BED 位于 `Meth_diff/common_resources/promoter_annotation/`，只需生成一次；两套阈值的步骤 07 都读取该文件。
6. 步骤 00 可以从 Meth_diff cell-group CSV 生成步骤 12 所需 DNA metadata；`group_A=IR`、`group_B=NR`，而标为 `-` 或空值的细胞会被排除。若同一 `cell_id` 出现冲突 response 或 cell type，脚本会停止。
7. 步骤 12 需要 DNA metadata（`cell_id`、`sample`、`response`、`cell_type`）以及 methylated-sites 和 total-sites 矩阵；DNA 与 RNA 必须使用可匹配的 sample ID。pooled correlation 仅作描述性结果，优先解读 IR 或 NR 组内相关性及 FDR。
8. 正式运行步骤 10 前，确认 `adata.raw.X` 是希望使用的表达尺度；脚本会写出 `RNA_pseudobulk_run_metadata.tsv`，但不会自行推断 raw matrix 是否已归一化。
