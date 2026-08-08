# 单细胞甲基化 DMR 的 GREAT 分析报告

## 一、完成情况

已完成 CD14 Monocytes 与 NK cells 的四组 DMR 预处理、质量核验和独立富集分析：

| 分析组 | 最终区域数 | 关联基因数 | 通过主筛选的 GO 条目 | Binomial 与 Hypergeometric 均显著 |
|---|---:|---:|---:|---:|
| CD14 IR-hypo | 849 | 1,308 | 8 | 3 |
| CD14 IR-hyper | 563 | 802 | 1 | 0 |
| NK IR-hypo | 9,938 | 6,571 | 125 | 11 |
| NK IR-hyper | 103 | 178 | 0 | 0 |

四个集合均不少于 20 个区域，因此都进行了富集分析。Hyper/Hypo 和细胞类型始终分开，没有混合。

## 二、数据预处理

处理顺序与预定方案一致：

1. 从原始 12 列 MethSCAn DMR 中保留第 12 列 `adjusted_p < 0.05`。
2. 仅保留 `chr1–chr22、chrX、chrY`。
3. 与对应 `*_clean.bed` 相交。
4. 第 10 列 `group_A` 归为 `IR_hypo`，`group_B` 归为 `IR_hyper`。
5. 同一细胞类型和方向内合并重叠区域。
6. 输出无表头 BED4，区域 ID 唯一。

各阶段区域数：

| 细胞类型 | 原始 DMR | q<0.05 | 常染色体/X/Y | 通过 clean 相交 | IR-hypo | IR-hyper |
|---|---:|---:|---:|---:|---:|---:|
| CD14 Monocytes | 88,500 | 1,996 | 1,984 | 1,412 | 849 | 563 |
| NK cells | 87,413 | 10,061 | 10,045 | 10,041 | 9,938 | 103 |

这里的 Q 值为第 12 列 `adjusted_p`，原始 P 值为第 11 列
`raw_p`。以下均为限制到 `chr1–chr22、chrX、chrY` 后的区域数；
IR-hypo 和 IR-hyper 是在 Q<0.05 的集合中进行方向拆分：

| 细胞类型 | 原始P<0.05 | Q<0.05 | IR-hypo | IR-hyper |
|---|---:|---:|---:|---:|
| CD14 Monocytes | 87,976 | 1,984 | 1,165 | 819 |
| NK cells | 60,120 | 10,045 | 9,941 | 104 |
| CD4 T cells | 89,043 | 26 | 16 | 10 |
| Cycling cells | 22,193 | 3 | 0 | 3 |
| MAIT cells | 60,846 | 3 | 0 | 3 |
| CD16 Monocytes | 86,385 | 2 | 1 | 1 |
| CD8 T cells | 93,466 | 2 | 1 | 1 |
| Treg cells | 44,151 | 2 | 0 | 2 |
| Plasma cells | 52,652 | 3 | 2 | 1 |
| B cells | 92,944 | 1 | 1 | 0 |
| HLAII-high APCs | 38,669 | 1 | 0 | 1 |
| B cells unresolved | 6,401 | 0 | 0 | 0 |
| Gamma-delta T cells | 7,471 | 0 | 0 | 0 |
| cDCs | 4 | 0 | 0 | 0 |
| pDCs | 4,649 | 0 | 0 | 0 |

该五列表保存在
`GREAT_inputs/all_celltypes_p_q_hypo_hyper.tsv`。另外，
`raw_p_counts_all_celltypes.tsv` 保留了原始 P 条件下的方向拆分，
`adjusted_p_counts_all_celltypes.tsv` 保留了染色体过滤前后的 Q 值计数。

## 三、富集方法

Stanford GREAT 网页在执行时不可访问，因此采用本地 `rGREAT` 完成分析。四组使用相同设置：

- hg38；
- whole-genome background；
- basal plus extension；
- upstream 5 kb；
- downstream 1 kb；
- maximum extension 1,000 kb；
- UCSC hg38 known-gene TSS；
- 排除 UCSC assembly gaps；
- GO Biological Process、Cellular Component 和 Molecular Function。

本地模式不包含 Stanford 服务器的 curated regulatory domains，这是与原网页方案的唯一重要参数差异之一；网页恢复后应使用四个 BED 再做一次严格复现。

主结果筛选标准：

- Binomial FDR < 0.05；
- region fold enrichment ≥ 2；
- observed region hits ≥ 5；
- observed gene hits ≥ 3。

优先结果进一步要求 Hypergeometric FDR < 0.05。

## 四、主要结果

### 1. CD14 IR-hypo

共有 8 个条目通过主筛选，其中 3 个同时通过两种 FDR：

| 代表条目 | Binomial FDR | Hypergeometric FDR | 区域富集倍数 | Region hits | Gene hits | 集中度提示 |
|---|---:|---:|---:|---:|---:|---|
| embryonic placenta morphogenesis | 0.0445 | 0.0266 | 3.78 | 13 | 8 | 否 |
| cell-cell signaling involved in cardiac conduction | 0.0445 | 0.00554 | 3.16 | 16 | 10 | 否 |
| chorio-allantoic fusion | 0.0131 | 0.00491 | 6.59 | 10 | 5 | 是，DNAJB6 占 60% |

这些条目可归为“胚外/胎盘发育”和“心脏传导信号”两个主题，但并不是典型的单核细胞免疫功能。尤其 `chorio-allantoic fusion` 有明显的 DNAJB6 位点集中，不能单独作为生物学结论。

### 2. CD14 IR-hyper

只有 `mitochondrial crista` 通过区域主筛选：

- Binomial FDR = 0.0129；
- region fold enrichment = 13.38；
- 5 个 region hits、3 个 gene hits；
- Hypergeometric FDR = 0.128，未显著；
- OPA1 相关区域占 60%。

因此不列为优先通路，只可作为待验证的区域型信号。

### 3. NK IR-hypo

共有 125 个条目通过主筛选，11 个同时通过两种 FDR，合并后形成 7 个主题：

- 神经与胚胎模式形成；
- 嘌呤能/腺苷 GPCR 信号；
- 组织修复与再生；
- 神经递质分泌；
- 内胚层分化；
- 心肌膜去极化；
- 线粒体凋亡。

其中嘌呤能和腺苷受体的两个 GO 条目统计量完全相同，按一个主题报告。`negative regulation of glutamate secretion` 主要由 GRM7 位点驱动（17/23，73.9%）；骨骼肌卫星细胞维持条目主要由 WNT7A 位点驱动（12/19，63.2%）。

NK IR-hypo 的区域数达到 9,938，关联 6,571 个基因；125 个主条目中有 48 个出现单基因位点贡献至少 50% 的集中度提示。虽然 11 个条目通过了双检验，但整体主题大多不是典型 NK 免疫通路，提示结果可能受大区域集、局部 DMR 簇和当前细胞级统计设计影响。

### 4. NK IR-hyper

没有 GO 条目达到主筛选标准，因此不报告通路富集。区域—基因关联表仍已输出，可用于查看候选基因。

## 五、总体判断

当前结果适合作为探索性结果，不支持直接下结论说 IR 与 CD14/NK 的经典免疫激活通路显著相关。最可靠的表述是：

> 在去除样本效应并按 adjusted p-value 筛选后，CD14 IR-hypo 与 NK IR-hypo 出现若干区域和基因双重显著的功能条目，但主题以发育、神经、心脏和组织修复为主，且部分由少数基因座附近的 DMR 簇驱动；需要供体层面 DMR 复算后验证。

不要把仅 Binomial FDR 显著的条目与双检验显著条目混在一起，也不要把这些非典型主题改写成“免疫激活”或“NK 功能增强”。

## 六、正式分析建议

论文主结果应按“样本 × 细胞类型”形成 5 个 IR 和 5 个 NR pseudobulk，再使用 DSS 或 bsseq 重新计算 DMR。建议正式阈值：

- FDR < 0.05；
- |Δ甲基化比例| ≥ 0.10；
- 至少 3 个 CpG；
- 每组至少 3/5 样本有效覆盖。

将新 DMR 继续按 IR-hypo/IR-hyper 分开，以相同 GREAT 参数运行。只有供体层结果与当前结果在方向和功能主题上都一致时，才作为稳定结论。

## 七、结果文件说明

- `GREAT_inputs/`：四个最终 BED 和源 DMR 审计记录。
- `GREAT_results/run_summary.tsv`：四组总体统计。
- `GREAT_results/primary_summary.tsv`：全部 134 个主筛选条目。
- `GREAT_results/significant_by_both_summary.tsv`：14 个双检验显著条目及命中基因。
- `GREAT_results/theme_summary.tsv`：去冗余后的 9 个主题。
- `GREAT_results/*/region_gene_associations.tsv`：区域到基因。
- `GREAT_results/*/gene_region_associations.tsv`：基因到区域。
- `GREAT_results/*/tss_distance_summary.tsv`：区域—TSS 距离分布。

`hit_count_translation_audit.tsv` 记录了 21 个条目在 Entrez ID 转基因符号后出现的轻微 region-hit 数量差异。筛选和报告中的 region hits 始终采用 `rGREAT` 原始统计值，基因符号表用于定位和解释。
