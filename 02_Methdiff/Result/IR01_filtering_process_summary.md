# IR01 MethSCAn 筛选与 DMR/Top200 热图流程数据汇总

## 1. 分析范围与最终口径

- 样本：`25110891_IR01_Met`
- 输入：7,981 个单细胞 Bismark cov 文件
- 最终细胞 QC：每细胞至少 300,000 个观测甲基化位点，全基因组平均甲基化率至少 55%
- DMR 范围：仅 `chr1–chr22、chrX、chrY`
- DMR 比较：IR01 样本内细胞类型两两比较
- `methscan diff --min-cells 10`：待检验区域在每组至少有 10 个细胞覆盖
- Top200：`raw p < 0.01`、`abs(meth_frac_A - meth_frac_B) >= 0.25`，每个目标细胞类型按甲基化差异绝对值排序取低甲基化 DMR top200
- 热图数值：每个 `cell × DMR` 内所有唯一 CpG ratio 的等权算术平均

30k 结果保留作为探索/对照；300k 结果为当前最终口径。

## 2. cov 重复 CpG 审计

### 2.1 原始 cov

| 指标 | 数值 |
|---|---:|
| cov 文件 | 7,981 |
| 读取正常文件 | 7,981 |
| 读取错误 | 0 |
| 总行数 | 1,193,992,170 |
| 列数异常行 | 0 |
| 含重复 CpG 的文件 | 7,974 |
| 重复 CpG 坐标 | 81,830,667 |
| 额外重复行 | 81,830,667 |
| 6 列完全一致的重复坐标 | 74,275,574 |
| 同坐标但至少一列不同 | 7,555,093 |
| 主染色体重复坐标 | 81,532,939 |
| 非主染色体重复坐标 | 297,728 |
| 坐标无序文件 | 0 |

### 2.2 最终去重规则

CpG 坐标定义为 cov 第 1–3 列：`chrom, start, end`。

1. 同坐标且甲基化概率（第 4 列）相同：任意保留第一行。
2. 同坐标但甲基化概率不同：删除该坐标的全部行。
3. 唯一坐标：原样保留。

| 去重结果 | 数值 |
|---|---:|
| 处理文件 | 7,981 |
| 失败文件 | 0 |
| 输入行 | 1,193,992,170 |
| 输出行 | 1,108,218,648 |
| 保留一行的同概率重复坐标 | 77,887,812 |
| 整个删除的不同概率重复坐标 | 3,942,855 |
| 删除行 | 85,773,522（7.18%） |
| 保留行比例 | 92.82% |

去重后 QC：

| 指标 | 数值 |
|---|---:|
| 总行数 | 1,108,218,648 |
| 重复 CpG 坐标 | 0 |
| 列数异常行 | 0 |
| 无序文件 | 0 |

去重 cov：

```text
/share/LCZX_Data/data/allcools/25110891_IR01_Met/cov_dedup_probability
```

## 3. MethSCAn prepare/profile 与全局 QC

- compact：`compact_data_dedup_probability`
- TSS BED：`/share/LCZX_Data/ref/human_hg38_TSS.bed`
- TSS BED SHA-256：`cabd65c85ce0db017d771744e9db6ea80b3ef741594763730d6917756914d631`
- profile 区域：42,017
- 在所有细胞中都未观测的 profile 区域：594（1.41%）

`min_meth=55` 是每个细胞的全基因组平均甲基化百分比下限，不是 CpG 数量或 DMR 内覆盖门槛。

| 阶段 | min_sites | min_meth | 输入细胞 | 保留细胞 | 相对 7,981 保留率 |
|---|---:|---:|---:|---:|---:|
| 30k 探索 | 30,000 | 55% | 7,981 | 6,855 | 85.89% |
| 300k 最终 | 300,000 | 55% | 7,981 | 583 | 7.30% |

300k 相对 30k 仅保留 8.50% 的细胞。30k 与 300k 分别写入 `filtered_data_single_30k` 和 `filtered_data_single_300k`，互不覆盖。

## 4. DMR 输入限制

- 仅使用 24 条主染色体：`chr1–chr22, chrX, chrY`
- 排除 `chrM`、`KI*`、`GL*` 及其他非主染色体
- 30k 输入视图审计：24 条主染色体纳入，159 条非主 contig 排除，空主染色体0
- 同一样本内细胞类型两两比较
- 组总细胞数小于10的比较跳过
- 已运行比较使用 `methscan diff --min-cells 10`

### 4.1 30k 探索 DMR

| 指标 | 数值 |
|---|---:|
| QC 后细胞 | 6,855 |
| 未匹配 Scanpy 注释 | 97 |
| 进入最终热图的已注释、非排除细胞 | 6,746 |
| 非排除细胞类型 | 15 |
| 全部两两组合 | 105 |
| 有效比较 | 78（74.29%） |
| 不合格比较 | 27 |
| 完成比较 | 78 |

### 4.2 300k 最终 DMR

| 指标 | 数值 |
|---|---:|
| QC 后细胞 | 583 |
| 进入热图的已注释、非排除细胞 | 560（96.05%） |
| 参与 DMR 组合的非排除细胞类型 | 13 |
| 总比较 | 78 |
| 达到至少10个细胞的类型 | 10 |
| 有效并完成的比较 | 45（57.69%） |
| 不合格比较 | 33 |

300k 中达到10细胞门槛的类型：

| 细胞类型 | 细胞数 |
|---|---:|
| CD4_T_cells | 136 |
| CD14_Monocytes | 120 |
| CD8_T_cells | 82 |
| NK_cells | 67 |
| CD16_Monocytes | 38 |
| B_cells | 31 |
| HLAII_high_APCs | 31 |
| MAIT_cells | 18 |
| Treg_cells | 12 |
| pDCs | 10 |

未达到10细胞门槛的类型：

| 细胞类型 | 细胞数 |
|---|---:|
| Cycling_cells | 6 |
| Plasma_cells | 6 |
| cDCs | 3 |

注：pDCs 只有10个细胞，而区域级 `min-cells=10`，因此涉及 pDCs 的区域检验最严格。

## 5. Top200 hypo-DMR 与热图筛选

Top200 采用原始 p 值而非 adjusted p/FDR：

```text
raw p < 0.01
abs(meth_frac_group_A - meth_frac_group_B) >= 0.25
仅主染色体
按低甲基化目标细胞类型分组
每个目标细胞类型最多200个唯一DMR区间
```

| 指标 | 30k 探索 | 300k 最终 |
|---|---:|---:|
| DMR BED 输入行 | 1,232,045 | 356,581 |
| 通过 raw-p/甲基化差异门槛的行 | 340,363（27.63%） | 85,198（23.89%） |
| 各目标细胞类型 top200 后的唯一区间 | 2,216 | 1,584 |
| pairwise union 来源区域 | 4,743 | 2,143 |
| 去重并合并真正重叠后的 DMR | 2,004 | 1,516 |

同一 DMR 可在多个 pairwise 比较中支持同一低甲基化细胞类型，因此 `pairwise union 来源区域` 可大于各类型 top200 唯一区间数。合并后热图列数不必等于 `细胞类型数 × 200`。

## 6. mean CpG ratio 定义与覆盖审计

对每个唯一 CpG：

```text
CpG ratio_i = meth_i / (meth_i + unmeth_i)
```

对每个 `cell × DMR`：

```text
mean CpG ratio = (1 / K) * sum(CpG ratio_i), i = 1..K
```

- 不同唯一 CpG 等权，不按 read depth 加权。
- DMR 内无 cov CpG 时记为 `NA`。
- 直接读取 `cov_dedup_probability`，不使用 NPZ 二元调用或 smoothed 值。
- DMR 使用 BED 半开区间：`start <= CpG position < end`。

| 指标 | 30k 探索 | 300k 最终 |
|---|---:|---:|
| 热图细胞 | 6,746 | 560 |
| 热图 DMR | 2,004 | 1,516 |
| 总 `cell × DMR` | 13,518,984 | 848,960 |
| 有至少1个 CpG 覆盖的值 | 739,892 | 180,840 |
| 覆盖率 | 5.47% | 21.30% |
| 缺失率 | 94.53% | 78.70% |
| 纳入计算的 cov 行 | 2,700,950 | 681,565 |
| 唯一 CpG | 2,700,950 | 681,565 |
| DMR 内重复 CpG 组 | 0 | 0 |
| DMR 内冲突重复 CpG 组 | 0 | 0 |
| 每个非 NA `cell × DMR` 平均 CpG 数 | 3.65 | 3.77 |

300k 的 `cell × DMR` 覆盖率是 30k 的约 3.9 倍，但细胞数从 6,746 下降到560。

## 7. 关键输出路径

### 30k 探索结果

```text
/share/LCZX_Data/data/allcools/25110891_IR01_Met/
  qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/
    methdiff_celltype_30k/
```

### 300k 最终结果

```text
/share/LCZX_Data/data/allcools/25110891_IR01_Met/
  qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/
    methdiff_celltype_300k/
```

### 300k 最终 Top200 mean-CpG-ratio 热图

```text
/share/LCZX_Data/data/allcools/25110891_IR01_Met/
  qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/
    methdiff_celltype_300k/
      heatmap_top200_rawp0p01_diff0p25/
        figures_top200_mean_of_unique_CpG_ratios/
          IR01/
            IR01__cells_by_all_DMRs_grouped_heatmap.png
```

Heatmap 目录软链接：

```text
/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/IR01_300k_top200_results
```

## 8. 主要解读注意事项

1. 300k 是高覆盖细胞子集，只保留原始 7,981 个细胞的 7.30%，可能引入细胞类型构成偏倚。
2. Top200 当前基于 raw p，而非 adjusted p/FDR；用于正式统计结论时需明确这一点。
3. pDCs 总数恰好为10，区域级 `min-cells=10` 对其比较尤其严格。
4. 300k 单细胞热图仍有 78.70% `NA`，这是原始单细胞 CpG 覆盖的稀疏性，不是去重或均值计算错误。
5. DMR 是基于组比较发现；单个 DMR 无需在热图的所有细胞中都有覆盖。
