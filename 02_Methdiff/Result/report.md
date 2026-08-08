# MethSCAn 细胞类型 hypo-DMR top1500 筛选与样本级合并报告

报告日期：2026-08-06

## 1. 分析目的

从 MethSCAn 样本内细胞类型两两比较结果中提取各细胞类型相对其他细胞类型的 hypo-DMR。每个样本、每种细胞最多保留 `|diff|` 最大的 1500 个唯一 DMR，再在样本内去重并合并重叠区间，用于构建单细胞甲基化热图。

## 2. hypo-DMR top1500 筛选条件

- MethSCAn 结果第 11 列 `raw p < 0.01`；
- 第 8、9 列甲基化比例差的绝对值 `>= 0.25`；
- 第 10 列 `lower_methylated_group` 指向的组作为 hypo 细胞类型；
- 当第 10 列为 `group_A` 时，A 为 hypo；为 `group_B` 时，B 为 hypo；
- 仅保留 `chr1–chr22、chrX、chrY`；
- 每个样本、每种 hypo 细胞汇总所有 pairwise DMR，完全相同的坐标先去重；
- 按 `|column 8-column 9|` 降序保留前 1500 个唯一 DMR；通过基础条件不足 1500 个时保留全部。

筛选脚本：`05_extract_celltype_hypo_dmrs_top1500.py`

## 3. 样本级 DMR 合并规则

合并脚本 `02_merge_sample_dmrs.py` 读取每个样本 `pairwise_union/` 中的所有 top1500 hypo-DMR，并按以下规则处理：

1. 完全相同的 DMR 坐标只保留一次。
2. 位于同一染色体且真正重叠的 DMR 递归合并为一个新 DMR。
3. BED 区间按左闭右开区间解释；仅首尾相接（`next.start == current.end`）的区间不属于重叠，因此不合并。
4. 各样本独立合并，不跨样本合并 DMR。
5. 保留每个合并 DMR 的来源区间数、来源比较文件数及支持的 hypo 细胞类型。

## 4. 计算方式

两个步骤均在集群计算节点上通过 `dsub` 运行：

- hypo-DMR top1500 筛选：10 个样本并行；
- 样本级 DMR 合并作业：Job ID `163282`；
- 样本级并行：10 个进程，每个样本由一个独立进程处理；
- 筛选作业申请 10 CPU、20 GiB 内存；合并作业申请 10 CPU、4 GiB 内存。

## 5. 合并结果

| 样本 | 非空 pairwise union 文件 | 合并前 DMR | 合并后 DMR | 减少数量 | 减少比例 |
|---|---:|---:|---:|---:|---:|
| IR01 | 116 | 23,992 | 11,936 | 12,056 | 50.25% |
| IR02 | 128 | 27,962 | 12,502 | 15,460 | 55.29% |
| IR03 | 109 | 25,170 | 11,254 | 13,916 | 55.29% |
| IR04 | 106 | 18,414 | 9,616 | 8,798 | 47.78% |
| IR05 | 107 | 17,797 | 9,366 | 8,431 | 47.37% |
| NR01 | 92 | 20,005 | 10,242 | 9,763 | 48.80% |
| NR02 | 68 | 17,951 | 9,610 | 8,341 | 46.47% |
| NR03 | 95 | 23,852 | 11,984 | 11,868 | 49.76% |
| NR04 | 117 | 25,215 | 11,857 | 13,358 | 52.98% |
| NR05 | 50 | 12,402 | 8,129 | 4,273 | 34.45% |
| **总计** | **988** | **212,760** | **106,496** | **106,264** | **49.95%** |

## 6. 结果概要

- 10 个样本均成功产生样本级合并 DMR。
- 合并前共有 212,760 个来源区间，合并后得到 106,496 个样本内互不重叠的 DMR。
- 去重和合并共减少 106,264 个区间，整体减少 49.95%。
- 合并后 DMR 数最多的样本为 IR02（12,502），最少的样本为 NR05（8,129）。
- 各样本的 pairwise union 文件数不同，反映了样本内细胞类型组成以及通过筛选的非空比较数存在差异。

## 7. 结果文件

node4 上的结果目录：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p25_top1500
```

每个样本产生两个主要文件：

- `<sample>__merged_DMRs.bed`：无表头的标准三列 BED，保存合并后互不重叠的 DMR 坐标；
- `<sample>__merged_DMRs_annotation.tsv`：保存 DMR ID、来源区间数、来源文件数和支持的 hypo 细胞类型。

总结文件：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p25_top1500/merge_summary.tsv
```

## 8. 后续单细胞甲基化矩阵

每个样本的 `<sample>__merged_DMRs.bed` 作为矩阵的 DMR 行坐标集合。`03_compute_dmr_mean_cpg_ratio.py` 计算每个单细胞在各合并 DMR 中的 CpG 甲基化比例，形成：

```text
行：合并后 DMR
列：单细胞
值：该细胞在该 DMR 中的甲基化 CpG 调用数 / 有覆盖 CpG 调用数
```

某细胞在某 DMR 内无覆盖时记为 `NA`。对稀疏矩阵中由正负链调用合并得到的 `-2,-1,0,1,2`，分别按 `0/2,0/1,1/2,1/1,2/2` 恢复甲基化与总调用数。另外生成 `cell_annotations.tsv` 保存单细胞与细胞类型的对应关系。当前每个样本包含 8,129–12,502 个合并 DMR，每个单细胞单独一列，输出为 gzip 压缩宽表。
