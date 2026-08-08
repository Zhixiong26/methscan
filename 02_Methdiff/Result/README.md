# MethSCAn 细胞类型 DMR 与单细胞热图分析

## 分析流程

本目录保留 4 个主流程脚本，以及 top-N DMR 和 CpG 等权平均的独立可选版本。

```text
MethSCAn 样本内细胞类型两两比较 BED
        │
        ▼
01_extract_celltype_hypo_dmrs.py
  筛选 raw p < 0.01、|col8-col9| >= 0.30
  根据第 10 列确定 hypo 细胞类型，生成 pairwise union
        │
        ▼
02_merge_sample_dmrs.py
  在每个样本内去重并合真正重叠的 DMR
        │
        ▼
03_compute_dmr_mean_cpg_ratio.py
  读取 MethSCAn chromosome CSR 稀疏矩阵
  计算每个 DMR × 每个单细胞的 mean CpG ratio
        │
        ▼
04_plot_single_cell_dmr_heatmaps.py
  y 轴：所有单细胞，按 cell type 分组
  x 轴：全部 DMR，按上游差异分析的 hypo 细胞类型分组
  组内保留原始基因组顺序，不筛选 DMR
```

四步的 node4 默认结果依次为：

```text
methdiff_30k/results/celltype_hypo_DMRs_diff0p30
methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30
methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30
/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/figures_diff0p30_hypo_DMRs_grouped_5000bins_blue_red
```

## 目录中的脚本

| 顺序 | 脚本 | 作用 | 主要并行单位 |
|---:|---|---|---|
| 01 | `01_extract_celltype_hypo_dmrs.py` | 筛选各样本细胞类型两两比较 DMR | 10 个样本 |
| 02 | `02_merge_sample_dmrs.py` | 每样本内 DMR 去重和重叠区间合并 | 10 个样本 |
| 03 | `03_compute_dmr_mean_cpg_ratio.py` | 计算 DMR × single-cell mean CpG ratio | 24 条染色体 |
| 04 | `04_plot_single_cell_dmr_heatmaps.py` | 绘制 cell × all-DMR 分组热图 | 10 个样本 |
| 05 | `05_extract_celltype_hypo_dmrs_top1500.py` | 可选：每样本每种细胞保留 abs(diff) 最大的 1500 个 hypo DMR | 10 个样本 |
| 06 | `06_compute_dmr_mean_of_cpg_ratios.py` | 可选：先按坐标合并 cov 重复 CpG，再在 DMR 内对唯一 CpG ratio 等权取平均 | 10 个样本 |

`report.md` 保留早期运行的结果统计，不是运行脚本。本地的 `IR*/NR*_sample_celltype` 和长名结果目录为历史数据，不会被四个脚本改写。

## 上传到 node4

以下命令必须在 Mac 终端执行，不是在 node4 内执行：

```bash
cd "/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/02_Methdiff/Result"

scp \
  01_extract_celltype_hypo_dmrs.py \
  02_merge_sample_dmrs.py \
  03_compute_dmr_mean_cpg_ratio.py \
  04_plot_single_cell_dmr_heatmaps.py \
  05_extract_celltype_hypo_dmrs_top1500.py \
  06_compute_dmr_mean_of_cpg_ratios.py \
  README.md \
  rzli@node-4:/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/
```

node4 代码目录：

```text
/share/home/rzli/METHSCAN/02_Methdiff/Heatmap
```

## 01. 筛选细胞类型 hypo DMR

### 输入

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_celltype
```

输入文件例如：

```text
IR01__B_cells_vs_CD4_T_cells_DMRs.bed
```

MethSCAn DMR BED 的 12 列：

| 列 | 字段 |
|---:|---|
| 1 | `chrom` |
| 2 | `start` |
| 3 | `end` |
| 4 | `t_statistic` |
| 5 | `n_sites` |
| 6 | `n_cells_group_A` |
| 7 | `n_cells_group_B` |
| 8 | `meth_frac_group_A` |
| 9 | `meth_frac_group_B` |
| 10 | `lower_methylated_group` |
| 11 | `raw_p` |
| 12 | `adjusted_p` |

### 筛选条件

- 第 11 列 `raw p < 0.01`；
- 第 8、9 列 `abs(meth_frac_group_A - meth_frac_group_B) >= 0.30`；
- 第 10 列为 `group_A` 时，A 是 hypo；为 `group_B` 时，B 是 hypo；
- 默认仅保留 `chr1–chr22、chrX、chrY`。

该脚本只读原始 BED，不会修改或清空原始文件。

### dsub

```bash
cd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap
mkdir -p scheduler_logs

dsub \
  -n hypo_DMRs_diff0p30 \
  -R "cpu=10;mem=20G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_DMRs_diff0p30.out \
  -eo scheduler_logs/hypo_DMRs_diff0p30.err \
  python3 01_extract_celltype_hypo_dmrs.py \
  --jobs 10
```

### 输出

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/celltype_hypo_DMRs_diff0p30
```

主要内容：

- `by_sample/<sample>/pairwise/`：通过阈值的原始 12 列记录；
- `by_sample/<sample>/pairwise_union/`：同一 pairwise 文件内重叠区间的并集；
- `by_sample/<sample>/by_cell_type/`：按 hypo 细胞类型汇总；
- `overall_summary.tsv`：样本级统计；
- `parameters.tsv`：实际运行参数。

### 05 可选版本：|diff| >= 0.25 后每种细胞保留前 1500 个 hypo DMR

该版本使用 `raw p < 0.01`、`|column8-column9| >= 0.25` 和 hypo 方向筛选，再在每个样本的每种目标细胞内：

1. 汇总该细胞与所有其他细胞的 hypo DMR；
2. 对完全相同的 `chrom/start/end` 去重，排名时使用其最大 `abs(column8-column9)`；
3. 按 `abs(column8-column9)` 降序、`raw p` 升序和基因组坐标打破并列；
4. 保留前 1500 个唯一 DMR 区间。

`sample_summary.tsv` 中的 `selected_unique_hypo_DMR_intervals` 用于确认每种细胞的筛选数量。如果某种细胞通过基础条件的唯一 DMR 少于 1500，则保留全部。第 02 步进一步合并重叠区间后，最终独立区间数可能少于 1500。

第 01 步命令：

```bash
dsub \
  -n hypo_DMRs_diff0p25_top1500 \
  -R "cpu=10;mem=20G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_DMRs_diff0p25_top1500.out \
  -eo scheduler_logs/hypo_DMRs_diff0p25_top1500.err \
  python3 05_extract_celltype_hypo_dmrs_top1500.py \
  --jobs 10 \
  --top-dmrs-per-cell 1500
```

默认输出：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/celltype_hypo_DMRs_diff0p25_top1500
```

## 02. 每样本 DMR 去重与重叠合并

合并规则：

1. 完全相同的区间只保留一次；
2. 同一染色体上真正重叠的区间递归合并；
3. BED 按左闭右开理解，`next.start == current.end` 只是首尾相接，不合并；
4. 只在样本内合并，不跨样本合并。

```bash
dsub \
  -n merge_hypo_DMRs_diff0p30 \
  -R "cpu=10;mem=20G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/merge_hypo_DMRs_diff0p30.out \
  -eo scheduler_logs/merge_hypo_DMRs_diff0p30.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  02_merge_sample_dmrs.py \
  --jobs 10
```

输出：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30
```

- `<sample>__merged_DMRs.bed`：无表头三列 BED；
- `<sample>__merged_DMRs_annotation.tsv`：DMR ID、来源区间数、来源文件数和支持的 hypo 细胞类型；
- `merge_summary.tsv`：合并前后数量。

## 03. 计算每个单细胞在 DMR 上的 mean CpG ratio

### 输入

| 输入 | 默认路径 |
|---|---|
| chromosome CSR 稀疏矩阵 | `.../qc_minmeth55_maxmethnone_maxsites10000000/filtered_data_merged_30k/chr*.npz` |
| 单细胞 metadata | `.../methdiff_30k/metadata/cell_metadata.tsv` |
| 每样本合并 DMR | `.../methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30` |

MethSCAn 稀疏值解码：

| 稀疏值 | 甲基化 CpG | 总 CpG |
|---:|---:|---:|
| `-2` | 0 | 2 |
| `-1` | 0 | 1 |
| `0` | 1 | 2 |
| `1` | 1 | 1 |
| `2` | 2 | 2 |

对每个 `DMR × cell`：

```text
mean CpG ratio
    = DMR 内甲基化 CpG 数之和 / DMR 内总 CpG 数之和
```

当 DMR 内该细胞没有可用 CpG 调用时输出 `NA`。每个单细胞占一列，不对细胞类型取平均。

### dsub

```bash
dsub \
  -n hypo_DMR_cell_ratio_diff0p30 \
  -R "cpu=16;mem=80G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_DMR_cell_ratio_diff0p30.out \
  -eo scheduler_logs/hypo_DMR_cell_ratio_diff0p30.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  03_compute_dmr_mean_cpg_ratio.py \
  --jobs 16
```

输出：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30
```

- `<sample>__single_cell_DMR_mean_CpG_ratio.tsv.gz`：DMR × single-cell 宽表；
- `cell_annotations.tsv`：样本、单细胞和细胞类型映射；
- `matrix_summary.tsv`：DMR 数、单细胞数、覆盖值数和缺失比例；
- `parameters.tsv`：输入、解码和运行参数。

### 06 可选：合并重复记录后，对唯一 CpG ratio 等权平均

`06_compute_dmr_mean_of_cpg_ratios.py` 保留 `03` 的输出表格格式，但不读取 NPZ。脚本直接读取每个单细胞 `.cov.gz`，先按 `chrom + normalized CpG position` 合并同一 CpG 的正负链/重复记录：

```text
unique CpG ratio
    = sum(cov 第 5 列 methylated)
      / (sum(cov 第 5 列 methylated) + sum(cov 第 6 列 unmethylated))

DMR × cell
    = mean(DMR 内所有 unique CpG ratio)
```

因此，覆盖深度只用于整合同一 CpG 的多行记录；不同唯一 CpG 在 DMR 均值中仍各占一个相同权重。正负链 ratio 不一致时不会删除任一行，而是保留双方的 methylated/unmethylated 计数并合并成一个 CpG。脚本同时验证 cov 染色体和坐标顺序，并在 `matrix_summary.tsv` 中报告参与计算的原始 cov 行数、唯一 CpG 数、重复 CpG 组数及不一致重复组数。

默认读取 `diff0p25_top200` 合并 DMR，使用新的独立输出目录，不覆盖旧版逐行平均结果：

```text
.../methdiff_30k/results/single_cell_hypo_DMR_mean_of_unique_CpG_ratios_diff0p25_top200
```

node4 提交命令：

```bash
cd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap

dsub \
  -n hypo_top200_unique_CpG_ratios \
  -R "cpu=10;mem=40G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_top200_unique_CpG_ratios.out \
  -eo scheduler_logs/hypo_top200_unique_CpG_ratios.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  06_compute_dmr_mean_of_cpg_ratios.py \
  --jobs 10
```

## 04. 绘制 single-cell × all-DMR 分组热图

热图规则：

- y 轴每行为一个单细胞；
- 单细胞按 metadata `cell_type` 分组，组内保留原始顺序；
- x 轴使用全部合并 DMR，不做覆盖率、方差或 top-N 筛选；
- 从第 02 步 `<sample>__merged_DMRs_annotation.tsv` 读取 `supporting_hypo_cell_types`；
- 只有一个支持 hypo 类型时，DMR 直接归入该类型；
- 合并 DMR 同时支持多个 hypo 类型时，仅在这些上游支持类型中，选择实际 mean CpG ratio 最低者；并列或无观测值时归入 `Unresolved`；
- 同一 hypo 细胞类型的 DMR 放在连续区域，组内保留原始基因组顺序；
- 颜色只指定两个端点：`0=蓝色`、`1=红色`，中间由 Matplotlib 直接线性插值，没有白色中点；
- `NA=灰色`。

正式热图使用 5,000 个 display columns。4–12 万个 DMR 先按上游支持的 hypo 细胞类型分组，再在每个组内按原始基因组顺序分入显示列。每个 DMR 都参与一个显示列，不是筛选前 5,000 个 DMR。

```text
display_value(cell, DMR_group_bin)
    = mean(该单细胞在该 bin 内所有非 NA DMR ratio)
```

脚本仍保留 `--exact-dmr-columns` 作为可选模式，但数万个独立横向列在普通 PNG 和屏幕上无法分辨，不作为正式出图方式。

### dsub

```bash
dsub \
  -n hypo_DMR_heatmaps_5000bins \
  -R "cpu=10;mem=40G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_DMR_heatmaps_5000bins.out \
  -eo scheduler_logs/hypo_DMR_heatmaps_5000bins.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  04_plot_single_cell_dmr_heatmaps.py \
  --jobs 10 \
  --dmr-display-bins 5000 \
  --output-dir /share/home/rzli/METHSCAN/02_Methdiff/Heatmap/figures_diff0p30_hypo_DMRs_grouped_5000bins_blue_red
```

输出：

```text
/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/figures_diff0p30_hypo_DMRs_grouped_5000bins_blue_red
```

每个样本的主图：

```text
<sample>/<sample>__cells_by_all_DMRs_grouped_heatmap.png
```

其他文件：

- `<sample>__all_DMR_group_assignments.tsv.gz`：每个 DMR 的分组、各细胞类型平均值和原始行号；
- `<sample>__DMR_group_counts.tsv`：每个 DMR 组的原始 DMR 数和热图显示列数；
- `<sample>__DMR_heatmap_columns.tsv`：每个热图列与该列中原始 DMR 行范围的映射；
- `<sample>__heatmap_rows.tsv`：热图行、单细胞和细胞类型的映射；
- 默认不重复保存巨大的稠密绘图矩阵；只有显式加上 `--save-plot-matrix` 时才生成 `<sample>__all_DMRs_grouped_heatmap_matrix.npz`；
- `heatmap_summary.tsv`：每样本 DMR 总数、已分组数、`Unresolved` 数和热图列数。

## top1500 版本的第 02–04 步

第 05 步 `hypo_DMRs_diff0p25_top1500` 成功后，使用以下独立目录继续运行，不会覆盖全部 DMR 版本。

```bash
dsub \
  -n merge_hypo_diff0p25_top1500 \
  -R "cpu=10;mem=4G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/merge_hypo_diff0p25_top1500.out \
  -eo scheduler_logs/merge_hypo_diff0p25_top1500.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  02_merge_sample_dmrs.py \
  --jobs 10 \
  --input-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/celltype_hypo_DMRs_diff0p25_top1500 \
  --output-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p25_top1500
```

```bash
dsub \
  -n hypo_diff0p25_top1500_ratio \
  -R "cpu=16;mem=80G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_diff0p25_top1500_ratio.out \
  -eo scheduler_logs/hypo_diff0p25_top1500_ratio.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  03_compute_dmr_mean_cpg_ratio.py \
  --jobs 16 \
  --dmr-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p25_top1500 \
  --output-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p25_top1500
```

```bash
dsub \
  -n hypo_diff0p25_top1500_heatmaps \
  -R "cpu=10;mem=40G" \
  --cwd /share/home/rzli/METHSCAN/02_Methdiff/Heatmap \
  -oo scheduler_logs/hypo_diff0p25_top1500_heatmaps.out \
  -eo scheduler_logs/hypo_diff0p25_top1500_heatmaps.err \
  /share/home/rzli/miniconda3/envs/scDNAm/bin/python \
  04_plot_single_cell_dmr_heatmaps.py \
  --jobs 10 \
  --dmr-display-bins 5000 \
  --input-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p25_top1500 \
  --dmr-annotation-dir /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p25_top1500 \
  --output-dir /share/home/rzli/METHSCAN/02_Methdiff/Heatmap/figures_diff0p25_hypo_DMRs_top1500_grouped_5000bins_blue_red
```

## 运行检查

查看任务：

```bash
djob
```

查看日志：

```bash
tail -F scheduler_logs/<job>.out scheduler_logs/<job>.err
```

任务结束后，日志中应显示：

```text
EXIT_CODE: 0
Job execution succeeded
```

查看汇总：

```bash
column -t -s $'\t' \
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/celltype_hypo_DMRs_diff0p30/overall_summary.tsv

column -t -s $'\t' \
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30/merge_summary.tsv

column -t -s $'\t' \
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30/matrix_summary.tsv

column -t -s $'\t' \
/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/figures_diff0p30_hypo_DMRs_grouped_5000bins_blue_red/heatmap_summary.tsv
```

检查单细胞矩阵 gzip 文件：

```bash
for f in \
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30/*.tsv.gz
do
  gzip -t "$f" || echo "FAILED: $f"
done
```

## 重新运行与输出安全

四个脚本均不覆盖已有输出目录。如果目标目录已存在，脚本会停止并报错。重新运行时使用新目录，例如：

```bash
--output-dir /path/to/result_rerun
```

不要为了重跑而删除或清空原始 MethSCAn BED。
