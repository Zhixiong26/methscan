# 10样本独立 cov → DMR → matrix 工作流

当前流程不再合并10个样本的compact或染色体矩阵。每个样本独立完成QC、filter、smooth和样本内细胞类型DMR；样本之间只在调度层面并行。

## 脚本顺序

| 序号 | 脚本 | 作用 |
|---:|---|---|
| 01 | `01_check_cov_duplicates.sh` | 审计10个样本原始cov中的重复CpG |
| 02 | `02_deduplicate_cov_by_probability.sh` | 按甲基化概率规则对每个样本单独去重 |
| 03 | `03_run_upstream_pipeline.sh` | 通用单/多样本上游实现；cov → compact/profile/filter/smooth |
| 04 | `04_run_all_samples_to_smooth.sh` | 调用03，让10个样本独立运行到smooth |
| 05 | `05_run_all_samples_dmr.sh` | 每个样本内部进行细胞类型两两DMR |
| 06 | `06_select_top200_dmrs.sh` | 每个样本独立筛选Top200 hypo-DMR，并仅在该样本内去重重叠区间 |
| 07 | `07_compute_top200_dmr_matrix.sh` | 每个样本用自身去重cov计算单细胞 × DMR mean-ratio矩阵 |
| 08 | `08_plot_top200_dmr_heatmap.sh` | 可选：每个样本独立绘制mean-ratio热图 |
| 09 | `09_plot_top200_dmr_zscore_heatmap.sh` | 可选：基于07的同一矩阵绘制逐DMR列Z-score热图 |

`02_Methdiff/run_single_sample_dmr.sh` 是05调用的通用单样本DMR实现。`02_Methdiff/Result/`中的Python脚本是06–09的计算实现。

## 固定分析口径

- 样本：`25110891_IR01–IR05_Met`、`25110891_NR01–NR05_Met`
- 原始cov：`<sample>/cov`
- 去重cov：`<sample>/cov_dedup_probability`
- compact：`<sample>/compact_data_dedup_probability`
- QC根目录：`<sample>/qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob`
- 过滤：`min_sites=300000`、`max_sites=10000000`、`min_meth=55`
- filtered：`<QC根目录>/filtered_data_single_300k`
- smooth：`<filtered>/smoothed`
- DMR根目录：`<QC根目录>/methdiff_celltype_300k`
- DMR比较：仅同一样本内不同细胞类型两两比较，`min_cells=10`
- DMR染色体：仅`chr1–chr22, chrX, chrY`
- Top200：`raw p < 0.01`且`abs(methylation difference) >= 0.25`
- matrix：每个cell-DMR内唯一CpG ratio的等权算术平均
- mean-ratio热图（08）：直接显示07的matrix值，范围`0–1`，`NA`为灰色
- Z-score热图（09）：保留07的完全相同cell顺序、DMR顺序及所有DMR；对每个DMR列仅用非`NA`细胞计算`Z=(ratio-mean_DMR)/SD_DMR`
- Z-score资格：默认每个DMR至少30个观测细胞且`SD_DMR > 0`；不合格DMR整列显示为灰色`NA`
- Z-score显示：默认截断至`[-3, +3]`，蓝=低于该DMR的细胞均值、白=该DMR的细胞均值、红=高于该DMR的细胞均值
- 08与09：均保留左侧单细胞cell type色带、顶部DMR type色带；09的Z-score是matrix的可视化标准化，并非MethSCAn DMR统计检验的Z-score

不再使用：

```text
/share/LCZX_Data/data/allcools/merged_10samples_covdedupprob
```

## 运行顺序

### 1. 重复审计

```bash
bash 01_check_cov_duplicates.sh all 2 48
```

### 2. 概率去重

```bash
bash 02_deduplicate_cov_by_probability.sh all 2 48
```

### 3. 生成每个样本compact

```bash
bash 03_run_upstream_pipeline.sh run-to-compact 300k 10 1 all
```

### 4. 每个样本独立filter和smooth

已有合格compact会复用，不会重新prepare：

```bash
bash 04_run_all_samples_to_smooth.sh run 10
```

状态：

```bash
bash 04_run_all_samples_to_smooth.sh status
```

### 5. 每个样本独立DMR

先准备注释和比较组：

```bash
bash 05_run_all_samples_dmr.sh prepare 10
```

再计算DMR。下例最多同时使用`2 × 2 × 24 = 96`个MethSCAn线程：

```bash
bash 05_run_all_samples_dmr.sh run 2 2 24
```

状态：

```bash
bash 05_run_all_samples_dmr.sh status
```

### 6. 每个样本独立Top200

```bash
SAMPLE_JOBS=10 bash 06_select_top200_dmrs.sh
```

### 7. 每个样本独立matrix

默认逐样本运行，每个样本使用64个cell worker：

```bash
SAMPLE_JOBS=1 CELL_JOBS=64 bash 07_compute_top200_dmr_matrix.sh
```

### 8. 可选热图

```bash
PLOT_JOBS=2 bash 08_plot_top200_dmr_heatmap.sh
```

### 9. 可选逐DMR列Z-score热图

输入与08完全一致，仅转换显示数值；不会改写matrix或mean-ratio热图。

```bash
PLOT_JOBS=2 \
ZSCORE_MIN_OBSERVED_CELLS=30 \
ZSCORE_CLIP=3 \
bash 09_plot_top200_dmr_zscore_heatmap.sh
```

输出：

```text
<sample>/<QC根目录>/methdiff_celltype_300k/
heatmap_top200_rawp0p01_diff0p25/figures_top200_DMRwise_zscore/
```

## 关键解释

- 04中的10路并发是10个样本独立运行，不会构建联合pseudobulk。
- 每个样本的smooth背景只由该样本通过QC的细胞产生。
- 05不会进行IR与NR之间、或不同样本之间的DMR比较。
- 06只在同一个样本内对Top200 DMR重叠区间去重，不会合并样本数据。
- 07直接读取每个样本自己的`cov_dedup_probability`。
- 09的Z-score按DMR列在单个样本内部标准化；它不能用于比较不同样本的绝对甲基化水平。
