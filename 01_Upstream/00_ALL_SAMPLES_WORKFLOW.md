# 10样本 MethSCAn 已验证单样本工作流

本目录只保留2026-08已在服务器运行成功的单样本分析线：10个样本分别完成cov去重、filter、smooth、样本内cell-type DMR、Top200、单细胞DMR matrix和统一热图。

本流程不合并10个样本的compact，不构建联合pseudobulk，不执行跨样本DMR检验。

## 脚本顺序

| 序号 | 脚本 | 作用 |
|---:|---|---|
| 01 | `01_check_cov_duplicates.sh` | 审计10个样本原始cov的重复CpG |
| 02 | `02_deduplicate_cov_by_probability.sh` | 按甲基化概率规则逐样本去重 |
| 03 | `03_run_upstream_pipeline.sh` | cov → compact/profile/filter/smooth通用实现 |
| 04 | `04_run_all_samples_to_smooth.sh` | 10样本独立运行到smooth |
| 05 | `05_run_all_samples_dmr.sh` | 每个样本内部进行cell-type两两DMR |
| 06 | `06_select_top200_dmrs.sh` | 逐样本筛选Top200 hypo-DMR并去除重叠区间 |
| 07 | `07_compute_top200_dmr_matrix.sh` | 逐样本计算single-cell × DMR mean-ratio matrix |
| 08 | `08_plot_all_top200_heatmaps.sh` | 统一生成/复用8类Top200热图 |
| 09 | `09_rerun_rawp_no_null_fdr.sh` | 对null-FDR但raw-p有效的比较执行回退补跑 |

05调用`lib/methdiff/run_single_sample_dmr.sh`；06–08使用`lib/methdiff/python/`中的`02、04、05、06` Python脚本。调度与计算实现已合并在同一流程目录。

## 固定分析口径

- Scanpy DMR注释：`/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_all_cells.csv`
- smooth前 clean-cell名单：`/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_clean_cells.csv`
- `Monocytes`是旧`CD14_Monocytes`与`CD16_Monocytes`的合并类型
- 排除`Platelet_erythroid_contamination`
- 样本：`25110891_IR01–IR05_Met`、`25110891_NR01–NR05_Met`
- 过滤：先应用`min_sites=300000`、`max_sites=1200000`、`min_meth=55`，再与Scanpy clean singlets取交集后运行smooth
- DMR：仅同样本内不同cell type两两比较，`min_cells=10`
- 染色体：仅`chr1–chr22, chrX, chrY`
- Top200：`raw p < 0.01`且`abs(methylation difference) >= 0.25`
- matrix：每个cell-DMR内唯一CpG ratio的等权算术平均
- 热图行：若某cell type在该样本最终Top200合并DMR中没有分配给自身的特异DMR，该type的全部细胞不纳入图中
- Z-score的均值与标准差只用筛选后保留的细胞计算

## 运行顺序

```bash
cd /share/home/rzli/METHSCAN/20260815/01_Upstream

bash 01_check_cov_duplicates.sh all 2 48
bash 02_deduplicate_cov_by_probability.sh all 2 48
bash 03_run_upstream_pipeline.sh run-to-compact 300k 10 1 all
bash 04_run_all_samples_to_smooth.sh run 10

bash 05_run_all_samples_dmr.sh prepare 2
bash 05_run_all_samples_dmr.sh run 2 2 24

SAMPLE_JOBS=10 bash 06_select_top200_dmrs.sh
SAMPLE_JOBS=1 CELL_JOBS=64 bash 07_compute_top200_dmr_matrix.sh

SAMPLE_JOBS=2 \
ZSCORE_MIN_OBSERVED_CELLS=30 \
ZSCORE_STANDARD_CLIP=3 \
bash 08_plot_all_top200_heatmaps.sh all
```

null-FDR回退状态与补跑：

```bash
bash 09_rerun_rawp_no_null_fdr.sh status
bash 09_rerun_rawp_no_null_fdr.sh all 2 16
```

## 8类正式热图

| `result/`目录 | 数值含义 | 颜色范围 |
|---|---|---|
| `MeanRatio` | 原始mean CpG ratio | `0–1` |
| `DMRwise_Zscore` | 逐DMR Z-score，数值截断至±3 | `-3–3` |
| `DMRwise_Zscore_Clip1` | 逐DMR Z-score，数值截断至±1 | `-1–1` |
| `DMRwise_Zscore_ColorClip1` | 未截断的逐DMR Z-score | 仅颜色在±1饱和 |
| `DMRwise_Zscore_MaxAbs` | `Z/max(abs(Z))` | `-1–1` |
| `DMRtypeMean_Zscore` | 先按DMR type等权平均，再Z-score，数值截断至±3 | `-3–3` |
| `DMRtypeMean_Zscore_Clip1` | DMR-type mean Z-score，数值截断至±1 | `-1–1` |
| `DMRtypeMean_Zscore_MaxAbs` | DMR-type mean Z-score再max-abs归一化 | `-1–1` |

`ColorClip1`不改变matrix：`Z≤-1`都显示为最深蓝，`Z≥1`都显示为最深红，保存的Z-score数值仍保持原值。

状态与软链接：

```bash
bash 08_plot_all_top200_heatmaps.sh status
bash 08_plot_all_top200_heatmaps.sh links
```

`status`审计`10样本 × 8图型 = 80`项；`links`只在80张PNG都完整时替换`result/`。
