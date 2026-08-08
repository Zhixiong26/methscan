# 归档脚本

本目录保存当前主流程不再直接使用的补充分析、旧流程和兼容入口。

## 补充DMR分析

| 脚本 | 用途 |
|---|---|
| `01_meth_diff_pairwise_200k.sh` | 细胞类型内IR vs NR DMR |
| `02_generate_cell_groups.R` | 为01生成cell-group |
| `05_meth_diff_sample_pairwise_200k.sh` | 不区分细胞类型的组内样本DMR |
| `06_generate_sample_pairwise_groups.R` | 为05生成sample pair分组 |

## 旧流程

| 脚本 | 用途 |
|---|---|
| `07_merge_celltype_ir_vs_nr_dmrs.sh` | 合并细胞类型内IR vs NR DMR |
| `09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh` | 旧DMR相减流程 |
| `10_map_clean_dmrs_to_all_vmrs.sh` | 将旧Clean DMR映射到All VMR |
| `11_select_top_clean_vmrs.sh` | 兼容入口，转调主目录中的threshold流程 |
| `13_run_top_clean_vmr_reclustering.sh` | 兼容入口，转调主目录中的threshold聚类 |
| `run_steps_07_09.sh` | 旧07–09串联入口 |
| `run_steps_10_11.sh` | 旧10–11包装入口 |
| `run_steps_12_13.sh` | 旧12–13包装入口 |
| `legacy/` | 二次截取Top Clean VMR的更早版本 |

这些文件仅用于追溯既往分析，不属于README所述的当前主链。
