# Shared individual-effect mask去除个体效应后的threshold VMR聚类

服务器目录：

```bash
cd /share/home/rzli/METHSCAN/Meth_diff/20260716
```

## 完整流程

```text
MethSCAn filtered data（16,241个细胞）
└─ ALL_annotation_200k.csv（匹配15,500个细胞）
   ↓
03/04：同细胞类型、同response组内的样本两两DMR
   ↓
08：按cell type × response筛选并合并q < 0.05 DMR
├─ IR：14个union BED
└─ NR：14个union BED
   ↓
10：合并IR和NR的组内样本DMR
└─ shared individual-effect mask（67,167个区域）
   ↓
11：分别scan threshold 0.05/0.02/0.01
├─ All VMR − shared individual-effect mask
├─ 3套Clean VMR
└─ 3套Clean VMR matrix（每套16,241个细胞）
   ↓
13：每套matrix分别PCA/UMAP/Leiden
└─ 匹配15,500个Scanpy注释细胞
   ↓
14：汇总3套聚类指标并比较threshold
```

`individual-effect mask`由组内样本DMR组成。任何与mask区间重叠的VMR都会被整体排除：

```text
Clean VMR = All VMR − 与individual_effect_union_q005.bed任意重叠的VMR
```

`01/02`、`05/06`是补充DMR分析；`07/09/12`属于旧流程，不进入上述主链。

## 脚本

| 脚本 | 输出 |
|---|---|
| `03_meth_diff_celltype_sample_pairwise_200k.sh` | 同细胞类型、同response组内的样本两两DMR |
| `04_generate_celltype_sample_pairwise_groups.R` | cell type × response × sample pair分组 |
| `08_merge_celltype_sample_pairwise_dmrs.sh` | IR和NR各14个q < 0.05 union BED |
| `10_prepare_individual_effect_mask.sh` | shared individual-effect DMR mask |
| `11_run_threshold_vmrs_remove_individual.sh` | 3套Clean VMR及matrix |
| `validate_threshold_matrix.py` | matrix细胞数、VMR数和特征集合验证 |
| `12_build_threshold_before_vmr_matrix.sh` | 3套去除前All VMR matrix |
| `13_run_threshold_clean_vmr_reclustering.sh` | 单套聚类入口 |
| `13_run_threshold_before_vmr_reclustering.sh` | 单套去除前聚类入口 |
| `13_recluster_threshold_clean_vmrs.R` | PCA/UMAP/Leiden及指标 |
| `14_collect_threshold_metrics.sh` | 汇总3套comparison metrics |

补充和旧流程脚本已移至`archive_old_workflow/`：

| 脚本 | 用途 |
|---|---|
| `archive_old_workflow/01_meth_diff_pairwise_200k.sh`、`02_generate_cell_groups.R` | 细胞类型内IR vs NR DMR |
| `archive_old_workflow/05_meth_diff_sample_pairwise_200k.sh`、`06_generate_sample_pairwise_groups.R` | 不区分细胞类型的组内样本DMR |
| `archive_old_workflow/07_merge_celltype_ir_vs_nr_dmrs.sh` | 合并细胞类型内IR vs NR DMR |
| `archive_old_workflow/09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh` | 旧DMR相减流程 |
| `archive_old_workflow/10_map_clean_dmrs_to_all_vmrs.sh`、`11_select_top_clean_vmrs.sh` | 旧Clean DMR映射和筛选流程 |
| `archive_old_workflow/legacy/12_*`、`13_run_top_clean_vmr_reclustering.sh` | 旧二次截取和聚类流程 |

## 已完成统计

### 输入与细胞

| 项目 | 统计 |
|---|---:|
| MethSCAn filtered data细胞 | 16,241 |
| 匹配Scanpy注释 | 15,500 |
| 未匹配 | 741 |
| 正式DMR有效细胞 | 15,466 |
| IR有效细胞 | 8,227 |
| NR有效细胞 | 7,239 |
| 排除污染细胞 | 34 |
| `MIN_CELLS` | 6 |

| 样本 | 细胞数 | 样本 | 细胞数 |
|---|---:|---|---:|
| IR01 | 1,523 | NR01 | 917 |
| IR02 | 2,012 | NR02 | 1,981 |
| IR03 | 2,173 | NR03 | 1,556 |
| IR04 | 1,132 | NR04 | 1,785 |
| IR05 | 1,387 | NR05 | 1,000 |

### Individual-effect DMR mask

| Response | Union BED | 非空BED |
|---|---:|---:|
| IR | 14 | 14 |
| NR | 14 | 12 |
| 合计 | 28 | 26 |

```text
IR 14个union BED + NR 14个union BED
→ 排除2个空文件
→ 从26个非空BED提取chr/start/end
→ 拼接为76,193个区域
→ 按chr/start/end排序并bedtools merge
→ shared individual-effect mask：67,167个非重叠区域
```

第10步输出：

```text
result/individual_effect_mask/individual_effect_union_q005.bed
result/individual_effect_mask/source_files.tsv
result/individual_effect_mask/mask_summary.tsv
```

### Clean VMR matrix

三套threshold均使用同一个shared individual-effect mask。

| Variant | Threshold | All VMR | Clean VMR | Removed | Removed % | Cells |
|---|---:|---:|---:|---:|---:|---:|
| `threshold005` | 0.05 | 159,457 | 133,130 | 26,327 | 16.51% | 16,241 |
| `threshold002` | 0.02 | 83,245 | 70,110 | 13,135 | 15.78% | 16,241 |
| `threshold001` | 0.01 | 45,119 | 37,867 | 7,252 | 16.07% | 16,241 |

矩阵验证：

| Variant | Matrix VMR | Missing | Extra | BED顺序一致 | 状态 |
|---|---:|---:|---:|---|---|
| `threshold005` | 133,130 | 0 | 0 | 否 | 通过 |
| `threshold002` | 70,110 | 0 | 0 | 否 | 通过 |
| `threshold001` | 37,867 | 0 | 0 | 否 | 通过 |

`BED顺序一致=否`仅表示MethSCAn重新排列了矩阵列，VMR坐标集合完全一致。

### 聚类评价指标

| 指标 | 简要含义 | 去个体效应后的理想方向 |
|---|---|---|
| Cell-type purity | Leiden cluster由同一种细胞类型构成的程度 | 保持不变或升高 |
| Sample purity | Leiden cluster被单个样本主导的程度 | 降低 |
| Sample mixing entropy | Leiden cluster内不同样本混合的均匀程度 | 升高 |

因此，`Sample purity下降 + Sample mixing entropy上升 + Cell-type purity基本保持`表示样本效应减弱，同时细胞类型结构得到保留。

统一聚类参数：

```text
矩阵细胞：16,241
匹配Scanpy注释：15,500
PCA：20 PCs
UMAP：30 neighbors，min_dist=0.05，seed=2
Leiden：resolution=0.001
```

原始200k参考：

| Cell-type purity | Sample purity | Sample mixing entropy |
|---:|---:|---:|
| 0.7203 | 0.5154 | 0.5688 |

Threshold结果：

| Variant | Cell-type purity | Sample purity | Sample mixing entropy | Response purity | Response entropy | ARI vs 200k |
|---|---:|---:|---:|---:|---:|---:|
| `threshold001` | 0.7126 | 0.3005 | 0.8443 | 0.6094 | 0.9167 | 0.4033 |
| `threshold002` | 0.7135 | 0.3369 | 0.8090 | 0.6312 | 0.8910 | 0.4932 |
| `threshold005` | 待完成 | 待完成 | 待完成 | 待完成 | 待完成 | 待完成 |

相对原始200k参考的变化：

| Variant | Cell-type purity Δ | Sample purity Δ | Sample entropy Δ |
|---|---:|---:|---:|
| `threshold001` | -0.0077 | -0.2149 | +0.2754 |
| `threshold002` | -0.0068 | -0.1785 | +0.2402 |
| `threshold005` | 待完成 | 待完成 | 待完成 |

当前判断：

- `threshold001`的样本效应去除最强。
- `threshold002`的聚类结构更接近原始200k结果。
- `threshold005`等待Job 162819完成后再比较。

### 运行记录

| Variant | Job ID | OpenBLAS线程 | 运行时间 | VMR | Leiden簇 | 状态 |
|---|---:|---:|---:|---:|---:|---|
| `threshold001` | 162815 | 32 | 41分25秒 | 37,867 | 59 | COMPLETE |
| `threshold002` | 162818 | 64 | 1时14分45秒 | 70,110 | 42 | COMPLETE |
| `threshold005` | 162819 | 32 | 待完成 | 133,130 | 待完成 | RUNNING |

| Variant | All-NA VMR | Zero-variance VMR | 缺失注释细胞 |
|---|---:|---:|---:|
| `threshold001` | 0 | 0 | 0 |
| `threshold002` | 0 | 0 | 0 |
| `threshold005` | 待完成 | 待完成 | 待完成 |

补充DMR结果：

- 细胞类型内IR vs NR：15个比较，Job 162782。
- 不区分细胞类型的组内样本比较：IR 10组、NR 10组，Job 162785。
- 全细胞样本比较得到IR DMR 801,366个、NR DMR 828,608个；该结果不进入当前individual-effect mask。

## 结果目录

```text
result/
├── celltype_sample_pairwise/merged_DMRs_200k/q005/
│   ├── IR/                                  # 14个IR union BED
│   └── NR/                                  # 14个NR union BED
├── individual_effect_mask/
│   ├── individual_effect_union_q005.bed     # shared mask
│   ├── source_files.tsv
│   └── mask_summary.tsv
├── threshold_VMR_remove_individual/
│   └── {threshold001,threshold002,threshold005}/
│       ├── all_VMRs.bed
│       ├── clean_VMRs.bed
│       ├── removed_individual_effect_VMRs.bed
│       ├── VMR_matrix/mean_shrunken_residuals.csv.gz
│       ├── matrix_validation.tsv
│       ├── run_metadata.tsv
│       └── .complete
├── threshold_VMR_before_individual/
│   └── {threshold001,threshold002,threshold005}/
│       ├── all_VMRs.bed
│       ├── VMR_matrix/mean_shrunken_residuals.csv.gz
│       ├── matrix_validation.tsv
│       ├── run_metadata.tsv
│       └── .complete
├── threshold_VMR_before_individual_reclustering/
│   └── {threshold001,threshold002,threshold005}/
│       ├── comparison_metrics.tsv
│       ├── plots/
│       └── .complete
└── threshold_VMR_remove_individual_reclustering/
    └── {threshold001,threshold002,threshold005}/
        ├── feature_qc.tsv
        ├── feature_qc_summary.tsv
        ├── comparison_metrics.tsv
        ├── plots/
        └── .complete
```

旧的threshold005单线程部分结果：

```text
result/threshold_VMR_remove_individual_reclustering/
└── threshold005.singlethread_partial_162814/
```

## 上传后检查

```bash
cd /share/home/rzli/METHSCAN/Meth_diff/20260716

bash -n \
  03_meth_diff_celltype_sample_pairwise_200k.sh \
  08_merge_celltype_sample_pairwise_dmrs.sh \
  10_prepare_individual_effect_mask.sh \
  11_run_threshold_vmrs_remove_individual.sh \
  13_run_threshold_clean_vmr_reclustering.sh \
  14_collect_threshold_metrics.sh

python -c \
  "import ast,pathlib; \
  ast.parse(pathlib.Path('validate_threshold_matrix.py').read_text()); \
  print('Python syntax OK')"

source /share/home/rzli/miniconda3/bin/activate scDNAm
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript -e \
  "parse(file='13_recluster_threshold_clean_vmrs.R'); cat('R syntax OK\n')"
```

## 10–11：构建shared mask和3套Clean VMR matrix

```bash
bash 10_prepare_individual_effect_mask.sh
bash 11_run_threshold_vmrs_remove_individual.sh
```

检查每套matrix：

```bash
for v in threshold001 threshold002 threshold005; do
  dir="result/threshold_VMR_remove_individual/${v}"
  test -f "${dir}/.complete" &&
    echo "${v} COMPLETE" ||
    echo "${v} NOT_COMPLETE"
  test -s "${dir}/matrix_validation.tsv" &&
    column -t -s $'\t' "${dir}/matrix_validation.tsv"
done
```

## 13：分别聚类3套Clean VMR matrix

```bash
v=threshold001

dsub \
  -n "MethScan_cluster_${v}" \
  -R "cpu=32;mem=90G" \
  --cwd "$PWD" \
  -oo "logs/cluster_${v}.%J.out" \
  -eo "logs/cluster_${v}.%J.err" \
  env \
    VARIANT="${v}" \
    OPENBLAS_NUM_THREADS=32 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
  bash 13_run_threshold_clean_vmr_reclustering.sh
```

完成后依次将`v`改为`threshold002`和`threshold005`。

## 14：汇总3套结果

3套聚类全部完成后运行：

```bash
bash 14_collect_threshold_metrics.sh

column -t -s $'\t' \
  result/threshold_VMR_remove_individual_reclustering/comparison_metrics_all.tsv
```

主要关注：

```text
sample purity降低
sample mixing entropy升高
cell-type purity不应明显下降
```
