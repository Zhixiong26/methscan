# IR/NR独立的VMR去个体效应前后聚类

## 完整流程

```text
All filtered data
├─ 01：IR/NR分别filter + smooth
├─ 02：IR/NR分别scan threshold 0.05/0.02/0.01
│      └─ 6套response-specific All VMR
├─ 03：旧第3/4/8步组内样本DMR
│      └─ IR和NR分别合并为individual-effect DMR union
├─ 04b：All VMR直接matrix
│      └─ 6套before matrix（去个体效应前）
├─ 04：All VMR − 同response的individual-effect DMR union
│      └─ 6套Clean VMR及after matrix（去个体效应后）
├─ 05：每套before/after分别PCA/UMAP/Leiden
│      └─ 12套独立聚类，不合并IR和NR，不叠加前后图片
└─ 06：汇总12套指标及6组after − before差值
```

## 脚本

| 脚本 | 输出 |
|---|---|
| `01_prepare_response_data.sh` | 2套response data |
| `02_scan_response_vmrs.sh` | 6套All VMR |
| `03_prepare_individual_effect_dmr_union.sh` | IR/NR individual-effect DMR union |
| `04b_build_response_before_vmr_matrix.sh` | 6套before All-VMR matrix |
| `04_build_response_clean_vmr_matrix.sh` | 6套after Clean-VMR matrix |
| `validate_response_matrix.py` | matrix细胞数、VMR数和特征集合验证 |
| `05_run_response_reclustering.sh` | 单套聚类入口 |
| `05_recluster_response_clean_vmrs.R` | PCA/UMAP/Leiden及指标 |
| `06_collect_response_metrics.sh` | 12套总表和前后差值表 |
| `06_compare_before_after_metrics.py` | 配对计算after − before |
| `07_prepare_filtered_celltype_dmrs.sh` | 筛选15套cell-type IR-vs-NR DMR并生成matrix regions |
| `07_filter_celltype_dmrs.py` | DMR格式检查、筛选及精确坐标去重 |
| `08_build_filtered_celltype_dmr_matrix.sh` | 在All filtered data上构建DMR matrix |
| `09_run_filtered_celltype_dmr_reclustering.sh` | 筛选DMR matrix的All-cell聚类入口 |
| `09_recluster_filtered_celltype_dmrs.R` | PCA/UMAP/Leiden、图和聚类指标 |
| `10_prepare_matched_individual_effect_clean_dmrs.sh` | 使用相同筛选条件构建同细胞类型individual-effect DMR union并从response DMR中相减 |

## 已完成统计

### Response data与All VMR

| Threshold | IR cells | IR All VMR | NR cells | NR All VMR |
|---|---:|---:|---:|---:|
| 0.05 | 8,644 | 162,322 | 7,597 | 158,850 |
| 0.02 | 8,644 | 85,636 | 7,597 | 80,679 |
| 0.01 | 8,644 | 47,659 | 7,597 | 43,184 |

### Individual-effect DMR union

| Response | Raw DMR | Merged DMR |
|---|---:|---:|
| IR | 60,922 | 57,511 |
| NR | 15,271 | 14,658 |
| Shared（只记录） | 76,193 | 67,167 |

### After Clean VMR matrix

| Threshold | Response | All VMR | Clean VMR | Removed | Removed % | Cells |
|---|---|---:|---:|---:|---:|---:|
| 0.05 | IR | 162,322 | 140,985 | 21,337 | 13.14% | 8,644 |
| 0.05 | NR | 158,850 | 151,121 | 7,729 | 4.87% | 7,597 |
| 0.02 | IR | 85,636 | 75,179 | 10,457 | 12.21% | 8,644 |
| 0.02 | NR | 80,679 | 76,179 | 4,500 | 5.58% | 7,597 |
| 0.01 | IR | 47,659 | 41,938 | 5,721 | 12.00% | 8,644 |
| 0.01 | NR | 43,184 | 40,316 | 2,868 | 6.64% | 7,597 |

聚类评价指标：

| 指标 | 简要含义 | 去个体效应后的理想方向 |
|---|---|---|
| Cell-type purity | Leiden cluster由同一种细胞类型构成的程度 | 保持不变或升高 |
| Sample purity | Leiden cluster被单个样本主导的程度 | 降低 |
| Sample mixing entropy | Leiden cluster内不同样本混合的均匀程度 | 升高 |

因此，`Sample purity下降 + Sample mixing entropy上升 + Cell-type purity保持或升高`，表示个体/样本效应减弱，同时细胞类型结构得到保留。

### Threshold 0.01聚类：Before与After

| Group | 指标 | Before | After | After − Before |
|---|---|---:|---:|---:|
| IR | Cell-type purity | 0.6913 | 0.7036 | +0.0122 |
| IR | Sample purity | 0.6821 | 0.4423 | -0.2398 |
| IR | Sample mixing entropy | 0.4759 | 0.8165 | +0.3406 |
| NR | Cell-type purity | 0.7124 | 0.7153 | +0.0029 |
| NR | Sample purity | 0.6066 | 0.4029 | -0.2037 |
| NR | Sample mixing entropy | 0.5748 | 0.8287 | +0.2539 |

统计结果：

- IR和NR的sample purity均明显下降，sample mixing entropy均明显上升。
- Cell-type purity未下降，IR和NR分别提高0.0122和0.0029。
- Before与After使用相同注释细胞：IR 8,248个，NR 7,252个。
- IR Leiden clusters由21变为24；NR由16变为25。
- ARI下降表示聚类结构发生重新组织，需结合上述sample和cell-type指标评价。

| Threshold | 去个体效果 | Cell-type保持 | 建议 |
|---|---|---|---|
| 0.01 | 最强 | 略微提高 | **主结果** |
| 0.02 | 中等 | 基本不变 | 稳健性结果 |
| 0.05 | 最弱 | 轻微下降 | 补充结果 |

判断依据是sample purity下降幅度、sample mixing entropy上升幅度，以及cell-type purity是否保持。三套threshold均改善样本混合，其中0.01改善最强且未损害细胞类型结构，因此作为主结果；0.02用于支持结果稳健性，0.05作为较弱阈值效果的补充对照。

### Threshold 0.02聚类：Before与After

| Group | 指标 | Before | After | After − Before |
|---|---|---:|---:|---:|
| IR | Cell-type purity | 0.6910 | 0.6887 | -0.0023 |
| IR | Sample purity | 0.7359 | 0.5603 | -0.1757 |
| IR | Sample mixing entropy | 0.4160 | 0.6904 | +0.2744 |
| NR | Cell-type purity | 0.7101 | 0.7125 | +0.0023 |
| NR | Sample purity | 0.6389 | 0.4691 | -0.1697 |
| NR | Sample mixing entropy | 0.5432 | 0.7530 | +0.2098 |

统计结果：

- IR和NR的sample purity分别下降0.1757和0.1697。
- IR和NR的sample mixing entropy分别提高0.2744和0.2098。
- Cell-type purity基本不变：IR下降0.0023，NR提高0.0023。
- IR Before和After各去除1个零方差VMR；无全NA VMR进入PCA。
- IR Leiden clusters由29变为31；NR由22变为28。

| Threshold | 去个体效果 | Cell-type保持 | 建议 |
|---|---|---|---|
| 0.01 | 最强 | 略微提高 | **主结果** |
| 0.02 | 中等 | 基本不变 | 稳健性结果 |
| 0.05 | 最弱 | 轻微下降 | 补充结果 |

### Threshold 0.05聚类：Before与After

| Group | 指标 | Before | After | After − Before |
|---|---|---:|---:|---:|
| IR | Cell-type purity | 0.6862 | 0.6817 | -0.0045 |
| IR | Sample purity | 0.7482 | 0.6140 | -0.1342 |
| IR | Sample mixing entropy | 0.3867 | 0.5964 | +0.2097 |
| NR | Cell-type purity | 0.7183 | 0.7169 | -0.0014 |
| NR | Sample purity | 0.6791 | 0.5863 | -0.0928 |
| NR | Sample mixing entropy | 0.4704 | 0.6180 | +0.1476 |

统计结果：

- IR和NR的sample purity分别下降0.1342和0.0928。
- IR和NR的sample mixing entropy分别提高0.2097和0.1476。
- Cell-type purity仅轻微下降：IR下降0.0045，NR下降0.0014。
- IR Before和After各去除2个零方差VMR；无全NA VMR进入PCA。
- IR Leiden clusters由40变为38；NR由42变为43。

| Threshold | 去个体效果 | Cell-type保持 | 建议 |
|---|---|---|---|
| 0.01 | 最强 | 略微提高 | **主结果** |
| 0.02 | 中等 | 基本不变 | 稳健性结果 |
| 0.05 | 最弱 | 轻微下降 | 补充结果 |

### 三种threshold比较

| Threshold | IR sample purity Δ | IR entropy Δ | NR sample purity Δ | NR entropy Δ | Cell-type purity |
|---|---:|---:|---:|---:|---|
| 0.01 | -0.2398 | +0.3406 | -0.2037 | +0.2539 | 两组均略升 |
| 0.02 | -0.1757 | +0.2744 | -0.1697 | +0.2098 | 基本不变 |
| 0.05 | -0.1342 | +0.2097 | -0.0928 | +0.1476 | 两组轻微下降 |

三套均改善样本混合，其中threshold 0.01改善最强且未损害cell-type purity，作为当前主结果；0.02作为稳健性结果，0.05改善最弱。

## 结果目录

```text
result/
├── response_VMR_before_individual/{IR,NR}/{threshold...}/
│   └── VMR_matrix/                         # before
├── response_VMR_remove_individual/{IR,NR}/{threshold...}/
│   └── VMR_matrix/                         # after，已完成
└── response_VMR_reclustering/{IR,NR}/{threshold...}/
    ├── before/
    └── after/
```

## 上传后检查

```bash
cd /share/home/rzli/METHSCAN/Meth_diff/20260718

bash -n \
  04b_build_response_before_vmr_matrix.sh \
  04_build_response_clean_vmr_matrix.sh \
  05_run_response_reclustering.sh \
  06_collect_response_metrics.sh

python -c \
  "import ast,pathlib; \
  [ast.parse(pathlib.Path(x).read_text()) for x in \
  ['validate_response_matrix.py','06_compare_before_after_metrics.py']]; \
  print('Python syntax OK')"

source /share/home/rzli/miniconda3/bin/activate scDNAm
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript -e \
  "parse(file='05_recluster_response_clean_vmrs.R'); cat('R syntax OK\n')"
```

## 04b：补建6套before matrix

总内存200 GB。每次并行IR和NR，各90 GB；threshold依次运行，先001。

```bash
v=threshold001

for g in IR NR; do
  dsub \
    -n "MethScan_before_matrix_${g}_${v}" \
    -R "cpu=32;mem=90G" \
    --cwd "$PWD" \
    -oo "logs/before_matrix_${g}_${v}.%J.out" \
    -eo "logs/before_matrix_${g}_${v}.%J.err" \
    env \
      GROUP="${g}" \
      VARIANT="${v}" \
      METHSCAN_THREADS=32 \
      METHSCAN_BLAS_THREADS=1 \
    bash 04b_build_response_before_vmr_matrix.sh
done
```

检查完成：

```bash
for g in IR NR; do
  dir="result/response_VMR_before_individual/${g}/${v}"
  test -f "${dir}/.complete" && echo "${g} COMPLETE" || echo "${g} NOT_COMPLETE"
  test -s "${dir}/matrix_validation.tsv" &&
    column -t -s $'\t' "${dir}/matrix_validation.tsv"
done
```

完成后依次改为 `v=threshold002`、`v=threshold005`，重复提交。

## 05：分别聚类before和after

同一 `threshold + stage` 的IR和NR可各用90 GB并行。建议先跑001的before，再跑001的after。

```bash
v=threshold001
stage=before

for g in IR NR; do
  dsub \
    -n "MethScan_cluster_${g}_${v}_${stage}" \
    -R "cpu=32;mem=90G" \
    --cwd "$PWD" \
    -oo "logs/cluster_${g}_${v}_${stage}.%J.out" \
    -eo "logs/cluster_${g}_${v}_${stage}.%J.err" \
    env \
      GROUP="${g}" \
      VARIANT="${v}" \
      STAGE="${stage}" \
      OPENBLAS_NUM_THREADS=32 \
      OMP_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 \
      METHSCAN_UMAP_THREADS=32 \
    bash 05_run_response_reclustering.sh
done
```

实时日志：

```bash
tail -F \
  "logs/response_reclustering/IR_${v}_${stage}.log" \
  "logs/response_reclustering/NR_${v}_${stage}.log"
```

`before`完成后设 `stage=after` 重复。再对002、005分别重复。

## 06：汇总12套结果

```bash
bash 06_collect_response_metrics.sh

column -t -s $'\t' \
  result/response_VMR_reclustering/before_after_metric_deltas.tsv
```

主要看：

```text
sample_cluster_purity_after_minus_before < 0
sample_mixing_entropy_after_minus_before > 0
cell_type_cluster_purity_after_minus_before 不应明显下降
```

## 新增支线：筛选cell-type IR-vs-NR DMR构建matrix

该支线不修改上述VMR流程。只读复用：

```text
/share/home/rzli/METHSCAN/Meth_diff/20260716/result/
└── DMR_results_200k/3_same_cell_type_IR_vs_NR/
    └── 15个 *_IR_vs_NR_DMRs.bed，共837,229行
```

筛选条件：

```text
raw p（第11列）< 0.05
abs(meth_frac_group1（第8列）− meth_frac_group2（第9列）) > 0.3
```

所有新结果写入当前工作目录：

```text
result/supervised_celltype_DMR_p005_absdiff030/
├── by_cell_type/                     # 每个细胞类型筛选后的原始12列DMR
├── selected_DMRs_with_source.tsv     # 全部入选DMR及来源
├── matrix_regions.bed                # methscan matrix输入
├── region_sources.tsv                # region与来源细胞类型对应关系
├── filter_summary.tsv
├── selection_metadata.tsv
├── DMR_matrix/                       # 第08步生成
└── reclustering/                     # 第09步生成
```

不同DMR不执行`bedtools merge`，不映射回VMR，也不改变重叠区域边界。`matrix_regions.bed`只对完全相同的`chr/start/end`去重，防止matrix产生重复特征列。

### 07：筛选DMR

```bash
dsub \
  -n MethScan-filter-celltype-DMR \
  -R "cpu=4;mem=8G" \
  --cwd "$PWD" \
  -oo logs/filter_celltype_DMR.%J.out \
  -eo logs/filter_celltype_DMR.%J.err \
  bash 07_prepare_filtered_celltype_dmrs.sh
```

完成后检查：

```bash
column -t -s $'\t' \
  result/supervised_celltype_DMR_p005_absdiff030/selection_metadata.tsv

column -t -s $'\t' \
  result/supervised_celltype_DMR_p005_absdiff030/filter_summary.tsv
```

### 08：All-cell DMR matrix

第07步统计出最终region数量后再确定内存并提交：

```bash
dsub \
  -n MethScan-filtered-DMR-matrix \
  -R "cpu=32;mem=MEMORY_TO_SET" \
  --cwd "$PWD" \
  -oo logs/filtered_DMR_matrix.%J.out \
  -eo logs/filtered_DMR_matrix.%J.err \
  env METHSCAN_THREADS=32 METHSCAN_BLAS_THREADS=1 \
  bash 08_build_filtered_celltype_dmr_matrix.sh
```

### 09：All-cell PCA/UMAP/Leiden

第08步完成后，根据DMR matrix规模确定内存并提交：

```bash
dsub \
  -n MethScan-filtered-DMR-cluster \
  -R "cpu=32;mem=MEMORY_TO_SET" \
  --cwd "$PWD" \
  -oo logs/filtered_DMR_cluster.%J.out \
  -eo logs/filtered_DMR_cluster.%J.err \
  env \
    OPENBLAS_NUM_THREADS=32 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    METHSCAN_UMAP_THREADS=32 \
  bash 09_run_filtered_celltype_dmr_reclustering.sh
```

该分析使用已知`cell_type + response`标签筛选DMR，再用这些DMR进行降维，因此属于监督式特征选择。UMAP上的细胞类型或IR/NR分离不能作为独立验证，只能作为所选DMR特征的描述性展示。

### 10：用相同筛选条件去除individual-effect DMR

对同细胞类型、同response内的样本两两DMR也使用第07步条件：第11列`p<0.05`且`|col8-col9|>0.3`。IR和NR分别合并后，构建每种细胞类型的individual-effect DMR union，再从该细胞类型response DMR中删除任意重叠区域。

```bash
dsub \
  -n MethScan-matched-individual-DMR \
  -R "cpu=4;mem=8G" \
  --cwd "$PWD" \
  -oo logs/matched_individual_DMR.%J.out \
  -eo logs/matched_individual_DMR.%J.err \
  bash 10_prepare_matched_individual_effect_clean_dmrs.sh
```

新结果不覆盖旧union：

```text
result/supervised_celltype_DMR_p005_absdiff030_remove_matched_individual/
├── individual_effect_union_by_cell_type/  # IR、NR及合并union
├── clean_response_by_cell_type/           # 相减后12列response DMR
├── filter_summary.tsv                     # 每种细胞类删除统计
├── selection_metadata.tsv                 # 总统计
└── matrix_regions.bed                     # 精确坐标去重，不合并部分重叠
```

先检查删除比例，再决定是否为新`matrix_regions.bed`重跑matrix：

```bash
column -t -s $'\t' \
  result/supervised_celltype_DMR_p005_absdiff030_remove_matched_individual/selection_metadata.tsv

column -t -s $'\t' \
  result/supervised_celltype_DMR_p005_absdiff030_remove_matched_individual/filter_summary.tsv
```
