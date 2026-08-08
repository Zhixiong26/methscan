# Full merged VMR PCA–UMAP–Leiden

本目录直接对 `01_Upstream/run_merged_pipeline.sh` 生成的完整 merged-30k VMR 矩阵进行 PCA、UMAP 和 Leiden 聚类。

主入口：

```text
run_vmr_clustering.sh
        ↓
vmr_clustering.R
```

## 1. 分析口径

默认不对细胞或 VMR 下采样：

```text
52,561 个 filtered cells
×
88,261 个 VMRs
```

包括：

- 全部52,561个细胞，包括1,537个没有 Scanpy 注释的细胞；
- `chr1–chr22`、`chrX`、`chrY`、`chrM`、`GL*`、`KI*` 等矩阵中的全部VMR；
- 不设置 `MAX_CELLS`；
- 不设置 `MAX_FEATURES`；
- 不做覆盖率或方差预筛选。

Scanpy 注释仅在 PCA/UMAP/Leiden 完成后用于 cell type、sample 和 IR/NR 着色及审计，不参与特征选择或聚类。没有 Scanpy 注释的细胞显示为 `Unannotated`。

## 2. 与旧 Annotation 脚本的关系

方法参考：

```text
Annotation/20260716/scripts/02_All_200k_analysis.R
```

保留的核心实现：

1. `data.table::fread` 一次读取完整 `mean_shrunken_residuals.csv.gz`；
2. 转换为 cell × VMR dense numeric matrix；
3. `scale(center=TRUE, scale=FALSE)` 按VMR中心化，不做方差标准化；
4. 缺失位置首次填0；
5. `irlba::prcomp_irlba` 拟合PCA；
6. 使用 `pca$x %*% t(pca$rotation)` 重构缺失位置；
7. 最多迭代50次，或 relative gain `<0.001` 时提前结束；
8. 使用PC1–PC20做 UMAP；
9. 使用 UMAP 返回的 kNN 距离构建加权图；
10. 对该图执行 Leiden CPM 聚类。

新脚本与旧脚本的主要数据差异是：输入改为当前 merged-30k 的52,561个细胞和88,261个VMR，并且保留全部无注释细胞。

## 3. 默认输入

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/
├── metadata/sample_batch.tsv
└── qc_minmeth55_maxmethnone_maxsites10000000/
    ├── filtered_data_merged_30k/column_header.txt
    └── VMR_matrix_merged_30k/
        └── mean_shrunken_residuals.csv.gz
```

Scanpy 注释：

```text
/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv
```

脚本通过 `metadata/sample_batch.tsv` 将矩阵细胞ID转换为 `sample + barcode` 键，再与 Scanpy 注释对齐。

## 4. 默认参数

```text
EXPECTED_CELLS=52561
EXPECTED_VMRS=88261
N_PCS=20
PCA_ITERATIONS=50
PCA_MIN_GAIN=0.001
UMAP_N_NEIGHBORS=30
UMAP_MIN_DIST=0.05
LEIDEN_RESOLUTION=0.001
RANDOM_SEED=2
THREADS=32
UMAP_THREADS=32
```

与旧 Annotation 的主要分析参数一致：

```text
PCA: 20 PCs
UMAP: n_neighbors=30, min_dist=0.05, seed=2
Leiden: CPM resolution=0.001
```

如果使用 10k/20k/50k 矩阵，其细胞数或VMR数与30k不同，应设置：

```bash
THRESHOLD=20k \
EXPECTED_CELLS=auto \
EXPECTED_VMRS=auto \
ANALYSIS_LABEL=fullmatrix_20k \
bash run_vmr_clustering.sh run
```

## 5. 八个运行阶段

```text
[1/8] CHECK       完整矩阵、metadata、annotation和R依赖
[2/8] LOAD        一次读取52,561 × 88,261完整矩阵
[3/8] PCA         中心化、缺失值迭代重构和PCA20
[4/8] PCA-OUTPUT  PCA模型、坐标、方差和每轮MSE
[5/8] UMAP        PC1–PC20上的UMAP
[6/8] LEIDEN      UMAP kNN加权图上的Leiden CPM
[7/8] REPORT      PCA/UMAP图、聚类交叉表和ARI/NMI
[8/8] COMPLETE    完成标记
```

PCA 和 UMAP/Leiden 分为两个独立 R 进程。PCA 结束后大矩阵占用会随 R 进程退出而释放，之后 `cluster` 只读取PC坐标。

## 6. 上传后检查

服务器目录：

```text
/share/home/rzli/METHSCAN/04_VMR_clustering
```

检查：

```bash
cd /share/home/rzli/METHSCAN/04_VMR_clustering

chmod 750 run_vmr_clustering.sh vmr_clustering.R
bash -n run_vmr_clustering.sh

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

Rscript -e 'invisible(parse(file="vmr_clustering.R")); cat("R syntax OK\n")'
Rscript -e 'library(data.table); library(irlba); library(uwot); library(igraph); library(ggplot2); cat("R packages OK\n")'

bash run_vmr_clustering.sh --help
bash run_vmr_clustering.sh status
```

## 7. dsub 提交

不要在 node-4 登录节点直接执行 `pca/cluster/run`。

完整提交：

```bash
cd /share/home/rzli/METHSCAN/04_VMR_clustering
mkdir -p scheduler_logs

dsub \
  -n vmr_full_30k_pca_umap \
  -R "cpu=32;mem=194560MB" \
  --cwd /share/home/rzli/METHSCAN/04_VMR_clustering \
  -oo scheduler_logs/vmr_full_30k_pca_umap.out \
  -eo scheduler_logs/vmr_full_30k_pca_umap.err \
  bash run_vmr_clustering.sh run
```

这里 `194560MB = 190 GiB`，低于200G节点上限。

监控：

```bash
djob

tail -F \
  scheduler_logs/vmr_full_30k_pca_umap.out \
  scheduler_logs/vmr_full_30k_pca_umap.err
```

如果 PCA 已经成功，但 UMAP/Leiden 需要重跑，可单独提交：

```bash
dsub \
  -n vmr_full_30k_cluster \
  -R "cpu=32;mem=32768MB" \
  --cwd /share/home/rzli/METHSCAN/04_VMR_clustering \
  -oo scheduler_logs/vmr_full_30k_cluster.out \
  -eo scheduler_logs/vmr_full_30k_cluster.err \
  bash run_vmr_clustering.sh cluster
```

## 8. 状态和断点

```bash
bash run_vmr_clustering.sh status
```

状态表输出：

```text
cells  VMRs  completed_iterations  leiden_clusters  pca  clustered
```

支持的动作：

```text
pca      读取完整矩阵并执行迭代PCA
cluster  复用已有PCA，执行UMAP/Leiden和报告
run      先PCA，再UMAP/Leiden
status   只读状态，不启动计算
```

成功标记：

```text
.pca.ok
.cluster.ok
```

如果改变PC数、迭代数、UMAP或Leiden参数，必须使用新 `ANALYSIS_LABEL`，避免复用上一套成功标记。

## 9. 输出目录

默认：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/
└── vmr_clustering_30k/
    └── fullmatrix_allcells_iter50_pcs20_nn30_md0p05_lei0p001_seed2/
```

主要输出：

```text
cell_metadata.tsv.gz            52,561个细胞的sample/response/cell type和注释状态
matrix_features.tsv.gz          88,261个VMR及原矩阵顺序
iterative_pca_model.rds         irlba PCA模型、scores、loadings和迭代MSE
pca_coordinates.tsv.gz          PC1–PC20
pca_explained_variance.tsv      PCA方差信息
pca_imputation_mse.tsv          每轮缺失值重构MSE
pca_summary.tsv                 PCA参数、完成迭代数和软件版本
umap_coordinates.tsv.gz         UMAP1和UMAP2
cell_embeddings.tsv.gz          metadata + PC1–PC20 + UMAP + Leiden
leiden_by_response.tsv
leiden_by_sample.tsv
leiden_by_cell_type.tsv
clustering_summary.tsv          簇大小、ARI/NMI、参数和软件版本
plots/
```

`plots/` 包含：

```text
pca_scree.png
pca_by_response.png
pca_by_sample.png
pca_by_cell_type.png
pca_by_leiden.png
umap_by_response.png
umap_by_sample.png
umap_by_cell_type.png
umap_by_leiden.png
```

## 10. 结果解释边界

- VMR 是在这批52,561个细胞上由 `methscan scan` 发现的，因此属于当前数据驱动特征，不是独立验证集。
- PCA、UMAP 和 Leiden 不接收 cell type、sample 或 IR/NR 标签。
- UMAP 只用于低维可视化，不应在UMAP坐标上做正式差异检验。
- Leiden 簇不自动等于细胞类型，需要结合 Scanpy 标签和每个样本的细胞数解释。
- 默认不进行 batch correction。必须同时检查 `umap_by_sample.png`、`leiden_by_sample.tsv` 和 sample ARI，不能将 sample 效应直接解释为 IR/NR 差异。
