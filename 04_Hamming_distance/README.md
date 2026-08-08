# scWGBS DMR-Hamming distance clustering

本目录实现以 `02_Methdiff` 的“同一 cell type 内 IR vs NR”DMR为特征区域，在单细胞 CpG 调用层面计算 Hamming distance，并按 cell type 分别完成层次聚类及 `MDS → UMAP → Leiden` 分析。

本流程不包含博士论文 2.19.3 的全细胞类型扩展注释、双组学标签转移、mcall 合并或 IGV 验证。

入口：

```text
run_hamming_pipeline.sh
hamming_scwgbs.py
```

## 1. 方法定义

对一个 cell type（例如 B cells）：

1. 从 `B_cells__IR_vs_NR_DMRs.bed` 筛选 DMR；
2. 从 MethSCAn `filtered_data_merged_30k/chr*.npz` 直接提取这些 DMR 内的 CpG 位点；
3. MethSCAn 稀疏矩阵中 `+1` 转为甲基化 `1`，`-1` 转为未甲基化 `0`，稀疏零值作为未观测；
4. 两个细胞只在双方均观测到的 CpG 上计算：

```text
Hamming(i,j) = 不同甲基化调用数 / 共同观测CpG数
```

5. 没有足够共同 CpG 的细胞对不进行距离填补。脚本贪心移除造成最多无效配对的低覆盖细胞，直到保留细胞间的所有距离均可计算；
6. 默认使用 `average linkage` 进行层次聚类；
7. 使用 metric MDS 将预计算 Hamming distance 降至10维；
8. 在 MDS 坐标上构建10邻居 UMAP fuzzy graph，生成二维 UMAP；
9. 在同一邻接图上运行 Leiden（resolution=0.5）；
10. 输出 dendrogram、距离热图、MDS坐标、UMAP、Leiden assignment，以及 IR/NR 和 sample 审计表。

MethSCAn 官方说明其染色体矩阵是 `genomic position × cell` CSR 矩阵，甲基化为 `+1`、未甲基化为 `-1`、缺失为稀疏 `0`：

<https://anders-biostat.github.io/MethSCAn/commands.html>

因此，本流程不需要把原始单细胞文件逐个重新 `bedtools intersect`，也不把缺失值误当作未甲基化。

本流程与 Hammeth `matrix` 步骤的二值化、缺失值和 pairwise Hamming 定义一致，但跳过 BAM→PAT 和基于 read-pattern Hamming 重新发现高变 bins；特征区域直接使用现有 MethSCAn response DMR。因此应描述为“参考 Hammeth 矩阵定义的 DMR-restricted Hamming 分析”，而不是完整 Hammeth `prepare → hamdist → matrix` 复现。

Hammeth 资料：

- <https://github.com/2762038415/HamDis/blob/main/HAMMETH_TUTORIAL_GITHUB_CN.md>
- <https://github.com/2762038415/HamDis/blob/main/hammeth/hammeth/scripts/matrix.py>

## 2. 为什么不默认使用 Ward

截图中的流程写的是“Hamming distance + Ward”。Ward linkage 的目标函数建立在欧氏平方距离上；把任意的 pairwise-complete Hamming distance 直接传给 Ward，不能保持 Ward 方差最小化的统计含义。

本流程默认：

```text
LINKAGE_METHOD=average
```

可选 `complete` 或 `single`。为了复现截图，脚本保留显式 Ward 敏感性分析入口，但必须同时设置：

```bash
LINKAGE_METHOD=ward \
ALLOW_NON_EUCLIDEAN_WARD=1 \
bash run_hamming_pipeline.sh cluster B_cells__IR_vs_NR
```

Ward 结果必须标注为“reference-compatible sensitivity analysis”，不能作为默认主结果。

SciPy 官方同样明确说明 `ward` 只有在欧氏 pairwise metric 下才正确定义：

<https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html>

## 3. 默认输入

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/
├── qc_minmeth55_maxmethnone_maxsites10000000/
│   └── filtered_data_merged_30k/
│       ├── column_header.txt
│       └── chr*.npz
└── methdiff_30k/
    ├── metadata/cell_metadata.tsv
    └── results/response/
        ├── B_cells__IR_vs_NR_DMRs.bed
        ├── CD14_Monocytes__IR_vs_NR_DMRs.bed
        └── ...
```

## 4. 默认 DMR 特征集

按当前研究设计，主分析默认筛选：

```text
chr1–chr22
raw p < 0.01
不额外限制甲基化差异大小
n_sites >= 1
至少保留5个DMR
```

对应环境变量：

```text
P_COLUMN=raw
P_CUTOFF=0.01
ABS_DIFF=0
MIN_DMR_SITES=1
MIN_DMRS=5
CHROMOSOMES=autosomes
```

若要做 adjusted-p 严格敏感性分析，必须使用独立标签：

```bash
P_COLUMN=adjusted \
P_CUTOFF=0.05 \
ABS_DIFF=0.10 \
MIN_DMR_SITES=3 \
ANALYSIS_LABEL=adjustedp005_diff010 \
bash run_hamming_pipeline.sh run B_cells__IR_vs_NR
```

raw-p 与 adjusted-p 特征集不能混在同一目录。raw p 主分析的统计局限必须在结果解释中明确说明。

## 5. 细胞规模和缺失值

Hamming distance 矩阵的内存和计算量为 `O(n_cells²)`。不能直接对 51,024 个已注释细胞构建全距离矩阵。

默认对每个 cell type 最多选择 2,000 个细胞，并在 `IR/NR × sample` 十个 strata 间轮流抽样，避免大样本或单个 donor 主导：

```text
MAX_CELLS=2000
RANDOM_SEED=20260804
```

CpG和细胞覆盖默认门槛：

```text
MIN_SITE_CELLS=2    # 一个CpG至少在2个入选细胞中观测
MIN_CELL_SITES=5    # 一个细胞至少观测5个入选CpG
MIN_SHARED_SITES=1  # 任意保留细胞对至少共享1个CpG
```

若 overlap pruning 删除大量细胞，说明入选 DMR 提供的共同覆盖不足，不能通过把“无重叠”距离强行填成 0 或 1 来解决。应检查 DMR 特征数量、每细胞覆盖，或把区域级距离作为另一个明确标注的分析。

## 6. 运行逻辑

```text
[1/8] CHECK     上游DMR、metadata和filtered NPZ
[2/8] PREPARE   筛选raw p < 0.01 DMR，按response+sample平衡抽样
[3/8] EXTRACT   提取DMR内CpG的observed/methylated稀疏矩阵
[4/8] DISTANCE  共同观测CpG上的pairwise Hamming distance
[5/8] CLUSTER   average-linkage层次聚类、图和审计表
[6/8] MDS       预计算Hamming距离降至10维
[7/8] UMAP      MDS坐标上的10邻居图，min_dist=0.1
[8/8] LEIDEN    在同一邻接图上聚类，resolution=0.5
```

帮助和状态不会启动计算：

```bash
cd /share/home/rzli/METHSCAN/03_Hamming_distance

bash run_hamming_pipeline.sh --help
bash run_hamming_pipeline.sh status
bash run_hamming_pipeline.sh status B_cells__IR_vs_NR
```

## 7. B-cell pilot

上传脚本后先检查：

```bash
cd /share/home/rzli/METHSCAN/03_Hamming_distance

chmod 750 run_hamming_pipeline.sh hamming_scwgbs.py
bash -n run_hamming_pipeline.sh
python -m py_compile hamming_scwgbs.py
python - <<'PY'
import sklearn
import umap
import igraph
import leidenalg
print("sklearn", sklearn.__version__)
print("umap", umap.__version__)
print("igraph", igraph.__version__)
print("leidenalg: OK")
PY
bash run_hamming_pipeline.sh --help
```

提交完整 pilot：

```bash
mkdir -p scheduler_logs

dsub \
  -n hamming_30k_Bcell \
  -R "cpu=16;mem=65536MB" \
  --cwd /share/home/rzli/METHSCAN/03_Hamming_distance \
  -oo scheduler_logs/hamming_30k_Bcell.out \
  -eo scheduler_logs/hamming_30k_Bcell.err \
  bash run_hamming_pipeline.sh run B_cells__IR_vs_NR
```

监控：

```bash
tail -F \
  scheduler_logs/hamming_30k_Bcell.out \
  scheduler_logs/hamming_30k_Bcell.err
```

不要在登录节点直接运行 `prepare/extract/cluster/reduce/run`。`run` 会依次完成八步；已有 Hamming 和层次聚类结果时，也可单独提交 `reduce`。

若层次聚类已经完成，只提交 MDS–UMAP–Leiden：

```bash
dsub \
  -n hamming_30k_Bcell_reduce \
  -R "cpu=4;mem=16384MB" \
  --cwd /share/home/rzli/METHSCAN/03_Hamming_distance \
  -oo scheduler_logs/hamming_30k_Bcell_reduce.out \
  -eo scheduler_logs/hamming_30k_Bcell_reduce.err \
  bash run_hamming_pipeline.sh reduce B_cells__IR_vs_NR
```

## 8. 输出结构

默认分析标签由筛选参数自动生成：

```text
/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/
└── hamming_distance_30k/
    └── raw_p0p01_diff0_dmrsites1_autosomes_maxcells2000_sitecells2/
        └── B_cells__IR_vs_NR/
            ├── features/
            │   ├── selected_dmrs.bed
            │   ├── selected_dmrs.tsv
            │   ├── selected_cells.tsv
            │   ├── selected_cells.txt
            │   └── prepare_summary.tsv
            ├── cpg_calls/
            │   ├── observed_calls.npz
            │   ├── methylated_calls.npz
            │   ├── cells.tsv
            │   ├── sites.tsv.gz
            │   └── extraction_summary.tsv
            ├── clustering_average_shared1_cellsites5_k2/
            │   ├── hamming_distance.npy
            │   ├── shared_CpG_counts.npy
            │   ├── linkage.npy
            │   ├── cell_clusters.tsv
            │   ├── removed_cells.tsv
            │   ├── cluster_by_response.tsv
            │   ├── cluster_by_sample.tsv
            │   ├── clustering_summary.tsv
            │   ├── dendrogram.png
            │   ├── ordered_distance_heatmap.png
            │   └── mds10_umap_neighbors10_mindist0p10_leiden0p5_seed20260804/
            │       ├── mds_coordinates.npy
            │       ├── mds_coordinates.tsv
            │       ├── umap_coordinates.tsv
            │       ├── reduction_summary.tsv
            │       ├── leiden_by_response.tsv
            │       ├── leiden_by_sample.tsv
            │       ├── umap_by_response.png
            │       ├── umap_by_sample.png
            │       ├── umap_by_leiden.png
            │       └── umap_by_cluster.png
            └── markers/
```

## 9. 结果解释边界

这里使用的 DMR 是在同一批 IR/NR 细胞上发现的，再用这些 DMR 聚类相同细胞，因此属于监督式特征选择。若聚类分开 IR 和 NR，只能说明这些被预选的差异区域能够重现标签，不能作为独立、无监督验证。

MDS、UMAP 和 Leiden 都没有接收 response 标签；但输入距离来自按 IR/NR 筛出的 DMR，因此 `umap_by_response.png` 仍是监督式特征选择后的可视化，不能据此声称“无需去 batch”或完成独立验证。未观测 CpG 只在 Hamming 阶段通过 observed mask 处理，不会被 MDS 当作甲基化0。

`reduction_summary.tsv` 中的 `mds_normalized_stress` 衡量10维 MDS 对原始 Hamming 距离的保持程度，越小越好。UMAP 和 Leiden 使用同一张基于 MDS 坐标的 fuzzy neighbor graph；Leiden不是在二维 UMAP 坐标上重新计算距离。

至少同时报告：

- 聚类保留/删除的细胞数；
- 每对细胞共同 CpG 数的 min/median/max；
- hierarchical cluster 与 response 的 ARI/NMI；
- hierarchical cluster 与 sample 的 ARI；
- Leiden 与 response 的 ARI/NMI；
- Leiden 与 sample 的 ARI；
- MDS normalized stress；
- 各 donor 在 cluster 中的分布；
- leave-one-donor-out 或独立队列验证。

当 response ARI 高但 sample ARI 也高时，需要优先排查 donor/batch 驱动，而不能直接解释为疗效相关甲基化结构。
