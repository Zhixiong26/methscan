# clean cell-type DMR q<0.05 去 sample/response component 分析流程与脚本整理

## 0. 分析目的

老板说的步骤是：

```text
Meth_diff 先找细胞类型间的 DMRs
再找到样本组件的 DMRs
用细胞类型间的 DMRs 减去样本组件的 DMRs
得到去过 batch / sample component 的 DMRs
```

你已经完成的对应分析是：

```text
cell-type pairwise DMRs q<0.05
-
same-cell-type IR vs NR DMRs q<0.05
=
clean cell-type DMRs q<0.05
```

这里的“去 batch”更准确地说是：

```text
在 DMR feature 层面去除 sample / response component
```

不是 Harmony / MethylVI 那种 embedding 或 matrix 层面的 batch correction。

---

# 1. 输入数据路径

## 1.1 细胞类型间 DMRs 输入

路径：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/1_all_cells_cell_type_pairwise
```

示例文件：

```text
B_cells_vs_CD4_T_cells_DMRs.bed
B_cells_vs_CD8_T_cells_DMRs.bed
CD4_T_cells_vs_NK_cells_DMRs.bed
Monocytes_CD14_vs_NK_cells_DMRs.bed
...
```

这些文件代表：

```text
所有细胞中，不同 cell type 之间的 pairwise DMRs
```

---

## 1.2 sample / response component DMRs 输入

路径：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR
```

示例文件：

```text
B_cells_IR_vs_NR_DMRs.bed
CD4_T_cells_IR_vs_NR_DMRs.bed
CD8_T_cells_IR_vs_NR_DMRs.bed
Monocytes_CD14_IR_vs_NR_DMRs.bed
Monocytes_CD16_IR_vs_NR_DMRs.bed
NK_cells_IR_vs_NR_DMRs.bed
pDCs_IR_vs_NR_DMRs.bed
Plasma_cells_IR_vs_NR_DMRs.bed
```

这些文件代表：

```text
同一 cell type 内 IR vs NR 的 DMRs
```

在老板这一步中，它们被当作 sample / response component DMRs 来减掉。

---

## 1.3 DMR 文件关键列

Meth_diff 输出 BED 为 12 列。已经确认：

```text
第 10 列：direction / hypomethylated group
第 11 列：p value
第 12 列：q value / adjusted p value
```

q 值筛选使用：

```bash
$12 < 0.05
```

确认命令：

```bash
awk 'BEGIN{OFS="\t"} {print "NF="NF, $1,$2,$3,"direction="$10,"p="$11,"q="$12; exit}' \
/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR/CD4_T_cells_IR_vs_NR_DMRs.bed
```

之前输出示例：

```text
NF=12 chr9 113058368 113065368 direction=group_B p=0.0 q=0.0
```

---

## 1.4 All VMR matrix region 输入

All VMR matrix region BED：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.bed
```

All VMR matrix region map：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.map.tsv
```

原始 All VMR matrix：

```bash
/share/LCZX_Data/data/All/VMR_matrix
```

其中关键矩阵：

```bash
/share/LCZX_Data/data/All/VMR_matrix/mean_shrunken_residuals.csv.gz
```

---

# 2. Step 1：q<0.05 筛选、merge、subtract、映射回 All VMR

## 2.1 脚本路径

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/1_merge_subtract_DMR_200k_q005.sh
```

## 2.2 脚本内容

```bash
#!/bin/bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

BASE_DIR="/share/home/rzli/METHSCAN/Meth_diff"
CELLTYPE_DMR_DIR="${BASE_DIR}/DMR_results_200k/1_all_cells_cell_type_pairwise"
SAMPLE_DMR_DIR="${BASE_DIR}/DMR_results_200k/3_same_cell_type_IR_vs_NR"
OUT_DIR="${BASE_DIR}/DMR_clean_200k_q005"
MAP_DIR="${OUT_DIR}/matrix_mapping"

ALL_VMR_BED="/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.bed"
TOP5="/share/home/rzli/METHSCAN/TopVMR_individual_analysis/results/top_vmr_from_all_scan/VMRs_top5pct.bed"
TOP2="/share/home/rzli/METHSCAN/TopVMR_individual_analysis/results/top_vmr_from_all_scan/VMRs_top2pct.bed"

mkdir -p "${OUT_DIR}" "${MAP_DIR}"

# 1. cell-type pairwise DMRs: q<0.05 后 merge
cat "${CELLTYPE_DMR_DIR}"/*.bed \
  | awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0) < 0.05 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${OUT_DIR}/cell_type_pairwise_DMRs_q005_merged.bed"

# 2. sample / response component DMRs: same-cell-type IR vs NR q<0.05 后 merge
cat "${SAMPLE_DMR_DIR}"/*.bed \
  | awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0) < 0.05 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q005_merged.bed"

# 3. cell-type DMRs - sample component DMRs
bedtools intersect -v \
  -a "${OUT_DIR}/cell_type_pairwise_DMRs_q005_merged.bed" \
  -b "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q005_merged.bed" \
  > "${OUT_DIR}/cell_type_DMRs_without_sample_component_q005.bed"

# 4. 映射回 All VMR matrix regions
bedtools intersect -u \
  -a "${ALL_VMR_BED}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q005.bed" \
  > "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q005.bed"

cut -f4 "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q005.bed" \
  > "${MAP_DIR}/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt"

# 5. 同时检查 top5 / top2 VMR overlap
bedtools intersect -u \
  -a "${TOP5}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q005.bed" \
  > "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q005.bed"

bedtools intersect -u \
  -a "${TOP2}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q005.bed" \
  > "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q005.bed"

awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' \
  "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q005.bed" \
  > "${MAP_DIR}/top5_clean_cell_type_DMR_q005_regions.txt"

awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' \
  "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q005.bed" \
  > "${MAP_DIR}/top2_clean_cell_type_DMR_q005_regions.txt"

# 6. summary
{
  echo "cell-type DMRs q<0.05 merged:"
  wc -l "${OUT_DIR}/cell_type_pairwise_DMRs_q005_merged.bed"
  echo "sample-component DMRs q<0.05 merged:"
  wc -l "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q005_merged.bed"
  echo "clean cell-type DMRs q<0.05:"
  wc -l "${OUT_DIR}/cell_type_DMRs_without_sample_component_q005.bed"
  echo "All VMR overlap clean q<0.05 DMR:"
  wc -l "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q005.bed"
  echo "clean q<0.05 region names:"
  wc -l "${MAP_DIR}/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt"
  echo "top5 overlap clean q<0.05 DMR:"
  wc -l "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q005.bed"
  echo "top2 overlap clean q<0.05 DMR:"
  wc -l "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q005.bed"
}
```

## 2.3 提交命令

```bash
dsub -n 1_merge_subtract_DMR_200k_q005 \
  -R "cpu=64;mem=180G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/merge_subtract_DMR_200k_q005.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/merge_subtract_DMR_200k_q005.%J.err \
  bash /share/home/rzli/METHSCAN/Meth_diff/scripts/merge_subtract_DMR_200k_q005.sh
```

JobID：

```text
162594
```

## 2.4 日志路径

```bash
/share/home/rzli/METHSCAN/Meth_diff/logs/merge_subtract_DMR_200k_q005.162594.out
/share/home/rzli/METHSCAN/Meth_diff/logs/merge_subtract_DMR_200k_q005.162594.err
```

## 2.5 输出结果路径

主结果目录：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005
```

细胞类型间 q<0.05 merged DMR：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/cell_type_pairwise_DMRs_q005_merged.bed
```

sample / response component q<0.05 merged DMR：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/sample_component_IR_vs_NR_DMRs_q005_merged.bed
```

clean cell-type DMR：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/cell_type_DMRs_without_sample_component_q005.bed
```

映射回 All VMR 的 region list：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/matrix_mapping/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt
```

## 2.6 数量结果

```text
cell-type DMRs q<0.05 merged:          149867
sample-component DMRs q<0.05 merged:     4353
clean cell-type DMRs q<0.05:           147974
All VMR overlap clean q<0.05 DMR:       29863
top5 overlap clean q<0.05 DMR:           1254
top2 overlap clean q<0.05 DMR:            543
```

---

# 3. Step 2：根据 clean DMR region list 抽取 All VMR matrix

## 3.1 Python 抽矩阵脚本路径

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/2_subset_All_VMR_matrix_by_region_list.py
```

## 3.2 Python 脚本内容

```python
#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

MATRIX_FILES = [
    "methylation_fractions.csv.gz",
    "methylated_sites.csv.gz",
    "total_sites.csv.gz",
    "mean_shrunken_residuals.csv.gz",
]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--region-list", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    region_list = Path(args.region_list)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(region_list, header=None)[0].astype(str).tolist()
    region_set = set(regions)
    print(f"requested regions: {len(regions)}")

    for fname in MATRIX_FILES:
        infile = input_dir / fname
        outfile = output_dir / fname

        print(f"Reading: {infile}")
        df = pd.read_csv(infile, index_col=0)

        keep = [c for c in df.columns if c in region_set]
        missing = len(region_set) - len(keep)

        print(f"{fname} kept regions: {len(keep)} missing: {missing}")

        df_sub = df.loc[:, keep]
        df_sub.to_csv(outfile, compression="gzip")
        print(f"Written: {outfile}")

if __name__ == "__main__":
    main()
```

> 注：这是该脚本的核心逻辑版本。实际脚本路径如上，若需要以服务器真实内容为准，可直接 `cat` 该文件确认。

---

## 3.3 wrapper 脚本路径

```bash
/share/home/rzli/METHSCAN/Meth_diff/scripts/3_run_subset_all_clean_DMR_matrix_200k_q005.sh
```

## 3.4 wrapper 脚本内容

```bash
#!/bin/bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

export PYTHONNOUSERSITE=1
unset PYTHONPATH
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

python /share/home/rzli/METHSCAN/Meth_diff/scripts/3_subset_All_VMR_matrix_by_region_list.py \
  --input-dir "/share/LCZX_Data/data/All/VMR_matrix" \
  --region-list "/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/matrix_mapping/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt" \
  --output-dir "/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/VMR_matrix_all_clean_cell_type_DMR_q005"
```

## 3.5 提交命令

```bash
dsub -n 3_subset_all_clean_DMR_matrix_200k_q005 \
  -R "cpu=64;mem=120G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/3_subset_all_clean_DMR_matrix_200k_q005.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/3_subset_all_clean_DMR_matrix_200k_q005.%J.err \
  bash /share/home/rzli/METHSCAN/Meth_diff/scripts/3_run_subset_all_clean_DMR_matrix_200k_q005.sh
```

JobID：

```text
162595
```

## 3.6 输出结果路径

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/VMR_matrix_all_clean_cell_type_DMR_q005
```

包含：

```text
methylation_fractions.csv.gz
methylated_sites.csv.gz
total_sites.csv.gz
mean_shrunken_residuals.csv.gz
```

关键矩阵：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/VMR_matrix_all_clean_cell_type_DMR_q005/mean_shrunken_residuals.csv.gz
```

## 3.7 抽取结果

```text
requested regions: 29863
methylation_fractions.csv.gz kept regions: 29863 missing: 0
methylated_sites.csv.gz kept regions: 29863 missing: 0
total_sites.csv.gz kept regions: 29863 missing: 0
mean_shrunken_residuals.csv.gz kept regions: 29863 missing: 0
```

---

# 4. Step 3：用 clean DMR matrix 重新做 annotation / PCA / UMAP

## 4.1 clean annotation R 脚本路径

```bash
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/Downstream_annotation_clean_DMR_by_threshold.R
```

这个脚本由原始 annotation 脚本复制修改而来：

原脚本：

```bash
/share/home/rzli/METHSCAN/Annotation/common/Downstream_annotation_by_threshold.R
```

clean 版主要修改为支持：

```text
METHSCAN_MATRIX_DIR
METHSCAN_OUTPUT_DIR
```

以避免覆盖原始 MethSCAn 结果。

---

## 4.2 clean annotation wrapper：resolution=0.0001

脚本路径：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005.sh
```

脚本核心内容：

```bash
#!/bin/bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

export PYTHONNOUSERSITE=1
unset PYTHONPATH
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

SCRIPT="/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/Downstream_annotation_clean_DMR_by_threshold.R"
BASE="/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005"

export METHSCAN_MATRIX_DIR="${BASE}/VMR_matrix_all_clean_cell_type_DMR_q005"
export METHSCAN_OUTPUT_DIR="/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005/result"

export METHSCAN_THRESHOLD="200k"
export METHSCAN_N_PCS="20"
export METHSCAN_UMAP_N_NEIGHBORS="50"
export METHSCAN_UMAP_MIN_DIST="0.3"
export METHSCAN_LEIDEN_RESOLUTION="0.0001"

Rscript "${SCRIPT}"
```

提交命令：

```bash
dsub -n clean_all_DMR_annotation_200k_q005 \
  -R "cpu=16;mem=96G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_annotation_200k_q005.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_annotation_200k_q005.%J.err \
  bash /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005.sh
```

JobID：

```text
162596
```

输出目录：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005/result
```

结果说明：

```text
resolution=0.0001 时 Leiden 只有 1 个 cluster。
```

---

## 4.3 clean annotation wrapper：resolution=0.001

脚本路径：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0001.sh
```

生成命令：

```bash
cp /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005.sh \
   /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0001.sh

sed -i 's#clean_all_DMR_200k_q005/result#clean_all_DMR_200k_q005_res0001/result#g' \
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0001.sh

sed -i 's/METHSCAN_LEIDEN_RESOLUTION="0.0001"/METHSCAN_LEIDEN_RESOLUTION="0.001"/g' \
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0001.sh
```

提交命令：

```bash
dsub -n clean_all_DMR_q005_res0001 \
  -R "cpu=16;mem=96G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_q005_res0001.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_q005_res0001.%J.err \
  bash /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0001.sh
```

JobID：

```text
162597
```

输出目录：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005_res0001/result
```

结果说明：

```text
resolution=0.001 时 Leiden cluster 过多，约 100+ clusters。
```

---

## 4.4 clean annotation wrapper：resolution=0.0005

脚本路径：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0005.sh
```

生成命令：

```bash
cp /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005.sh \
   /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0005.sh

sed -i 's#clean_all_DMR_200k_q005/result#clean_all_DMR_200k_q005_res0005/result#g' \
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0005.sh

sed -i 's/METHSCAN_LEIDEN_RESOLUTION="0.0001"/METHSCAN_LEIDEN_RESOLUTION="0.0005"/g' \
/share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0005.sh
```

提交命令：

```bash
dsub -n clean_all_DMR_q005_res0005 \
  -R "cpu=64;mem=128G" \
  -o /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_q005_res0005.%J.out \
  -e /share/home/rzli/METHSCAN/Meth_diff/logs/clean_all_DMR_q005_res0005.%J.err \
  bash /share/home/rzli/METHSCAN/Annotation/clean_DMR_200k/run_all_clean_DMR_annotation_200k_q005_res0005.sh
```

JobID：

```text
162598
```

输出目录：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005_res0005/result
```

结果说明：

```text
0.0005 和 0.001 的 UMAP/PCA 形态差不多，主要区别只在 Leiden cluster 数量。
Leiden 不建议作为主结果。
```

---

# 5. Step 4：检查 q<0.05 clean DMR annotation 结果

## 5.1 主要输出图目录

resolution=0.0001：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005/result/plots
```

resolution=0.001：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005_res0001/result/plots
```

resolution=0.0005：

```bash
/share/home/rzli/METHSCAN/Annotation/clean_all_DMR_200k_q005_res0005/result/plots
```

## 5.2 主要图类型

```text
ALL_PCA_200k.png
ALL_umap_plot_by_cell_type_200k.jpeg
ALL_umap_plot_by_response_200k.jpeg
ALL_umap_plot_by_sample_200k.jpeg
ALL_umap_plot_by_leiden_200k.jpeg
ALL_umap_plot_response_by_cell_type_200k.png
```

## 5.3 结果解释

当前 clean DMR feature set 的图显示：

```text
1. UMAP by cell_type：大类 cell type 有一定结构。
2. UMAP by response：IR/NR 整体大体混合，不是主要分离轴。
3. UMAP by sample：仍存在局部 sample 富集，说明 DMR subtract 降低 sample/response component，但不是完整 matrix-level batch correction。
4. Leiden：resolution 很敏感，不适合作为主结果。
```

因此建议主图使用：

```text
UMAP by cell_type
UMAP by response
UMAP by sample
response within each cell type
```

不建议把 Leiden cluster 作为核心结果。

---

# 6. 最终结果总结

## 6.1 老板要求步骤是否完成

已完成。

对应关系：

```text
老板说：细胞类型间 DMRs
对应：/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/1_all_cells_cell_type_pairwise

老板说：样本组件 DMRs
对应：/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR

老板说：细胞类型间 DMRs - 样本组件 DMRs
对应：/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/cell_type_DMRs_without_sample_component_q005.bed
```

## 6.2 最重要结果文件

clean cell-type DMR：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/cell_type_DMRs_without_sample_component_q005.bed
```

映射回 All VMR 的 clean feature list：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/matrix_mapping/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt
```

抽取后的 clean DMR matrix：

```bash
/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005/VMR_matrix_all_clean_cell_type_DMR_q005/mean_shrunken_residuals.csv.gz
```

## 6.3 关键数量

```text
cell-type DMRs q<0.05 merged:          149867
sample-component DMRs q<0.05 merged:     4353
clean cell-type DMRs q<0.05:           147974
mapped All VMR features:                29863
top5 overlap clean DMR:                  1254
top2 overlap clean DMR:                   543
```

---

# 7. 可以回复老板的文字

```text
这一步已经完成。我先用 Meth_diff 的 cell-type pairwise DMRs 作为细胞类型间 DMRs，并按 q value < 0.05 筛选；然后用 same-cell-type IR vs NR DMRs 作为样本/response component DMRs，同样按 q value < 0.05 筛选；最后用 bedtools intersect -v 做 cell-type DMRs 减去 sample-component DMRs，得到 clean cell-type DMRs。

结果是：
cell-type DMRs q<0.05 merged: 149867
sample-component DMRs q<0.05 merged: 4353
subtract 后 clean cell-type DMRs q<0.05: 147974
这些 clean DMRs 映射回 All VMR matrix 后对应 29863 个 VMR features。
```

---

# 8. 注意事项

1. 这里的“去 batch”是 DMR feature 层面去除 sample / response component，不是 MethylVI/Harmony 这种 latent embedding 层面的 batch correction。

2. 这个 clean cell-type DMR 结果主要用于验证 cell-type structure 和构建 clean feature set。

3. 后续 disease promoter DMR / RNA expression integration 是另一条分析线，输入是 same-cell-type IR vs NR DMRs，不能再把这部分减掉，否则会把疾病差异本身删掉。

