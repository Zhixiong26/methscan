# MethScan ALL 200K 下游注释流程

## 1. 用途

本目录保存 2026-07-16 版本的 MethScan ALL 200K 下游分析脚本。流程将 IR01–IR05、NR01–NR05 共 10 个样本的合并甲基化矩阵统一执行：

1. 与新版 Scanpy 注释匹配，并在 PCA 前排除未匹配细胞；
2. 缺失值迭代填补和 PCA；
3. UMAP 与 Leiden 聚类；
4. 合并样本、IR/NR 分组和细胞类型注释；
5. 导出坐标、统计表、图片及各细胞类型的 IR/NR cell-group 文件。

这是 ALL 合并分析，不是 10 个样本分别运行的单样本流程。

## 2. 脚本关系

~~~text
01_run_All_200k.sh
        ↓
02_All_200k_analysis.R 200k
~~~

- 01_run_All_200k.sh：激活 scDNAm 环境、设置线程和分析参数、调用 R 脚本、写入日志。
- 02_All_200k_analysis.R：读取和过滤数据，运行 PCA、UMAP、Leiden，合并注释并保存结果。

服务器运行目录：

~~~text
/share/home/rzli/METHSCAN/Annotation/20260716
~~~

本地归档目录：

~~~text
/Users/luozhixiong/Documents/PHD/脚本/Methscan/20260716/scripts
~~~

服务器最终提交版本使用 32 线程。如果本地脚本仍显示 1 线程，说明本地副本早于服务器最终版本，应重新下载后归档。

## 3. 输入文件

### ALL 200K 甲基化矩阵

~~~text
/share/LCZX_Data/data/All/VMR_matrix/mean_shrunken_residuals.csv.gz
~~~

说明：

- VMR_matrix 是早期 20W/200K 配置对应的目录；
- 最终矩阵约有 83,245 个 VMR 特征；
- 矩阵细胞 ID 格式为 IR01__barcode。

### 新版 Scanpy 注释

~~~text
/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv
~~~

脚本使用字段：

~~~text
cell_id
sample
group
cell_type_integrated
~~~

Scanpy ID 格式为 IR01_barcode。脚本将其转换为 IR01__barcode 后与甲基化矩阵匹配。

提交前统计：

~~~text
甲基化矩阵细胞：16,241
匹配 Scanpy 注释：15,500
未匹配：741
总体匹配率：95.44%
~~~

741 个未匹配细胞在 PCA 前排除，后续分析使用约 15,500 个细胞。

## 4. 纳入排除口径

当前输入为 02_cell_annotation_all_cells.csv，因此：

- 排除没有新版 Scanpy 注释的细胞；
- 保留 Platelet_erythroid_contamination；
- 污染类型在细胞类型图中显示为灰色；
- 污染细胞仍参与 PCA、UMAP、Leiden 和 cell-group 文件生成。

如果要做严格 clean 分析，应改用：

~~~text
/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_clean_cells.csv
~~~

切换口径后应使用新的结果目录，避免与 all-cells 结果混合。

## 5. 参数与资源

分析参数：

~~~text
阈值：200k
PCA：20 PCs
UMAP n_neighbors：30
UMAP min_dist：0.05
UMAP seed：2
Leiden resolution：0.001
UMAP SGD：1线程（保证复现性）
~~~

2026-07-16 正式提交配置：

~~~text
OpenBLAS/OpenMP/MKL：32线程
UMAP邻居搜索：32线程
UMAP SGD：1线程
dsub：32 CPU，180G内存
~~~

Conda 环境：

~~~text
/share/home/rzli/miniconda3
环境名：scDNAm
~~~

主要 R 包：dplyr、tibble、ggplot2、irlba、uwot、igraph、data.table。

## 6. 提交前检查

进入目录：

~~~bash
cd /share/home/rzli/METHSCAN/Annotation/20260716
~~~

检查语法：

~~~bash
bash -n 01_run_All_200k.sh
~~~

~~~bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript -e 'invisible(parse(file="02_All_200k_analysis.R")); cat("Syntax OK\n")'
~~~

检查线程：

~~~bash
grep -nE 'THREADS|n_threads|n_sgd_threads' \
  01_run_All_200k.sh 02_All_200k_analysis.R
~~~

检查输入：

~~~bash
ls -lh /share/LCZX_Data/data/All/VMR_matrix/mean_shrunken_residuals.csv.gz
ls -lh /share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv
~~~

检查包接口：

~~~bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript -e 'library(uwot); library(igraph); cat("UMAP args:",all(c("seed","n_threads","n_sgd_threads","ret_nn") %in% names(formals(uwot::umap))),"\n"); cat("Leiden arg:","resolution_parameter" %in% names(formals(igraph::cluster_leiden)),"\n")'
~~~

UMAP args 和 Leiden arg 均应为 TRUE。

检查是否已有输出：

~~~bash
find result logs -maxdepth 3 -type f -ls
~~~

同名日志和结果在重跑时会被覆盖。需要保留旧结果时，应先归档或新建日期目录。

## 7. 使用 dsub 提交

先登录：

~~~bash
dlogin
~~~

提交：

~~~bash
dsub \
  -n ALL_200k_20260716 \
  -R "cpu=32;mem=180G" \
  --cwd /share/home/rzli/METHSCAN/Annotation/20260716 \
  -oo logs/dsub_ALL_200k_20260716.%J.out \
  -eo logs/dsub_ALL_200k_20260716.%J.err \
  bash 01_run_All_200k.sh
~~~

不要同时使用 nohup 或重复提交同一任务，否则多个进程会同时写入同一结果目录。

首次提交记录：

~~~text
日期：2026-07-16
Job ID：162776
任务名：ALL_200k_20260716
资源：32 CPU，180G内存
~~~

## 8. 监控

本服务器未安装 `bjobs`，不要用它查询 dsub 任务。任务是否已经开始，可通过核心日志判断；若需要查找站点提供的状态命令，可执行：

~~~bash
dsub -h
compgen -c | grep -E '^d(sub|job|jobs|stat|queue)' | sort -u
~~~

核心日志出现输入矩阵、匹配细胞数等信息，说明任务已进入实际运行阶段。最终以核心日志中的 `Done: 200k`、dsub 错误日志和结果文件完整性判断是否成功。

核心 R 日志：

~~~bash
tail -f /share/home/rzli/METHSCAN/Annotation/20260716/logs/ALL_200k_20260716.log
~~~

dsub 日志：

~~~bash
tail -f /share/home/rzli/METHSCAN/Annotation/20260716/logs/dsub_ALL_200k_20260716.162776.out
tail -f /share/home/rzli/METHSCAN/Annotation/20260716/logs/dsub_ALL_200k_20260716.162776.err
~~~

任务排队时日志可能不存在或为空。正常核心日志开头应包含：

~~~text
Threshold: 200k
Input matrix: /share/LCZX_Data/data/All/VMR_matrix/mean_shrunken_residuals.csv.gz
Output dir: /share/home/rzli/METHSCAN/Annotation/20260716/result
Scanpy-matched cells: 15500/16241
Loaded methylation matrix: 15500 cells x 83245 VMRs
~~~

## 9. 输出

输出目录：

~~~text
/share/home/rzli/METHSCAN/Annotation/20260716/result
~~~

主要文件：

~~~text
ALL_PCA_200k.RData
ALL_PCA_coordinates_200k.csv
ALL_UMAP_coordinates_200k.csv
ALL_annotation_200k.csv
ALL_cell_count_by_response_cell_type_200k.csv
ALL_cell_count_by_sample_cell_type_200k.csv
~~~

图片位于 result/plots/：

~~~text
ALL_PCA_200k.png
ALL_umap_plot_by_cell_type_200k.png
ALL_umap_plot_by_response_200k.png
ALL_umap_plot_response_by_cell_type_200k.png
ALL_umap_plot_by_sample_200k.png
ALL_umap_plot_by_leiden_200k.png
~~~

各细胞类型 IR/NR 分组文件位于：

~~~text
result/cell_groups_IR_vs_NR_by_cell_type/
~~~

## 10. 完成后核验

查看核心日志尾部：

~~~bash
tail -n 50 logs/ALL_200k_20260716.log
~~~

成功时应出现：

~~~text
Done: 200k
~~~

列出结果：

~~~bash
find result -maxdepth 3 -type f -ls
~~~

检查注释结果：

~~~bash
wc -l result/ALL_annotation_200k.csv
~~~

若 15,500 个匹配细胞全部进入结果，预期约为：

~~~text
15,501 行 = 15,500 个细胞 + 1 行表头
~~~

检查缺失注释：

~~~bash
grep 'Cells without cell type annotation' logs/ALL_200k_20260716.log
~~~

预期为 0。

## 11. 常见问题

### OpenBLAS 无法创建线程

若出现：

~~~text
OpenBLAS blas_thread_init: pthread_create failed
~~~

确认 dsub 申请 CPU 数与 OPENBLAS_NUM_THREADS、OMP_NUM_THREADS、MKL_NUM_THREADS 和 n_threads 一致。必要时降到 16 或 8 后重新提交。

### 内存不足

过滤后矩阵约为 15,500 × 83,245，单份 double 矩阵约 10.3 GB。读取、中心化、NA 位置矩阵、PCA 重建和临时对象会显著增加峰值内存。本次申请 180G。

如果任务被杀死，应检查 dsub 错误日志和调度详情，并优先优化缺失值填补过程，而不是继续增加线程。

### 细胞 ID 无法匹配

矩阵使用 IR01__barcode，Scanpy 使用 IR01_barcode。脚本通过 sample 和 substring() 构造双下划线 ID。若注释文件格式改变，需要同步修改转换逻辑。

### DONE 但结果不存在

检查：

~~~bash
cat logs/dsub_ALL_200k_20260716.162776.err
tail -n 100 logs/ALL_200k_20260716.log
~~~

启动脚本使用 set -euo pipefail。R 非零退出时不会打印完成信息。

## 12. 归档建议

每次正式分析建议使用独立日期目录，并保存：

- 实际运行的 Shell 和 R 脚本；
- Scanpy 注释文件路径和版本；
- dsub Job ID、CPU与内存配置；
- Conda 环境及主要 R 包版本；
- 全部日志；
- 参数与纳入排除标准；
- 结果文件校验信息。

不要修改历史日期目录后直接覆盖旧结果。若注释文件、参数或分析口径变化，应新建日期目录重新运行。
