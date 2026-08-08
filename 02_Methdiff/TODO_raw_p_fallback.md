# Meth diff 失败比较补跑 TODO

更新日期：2026-08-05

## 1. 目标

对 MethSCAn `diff` 因置换 DMR 数为 0 而在 FDR 计算阶段触发以下异常的比较进行定向补跑：

```text
ZeroDivisionError: division by zero
methscan/diff.py: calc_fdr(output_final[11] == "real")
```

补跑时保留真实分组 DMR 的 raw p，将无法估计的 adjusted p/FDR 明确记为 `NA`，不将其伪造为 0 或 1。

## 2. 当前状态

### `sample_celltype`

- [x] 总比较数：1,036
- [x] `eligible=yes`：774
- [x] 成功完成：723
- [x] 失败：51
- [x] 51/51 失败日志均包含 `ZeroDivisionError: division by zero`
- [x] `partial=0`，没有已写入但未验证的 BED
- [x] 调度作业 `163085` 已结束；因存在失败子比较，总作业状态为 `FAILED`

### `cross_response`

- [x] 总比较数：370
- [x] `eligible=yes`：303
- [x] 调度作业 `163087` 已提交
- [ ] 等待 `163087` 结束
- [ ] 记录 `complete`、`partial`、失败总数
- [ ] 统计其中有多少失败属于同一 `ZeroDivisionError`

## 3. 当前不做的操作

- [ ] **不**直接重提原始 `sample_celltype` 全量命令，否则51个比较会再次扫描后仍在 FDR 阶段失败。
- [ ] **不**删除723个已成功的 DMR BED 和 `.ok` 标记。
- [ ] **不**为了规避软件异常而降低 `MIN_CELLS=10`。
- [ ] **不**在 `163087` 运行期间修改用户 site-packages 中的 MethSCAn 源码。
- [ ] **不**把无法估计的 adjusted p 写成显著的 `0`。

## 4. 等 `cross_response` 结束后的审计

- [ ] 查看两个 mode 的状态：

```bash
cd /share/home/rzli/METHSCAN/02_Methdiff

bash run_methdiff_pipeline.sh status sample_celltype
bash run_methdiff_pipeline.sh status cross_response
```

- [ ] 统计 `sample_celltype` 失败比较：

```bash
grep '\[5/6 FAIL\] sample_celltype/' \
  scheduler_logs/methdiff_30k_sample_celltype_6x16.err |
sort -u |
wc -l
```

- [ ] 统计 `cross_response` 失败比较：

```bash
grep '\[5/6 FAIL\] cross_response/' \
  scheduler_logs/methdiff_30k_cross_response_6x16.err |
sort -u |
wc -l
```

- [ ] 确认两个 mode 的 FDR 除零日志数：

```bash
grep -l 'ZeroDivisionError: division by zero' \
  /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/logs/sample_celltype/*.log |
wc -l

grep -l 'ZeroDivisionError: division by zero' \
  /share/LCZX_Data/data/allcools/merged_10samples_upstream_v2/methdiff_30k/logs/cross_response/*.log |
wc -l
```

- [ ] 检查是否还存在其他错误类型，不将内存、输入文件或其他异常误归类为 raw-p fallback。

## 5. 开发 raw-p fallback

- [ ] 为 `run_methdiff_pipeline.sh` 增加独立 fallback 功能，只处理日志明确匹配下列条件的失败比较：

```text
calc_fdr(output_final[11] == "real")
ZeroDivisionError: division by zero
```

- [ ] fallback 不全局修改当前 Conda/用户 Python 环境中的 MethSCAn。
- [ ] fallback 使用与原运行完全相同的数据、group CSV、`MIN_CELLS=10`、窗口、阈值、随机种子和线程参数。
- [ ] 在置换 DMR 数为 0 时，保存真实分组 DMR 的 raw p，并将 adjusted p 写为 `NA`。
- [ ] 保持输出 BED 为12列，与已成功 MethSCAn DMR BED 的前11列语义一致。
- [ ] 修改 BED 校验：仅允许经 fallback 记录的结果在第12列使用 `NA`。
- [ ] 修改 summary：raw p 仍正常统计；adjusted-p 显著数对 `NA` 行不计数，并新增 `fdr_status`。
- [ ] fallback 完成标记至少记录：

```text
mode
comparison
group_A_label / group_B_label
group_A_n / group_B_n
group_file_sha256
DMR_rows
raw_p_fallback=yes
fdr_status=no_null_dmrs
adjusted_p=NA
```

- [ ] 在小型模拟数据上验证：常规路径不变，仅 `n_t_null=0` 时进入 fallback。

## 6. 生成精确补跑清单

- [ ] 从调度日志提取失败 comparison 名称。
- [ ] 与各 mode 的 `comparisons.tsv` 校验，确保每个名称唯一、`eligible=yes`且 group CSV 存在。
- [ ] 排除已经有有效 BED + `.ok` 标记的比较。
- [ ] 分别生成：

```text
sample_celltype_failed_no_null_dmrs.tsv
cross_response_failed_no_null_dmrs.tsv
```

- [ ] 补跑清单应包含 mode、comparison、group file、A/B label、A/B 细胞数和原始错误日志路径。

## 7. 定向补跑

- [ ] 先选1个 `sample_celltype` 失败比较运行 pilot。
- [ ] 验证 pilot 产生12列 BED，第11列 raw p 为数值，第12列为 `NA`。
- [ ] 验证 pilot 完成标记中含 `fdr_status=no_null_dmrs`。
- [ ] pilot 通过后，仅并行补跑失败清单中的比较。
- [ ] 先完成 `sample_celltype`，再补跑 `cross_response`。
- [ ] 保留补跑的独立 scheduler log，不覆盖原始失败日志。

> 补跑命令将在 fallback 脚本完成和 pilot 验证后填写，不提前给出未验证的命令。

## 8. 验收标准

- [ ] `sample_celltype` 的774个 eligible 比较全部有可审计结果：723个原始成功结果 + 51个 raw-p fallback 结果。
- [ ] `cross_response` 的303个 eligible 比较全部有可审计结果，或对其他错误类型有明确排除记录。
- [ ] 已有成功 BED 的 SHA-256 在补跑前后保持不变。
- [ ] 所有 fallback BED 的前11列通过格式和数值检查。
- [ ] raw p 汇总至少输出 `raw_p_lt_0.05` 和 `raw_p_lt_0.01`。
- [ ] adjusted-p 统计明确区分 `estimated` 与 `not_estimable_no_null_dmrs`。
- [ ] 最终 summary 不把 `NA` adjusted p 当成0、1或显著结果。
- [ ] README 记录错误原因、fallback 方法、结果限制和补跑作业资源。

## 9. 统计解释边界

- raw p 可用于探索性排序和后续外部统一多重检验，但不等于已经完成 FDR 控制。
- `adjusted_p=NA` 表示 MethSCAn 的置换空分布中没有 DMR，因而当前方法无法估计 FDR；它不表示 adjusted p 为0。
- 对 raw p 再做 BH/FDR 时，必须事先确定合理的检验家族（单个 comparison、单个 mode 或其他预定义范围），不在结果出来后为了获得更多显著结果任意切换。
