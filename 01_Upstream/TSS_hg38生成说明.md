# 人类 hg38 TSS 参考文件生成说明

本文记录 MethSCAn 上游流程所用人类 TSS 参考文件的来源、坐标转换、验证方法和校验和，确保后续分析能够复现。

## 1. 最终文件

服务器路径：

```text
/share/LCZX_Data/ref/human_hg38_TSS.bed
```

文件格式：

```text
chrom  TSS  TSS  gene  .  strand
```

示例：

```text
chr1    11874    11874    DDX11L1      .    +
chr1    17436    17436    MIR6859-1    .    -
```

MethSCAn 通过 `--strand-column 6` 将第 6 列解释为链方向：

```bash
methscan profile \
  --strand-column 6 \
  /share/LCZX_Data/ref/human_hg38_TSS.bed \
  compact_data \
  TSS_profile.csv
```

## 2. 数据来源

基因注释来自 UCSC Genome Browser 的 GRCh38/hg38 `refFlat` 表：

- 数据目录：<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/>
- 注释表：<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz>
- 表结构：<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.sql>
- MethSCAn 教程：<https://anders-biostat.github.io/MethSCAn/tutorial.html>

下载日期：2026-08-03。

下载命令：

```bash
cd /share/LCZX_Data/ref

wget -c \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz

wget -c \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.sql
```

## 3. refFlat 字段

`refFlat.txt.gz` 没有表头。生成 TSS 时使用以下字段：

| 列 | 字段 | 含义 |
|---:|---|---|
| 1 | `geneName` | 基因名称 |
| 2 | `transcriptName` | 转录本名称 |
| 3 | `chrom` | 染色体 |
| 4 | `strand` | 链方向，`+` 或 `-` |
| 5 | `txStart` | 转录本起点，UCSC 0-based |
| 6 | `txEnd` | 转录本终点，UCSC 右端不包含 |

## 4. TSS 坐标转换

UCSC 使用 0-based、半开区间的 transcript 坐标。为匹配 MethSCAn 教程提供的点坐标格式，转换规则为：

```text
正链：TSS = txStart + 1
负链：TSS = txEnd
```

不能对正负链都使用 `txStart`。对于负链基因，`txStart` 是基因组坐标较小的一端，而真正的转录起始位点位于 `txEnd` 一侧。

例如：

```text
MIR6859-1  chr1  -  txStart=17368  txEnd=17436
```

正确 TSS 为 `17436`，不是 `17368`。

## 5. 生成命令

以下命令：

- 仅保留主染色体 `chr1–chr22`、`chrX`、`chrY` 和 `chrM`；
- 正确处理正负链；
- 输出 MethSCAn 所需六列；
- 删除同一基因、同一 TSS、同一链方向的完全重复记录；
- 保留不同基因共享同一 TSS 坐标的记录。

```bash
cd /share/LCZX_Data/ref

zcat refFlat.txt.gz |
awk 'BEGIN {FS=OFS="\t"}
$3 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ {
    tss = ($4 == "+") ? $5 + 1 : $6
    print $3, tss, tss, $1, ".", $4
}' |
sort -k1,1V -k2,2n -k4,4 -k6,6 -u \
> human_hg38_TSS.bed
```

## 6. 格式验证

先检查下载文件是否完整：

```bash
gzip -t refFlat.txt.gz
```

检查文件大小、前几行和记录数：

```bash
ls -lh refFlat.txt.gz refFlat.sql human_hg38_TSS.bed
head -n 10 human_hg38_TSS.bed
wc -l human_hg38_TSS.bed
```

检查六列、坐标和链方向：

```bash
awk '
NF != 6 ||
$2 !~ /^[0-9]+$/ ||
$3 !~ /^[0-9]+$/ ||
$2 != $3 ||
$2 < 1 ||
($6 != "+" && $6 != "-") {
    print "bad line:", NR, $0
}' human_hg38_TSS.bed | head
```

正常情况下该命令没有输出。

统计各染色体记录数：

```bash
cut -f1 human_hg38_TSS.bed |
sort -V |
uniq -c
```

## 7. 本次验证结果

2026-08-03 实际生成结果：

- `human_hg38_TSS.bed`：1.5 MB；
- 总记录数：42,024；
- 六列格式检查：通过；
- 坐标检查：通过；
- 链方向检查：通过；
- 染色体范围：`chr1–chr22`、`chrX`、`chrY` 和 `chrM`；
- compact 数据使用 `chr*.npz`，与 TSS 文件命名一致。

各染色体记录数：

```text
chr1   4296    chr2   2736    chr3   2325    chr4   1642
chr5   1894    chr6   2136    chr7   1995    chr8   1531
chr9   1663    chr10  1701    chr11  2374    chr12  1995
chr13   889    chr14  1333    chr15  1517    chr16  1720
chr17  2244    chr18   634    chr19  2521    chr20  1081
chr21   685    chr22   938    chrX   1898    chrY    275
chrM      1
```

compact 中还存在 `KI*`、`GL*` 等替代序列或未定位 scaffold。TSS 文件只覆盖主染色体，因此这些 NPZ 不参与 TSS profile；这不会影响主染色体的 TSS 质量控制。

## 8. 文件校验和

```text
cabd65c85ce0db017d771744e9db6ea80b3ef741594763730d6917756914d631  human_hg38_TSS.bed
8d6bb2d024c1c0c466ed5f2cb05eeb1c9d1ccba03d8c525ffac5a0fab790da1e  refFlat.txt.gz
35ab64f2b0ad0b40ec249b8e4143030c08812ab6672c46c885b1b9558650131c  refFlat.sql
```

复核命令：

```bash
cd /share/LCZX_Data/ref

sha256sum \
  human_hg38_TSS.bed \
  refFlat.txt.gz \
  refFlat.sql
```

## 9. 与标准 BED/bedtools 的区别

本项目的 `human_hg38_TSS.bed` 用于匹配 MethSCAn 教程格式，采用：

```text
start = end = 1-based TSS
```

它不是标准的 0-based、半开区间、1 bp BED。若文件用于 `bedtools intersect`、`bedtools window` 等操作，应另行生成标准 BED：

```text
正链：[txStart, txStart + 1)
负链：[txEnd - 1, txEnd)
```

不要将 MethSCAn 点坐标文件与标准 bedtools 文件混用。

## 10. Profile 来源追踪

当前上游脚本会计算 `human_hg38_TSS.bed` 的 SHA-256，并将 profile 来源写入每个样本的：

```text
TSS_profile_common.meta.tsv
```

示例：

```text
tss_bed       /share/LCZX_Data/ref/human_hg38_TSS.bed
tss_sha256    cabd65c85ce0db017d771744e9db6ea80b3ef741594763730d6917756914d631
created_at    2026-08-03T...
```

只有 metadata 中的 SHA-256 与当前 TSS BED 一致，`TSS_profile_common.csv` 才允许复用。更换 TSS 文件后必须重新生成 profile，不能继续使用来源不明的旧 `TSS_profile_single_*.csv`。
