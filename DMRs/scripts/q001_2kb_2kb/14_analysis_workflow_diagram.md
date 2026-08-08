# q<0.01 分析流程示意图

```mermaid
flowchart TB
    A[研究问题：肺癌免疫治疗<br/>IR 应答者 vs NR 不应答者] --> B[输入数据]

    B --> B1[DNA 单细胞甲基化数据<br/>All VMR matrix + Meth_diff DMR]
    B --> B2[RNA 单细胞转录组数据<br/>ALL_batch_corrected_pbmc.h5ad]
    B --> B3[细胞分组文件<br/>sample、IR/NR、cell type]

    %% shared metadata
    B3 --> S00[00 生成 DNA cell metadata<br/>cell_id → sample + response + cell type]

    %% clean cell type branch
    subgraph C[分支 A：clean cell-type DMR / VMR（细胞类型特征）]
      direction TB
      C1[细胞类型 pairwise DMRs<br/>q<0.01] --> C3[01 合并 DMR]
      C2[same-cell-type IR vs NR<br/>sample/response component DMRs<br/>q<0.01] --> C3
      C3 --> C4[从 cell-type DMRs 中<br/>扣除 response component]
      C4 --> C5[clean cell-type DMRs]
      C5 --> C6[映射至 All VMR regions]
      C6 --> C7[02–03 提取对应 VMR matrix<br/>25,706 个 VMR features]
    end

    %% disease response branch
    subgraph D[分支 B：IR vs NR 应答相关 DMR（主分析）]
      direction TB
      D1[same-cell-type IR vs NR<br/>全部 Meth_diff DMRs] --> D2[04 汇总每种细胞类型]
      D2 --> D3[05 q<0.01 筛选并标注方向]
      D3 --> D4{Meth_diff 低甲基化组}
      D4 -->|NR 低甲基化| D5[IR-hyper DMR]
      D4 -->|IR 低甲基化| D6[IR-hypo DMR]
    end

    %% promoter and rna
    D5 --> P1
    D6 --> P1
    subgraph P[Promoter 注释与候选基因]
      direction TB
      P0[06 GENCODE v44 promoter<br/>TSS ±2 kb，共享注释] --> P1[07 DMR ∩ promoter]
      P1 --> P2[08 DMR-to-gene 表]
      P2 --> P3[09 保留 protein-coding genes<br/>并按 cell type × direction × gene 去重]
    end

    %% rna integration
    P3 --> R1
    B2 --> R1
    subgraph R[RNA pseudobulk 整合]
      direction TB
      R1[10–11 对每个 sample × cell type<br/>计算候选基因平均 RNA 表达] --> R2[计算 IR − NR expression delta]
      R2 --> R3{方向是否符合 promoter<br/>甲基化-表达负向模式？}
      R3 -->|IR-hyper + RNA down| R4[反向方向候选]
      R3 -->|IR-hypo + RNA up| R4
      R3 -->|其他| R5[其他 DMR-gene 候选]
    end

    %% correlation
    S00 --> K1
    P3 --> K1
    B1 --> K1
    R1 --> K1
    subgraph K[12 DNA-RNA 样本配对相关性]
      direction TB
      K1[将 promoter DMR 映射回 VMR feature] --> K2[DNA pseudobulk：同 sample × cell type<br/>sum methylated sites / sum total sites]
      K2 --> K3[与同 sample × cell type 的 RNA 平均表达配对]
      K3 --> K4[Pearson / Spearman + BH-FDR]
      K4 --> K5[pooled IR+NR：描述性相关性<br/>IR/NR 内部：需至少 6 个配对样本]
    end

    R4 --> O[输出：高优先级探索性候选]
    K5 --> O
    O --> O1[严格 q<0.01：CD14 Monocytes 为主要信号<br/>ZFP57：跨 q<0.01 / q<0.05 的候选]

    classDef input fill:#e8f1ff,stroke:#4a78b7,color:#102a43;
    classDef branch fill:#f4f4f5,stroke:#666,color:#222;
    classDef key fill:#e8f7ee,stroke:#3f8d5a,color:#143d24;
    classDef caution fill:#fff6df,stroke:#b7791f,color:#513b08;
    class A,B,B1,B2,B3 input;
    class C,D,P,R,K branch;
    class C5,D5,D6,P3,R4,O,O1 key;
    class K5 caution;
```

## 如何阅读这张图

- **分支 A（00–03）**：构建 clean cell-type DMR/VMR 特征；先扣除同一细胞类型中 IR 与 NR 的 component，避免将应答状态当作细胞类型差异。
- **分支 B（04–12）**：直接在相同细胞类型内比较 IR 与 NR，保留的正是与免疫治疗应答状态相关的 DMR；这是 promoter、RNA 和相关性整合所使用的主线。
- **关键区分**：clean 分支的“扣除”不能用于分支 B；若在疾病/应答 DMR 上再扣除 IR-vs-NR 信号，会把研究对象本身删掉。
- **相关性限制**：目前每个 response × cell type 只有 5 个样本。`min_paired_samples=6` 时只计算合并 IR+NR 的 pooled 相关性；它用于描述性优先级排序，不能作为组内关联或因果结论。
