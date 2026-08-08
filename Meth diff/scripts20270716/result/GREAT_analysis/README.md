# GREAT analysis handoff

## Status

The four hg38 input sets were prepared, audited, and analyzed separately on
2026-07-20. The Stanford GREAT web application was unavailable (expired TLS
certificate and HTTP 404 on the former public routes), so the analysis was
completed locally with `rGREAT`.

The local run reproduces the requested core model:

- whole-genome background;
- basal-plus-extension association;
- 5 kb upstream and 1 kb downstream basal domains;
- maximum extension of 1,000 kb;
- hg38 UCSC known-gene TSS annotation;
- UCSC assembly gaps excluded.

Important difference: Stanford's curated regulatory domains are not bundled
with local `rGREAT`, so they were not added. The local output covers GO
Biological Process, Cellular Component, and Molecular Function. When the web
service is restored, the four BED files should still be submitted separately
for a strict web-GREAT replication and export of all server ontologies.

## Prepared inputs

All inputs are headerless BED4 files. They contain canonical chromosomes only,
have unique region IDs, satisfy source MethSCAn adjusted p-value below 0.05,
and overlap the corresponding post-sample-effect clean BED.

| Input | Regions |
|---|---:|
| `GREAT_inputs/CD14_Monocytes_IR_hypo_q005_clean.bed` | 849 |
| `GREAT_inputs/CD14_Monocytes_IR_hyper_q005_clean.bed` | 563 |
| `GREAT_inputs/NK_cells_IR_hypo_q005_clean.bed` | 9,938 |
| `GREAT_inputs/NK_cells_IR_hyper_q005_clean.bed` | 103 |

`group_A` was defined as IR and `group_B` as NR. MethSCAn column 10 records
the lower-methylated group, so `group_A` maps to `IR_hypo`, while `group_B`
maps to `IR_hyper`.

## Cell counts and MethSCAn DMR counts

Cell counts were taken from
`ALL_cell_count_by_response_cell_type_200k.csv`. DMR counts below include only
`chr1`–`chr22`, `chrX`, and `chrY`. `Raw P < 0.05` uses MethSCAn column 11;
`adjusted P < 0.05` uses the permutation-derived adjusted P in column 12.
IR-hypo and IR-hyper are direction splits within the adjusted-P set.

| Cell type | IR cells | NR cells | Total cells | Raw P < 0.05 | Adjusted P < 0.05 | IR-hypo | IR-hyper |
|---|---:|---:|---:|---:|---:|---:|---:|
| CD14 Monocytes | 2,811 | 2,212 | 5,023 | 87,976 | 1,984 | 1,165 | 819 |
| NK cells | 737 | 745 | 1,482 | 60,120 | 10,045 | 9,941 | 104 |
| CD4 T cells | 1,829 | 1,263 | 3,092 | 89,043 | 26 | 16 | 10 |
| Cycling cells | 67 | 81 | 148 | 22,193 | 3 | 0 | 3 |
| MAIT cells | 128 | 129 | 257 | 60,846 | 3 | 0 | 3 |
| CD16 Monocytes | 693 | 555 | 1,248 | 86,385 | 2 | 1 | 1 |
| CD8 T cells | 708 | 1,154 | 1,862 | 93,466 | 2 | 1 | 1 |
| Treg cells | 119 | 126 | 245 | 44,151 | 2 | 0 | 2 |
| Plasma cells | 105 | 124 | 229 | 52,652 | 3 | 2 | 1 |
| B cells | 543 | 642 | 1,185 | 92,944 | 1 | 1 | 0 |
| HLAII-high APCs | 296 | 98 | 394 | 38,669 | 1 | 0 | 1 |
| B cells unresolved | 25 | 33 | 58 | 6,401 | 0 | 0 | 0 |
| Gamma-delta T cells | 57 | 36 | 93 | 7,471 | 0 | 0 | 0 |
| cDCs | 13 | 9 | 22 | 4 | 0 | 0 | 0 |
| pDCs | 96 | 32 | 128 | 4,649 | 0 | 0 | 0 |

The annotation additionally contains 21 IR and 13 NR cells labeled
`Platelet_erythroid_contamination` (34 total); this label had no corresponding
IR-vs-NR DMR output and is therefore excluded from the DMR table. Across all
annotation labels, there are 8,248 IR and 7,252 NR cells (15,500 total).
The per-DMR `n_cells_group1` and `n_cells_group2` columns are region-specific
coverage counts and should not be confused with the total cell counts above.

## Result overview

Primary criteria were binomial FDR below 0.05, region fold enrichment at least
2, at least 5 observed region hits, and at least 3 observed gene hits.

| Analysis | Regions | Associated genes | Primary terms | Significant by both tests | Concentration warnings |
|---|---:|---:|---:|---:|---:|
| CD14 IR-hypo | 849 | 1,308 | 8 | 3 | 3 |
| CD14 IR-hyper | 563 | 802 | 1 | 0 | 1 |
| NK IR-hypo | 9,938 | 6,571 | 125 | 11 | 48 |
| NK IR-hyper | 103 | 178 | 0 | 0 | 0 |

No input set had fewer than 20 regions. CD14 IR-hyper had one region-based
term (`mitochondrial crista`) but it failed the gene-based FDR threshold and
was driven mainly by OPA1-associated regions. NK IR-hyper had no term meeting
the primary criteria.

The 14 terms significant in both tests collapse to 9 descriptive themes in
`GREAT_results/theme_summary.tsv`. Most are developmental, neural, cardiac, or
tissue-repair terms rather than canonical CD14-monocyte or NK-cell immune
programs. They must therefore be treated as exploratory signals, not stable
immune mechanisms.

## Output guide

Top-level summaries:

- `GREAT_results/run_summary.tsv`: counts per analysis.
- `GREAT_results/primary_summary.tsv`: all 134 terms meeting primary criteria,
  including FDRs, fold enrichment, region/gene hits, hit genes, and locus
  concentration diagnostics.
- `GREAT_results/significant_by_both_summary.tsv`: the 14 terms passing both
  binomial and hypergeometric FDR below 0.05, with theme annotations.
- `GREAT_results/theme_summary.tsv`: redundant double-significant terms grouped
  into themes with representative terms.
- `GREAT_results/hit_count_translation_audit.tsv`: 21 small discrepancies found
  while translating Entrez gene sets to symbol-based region associations.
  Filtering always uses the authoritative `rGREAT` observed-region count.

Each analysis directory contains:

- `GO_all_terms.tsv`;
- `GO_primary_filtered.tsv`;
- `GO_significant_by_both.tsv`;
- `GO_primary_term_hits.tsv`;
- `region_gene_associations.tsv`;
- `gene_region_associations.tsv`;
- `tss_distance_summary.tsv`;
- `analysis_parameters.tsv`.

Input-audit files:

- `GREAT_inputs/preparation_summary.tsv` records every filtering stage.
- `GREAT_inputs/*_selected_source_rows.tsv` preserves the selected source rows
  and parsed adjusted p-values.

## Interpretation limits

The current MethSCAn comparison used cells as statistical units. In addition,
the NK IR-hypo set is very large and 48 of its 125 primary terms have at least
half of region hits concentrated at one associated gene locus. Only terms
passing both tests and surviving the concentration review should be discussed,
and even those require donor-level confirmation.

For a paper-level analysis, construct sample-by-cell-type pseudobulks, call DMRs
with DSS or bsseq, apply FDR below 0.05, absolute methylation difference at
least 0.10, at least 3 CpGs, and valid coverage in at least 3 of 5 samples per
group. Repeat the separated hyper/hypo GREAT analysis and retain only
directionally and functionally concordant findings.
