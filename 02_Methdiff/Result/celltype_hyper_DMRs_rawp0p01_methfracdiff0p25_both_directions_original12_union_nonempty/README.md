Cell-type hypermethylated DMR extraction
============================================

Selection
---------
* MethSCAn column 11 raw p < 0.01
* abs(column 8 - column 9) >= 0.25
* Direction follows column 10 (lower-methylated group); the opposite group is hyper.
* Primary chromosomes only: True

Outputs
-------
* pairwise/: hyper-DMRs for each target-cell versus other-cell comparison;
  rows retain the input BED's original 12 columns exactly, with no header added.
* pairwise_union/: one three-column BED per pairwise file, containing the union
  of overlapping intervals. Book-ended but non-overlapping intervals stay separate.
* by_cell_type/*__hyper_records.tsv: all standardized, unmerged pairwise records.
* by_cell_type/*__hyper_any_other.bed: merged union; supported by >=1 other cell type.
* by_cell_type/*__hyper_all_others.bed: strict genomic segments simultaneously
  covered by hyper-DMRs against every available other cell type.
* sample_summary.tsv and overall_summary.tsv: record and region counts.

Important
---------
These are DMR intervals, not individual CpG coordinates. To obtain individual CpG
sites, intersect these BED intervals with the original CpG-level methylation matrix.
Raw p selection does not control the false-discovery rate.
