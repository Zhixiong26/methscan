# Legacy scripts

These files implement the retired post-hoc top 5%/2%/1% matrix-column
subsetting workflow. They are retained only for reproducibility and are not
called by the current threshold-specific pipeline.

Current entry points are in the parent directory:

```text
10_prepare_individual_effect_mask.sh
11_run_threshold_vmrs_remove_individual.sh
13_run_threshold_clean_vmr_reclustering.sh
14_collect_threshold_metrics.sh
```
