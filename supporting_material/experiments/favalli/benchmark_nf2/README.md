# NF2 Enrichment Benchmark

This directory contains the curated Favalli NF-DEL benchmark used in
Supplementary Note 4.

## Contents

- `analysis.yaml`: Eight benchmark conditions for CAIX, CREBBP, wt-PI3K, and
  H1047R-PI3K. The failed CAIX ss count table `NF2_NatChem_fl1_1` is excluded
  because it contains only 488 total reads versus 335k-564k reads for the valid
  CAIX ss protein selections.
- `truth_set.csv`: Curated positive-control and contextual rows. The headline
  benchmark uses only `label == strict_affinity`.
- `strict14_method_ranks.csv`: Per-instance ranks for the six methods reported
  in the supplementary information.
- `zscore_outputs/benchmark_summary.csv`: Native DELT-Hit z-score ranks for the
  strict-affinity controls.
- `rerun_counts_edger_base_r.R`: Recomputes DELT-Hit simple-count and edgeR
  tables from generated `outputs/<experiment>/data.csv` and `samples.csv`
  without requiring tidyverse or GGally.

## Rerunning

From the repository root, rerun the native DELT-Hit z-score benchmark with:

```bash
pixi run python supporting_material/experiments/favalli/benchmark_zscore.py
```

This regenerates the `outputs/<experiment>/zscore/` directories and
`zscore_outputs/benchmark_summary.csv`.

After the z-score command has generated `outputs/<experiment>/data.csv` and
`samples.csv`, rerun simple-count and edgeR outputs for one condition with:

```bash
cd supporting_material/experiments/favalli/benchmark_nf2
Rscript rerun_counts_edger_base_r.R outputs/ss_caix
```

The compatibility script uses base R for I/O and reshaping, and uses the same
edgeR/limma model as the generated DELT-Hit edgeR script. It does not create the
optional GGally correlation plots.

## DELi Rank Provenance

The DELi ranks in `strict14_method_ranks.csv` were generated from DELi cube
CSVs and configs prepared in the DELexplore reviewer-response workflow
(`scripts/prepare_deli_nf2.py`, `scripts/run_deli_nf2.py`,
`scripts/reviewer_response/prepare_deli_fl2_pi3k.py`, and
`scripts/reviewer_response/extract_deli_pi3k_ranks.py`). The CAIX ss DELi
input was rechecked after excluding `NF2_NatChem_fl1_1` and rebuilding the cube
from the remaining two protein selections plus the three ss no-protein controls.
This check leaves the DELi NSC and MLE ranks unchanged for the two strict CAIX
ss positives and changes only the DELi norm-z per-instance ranks
(`71244 -> 71210` and `71245 -> 71211`). The summary metrics reported in the
supplementary information are unchanged.
