# Compound-level z-score enrichment — implementation reference

**Branch:** `feature/compound-zscore-enrichment`  
**Commit:** `2c6f232`  
**Date:** 2026-05-21  
**Status:** ready for pull request

---

## Overview

This document describes the compound-level normalized z-score enrichment added to DELT-Hit. It covers the statistical formula, CLI usage, all inputs and outputs, design choices, and the validation logic that guards against common input errors.

The z-score method is **pure Python** — unlike the `counts` and `edgeR` methods, it does not generate an R script. Results are written immediately as CSV files.

---

## Statistical background

### Normalized z-score (Faver et al. 2019)

The enrichment score for one compound in one selection is:

```
p_i        = 1 / library_size          # expected frequency under uniform sampling
p_observed = count / total_reads       # observed frequency
z_n        = (p_observed - p_i) / sqrt(p_i * (1 - p_i))
```

`z_n > 0` means the compound appears more often than expected by chance.  
`z_n = 0` means it appears exactly at the expected rate.  
`z_n < 0` means it appears less often than expected (depleted).

Under uniform sampling, `z_n` is approximately standard-normally distributed, so values > 2–3 indicate strong enrichment. However, **fixed z-score thresholds should not be used as a binary hit/no-hit gate**; fold enrichment (below) is a more interpretable companion metric for hit picking.

### Fold enrichment

```
expected_count    = total_reads / library_size
fold_enrichment   = count / expected_count
```

`fold_enrichment = 1` means observed equals expected.  
`fold_enrichment = 10` means the compound appeared ten times more often than expected.  
A fold enrichment ≥ 2 combined with a z-score > 0 is a reasonable starting threshold for hit picking.

### Confidence intervals

Confidence intervals on `z_n` are derived from Agresti-Coull intervals on `p_observed`, then transformed back to z-score units:

```
z_α      = scipy.stats.norm.ppf(1 - alpha/2)   # e.g. 1.96 for alpha=0.05
n'       = total_reads + z_α²
p_adj    = (count + z_α²/2) / n'
se_adj   = sqrt(p_adj * (1 - p_adj) / n')

ci_p_lower = p_adj - z_α * se_adj
ci_p_upper = p_adj + z_α * se_adj

zscore_ci_lower = (ci_p_lower - p_i) / sqrt(p_i * (1 - p_i))
zscore_ci_upper = (ci_p_upper - p_i) / sqrt(p_i * (1 - p_i))
```

Agresti-Coull was chosen because it has good coverage properties even at low counts (unlike the plain Wilson interval) and does not require count > 0.

---

## CLI usage

### Command

```bash
delt-hit analyse enrichment \
  --config_path analysis.yaml \
  --name experiment_name \
  --method zscore
```

### Configuration file (`analysis.yaml`)

```yaml
experiments:
  - name: protein_vs_no_protein
    save_dir: /path/to/results
    library_size: 122247          # optional when a DELT-Hit config.yaml is discoverable
    selections:
      - name: protein_1
        counts_path: /path/to/protein_1.txt
        group: protein
      - name: protein_2
        counts_path: /path/to/protein_2.txt
        group: protein
      - name: control_1
        counts_path: /path/to/control_1.txt
        group: no_protein
      - name: control_2
        counts_path: /path/to/control_2.txt
        group: no_protein
```

**`library_size`** is the total number of distinct compounds in the designed library. Setting it explicitly is always accepted. If it is omitted, DELT-Hit first tries to infer the exact designed diversity from the demultiplex `config.yaml` (`library.building_blocks` and the corresponding `whitelists`). The config can be given explicitly with `library_config_path`, or auto-discovered as a parent `config.yaml` of the configured count files. Only if no library config is available does the implementation fall back to the product of observed per-cycle code values, which is a warning-level approximation.

**`group`** names are free-form strings. The two canonical group names are `protein` and `no_protein`, which enable the protein–control contrast column in `stats.csv`. Any other group name (e.g. `naive`) is included in `stats.csv` as an extra column, and the `enrichment` column is always `protein - no_protein` (set to 0.0 if either is absent).

### Input count files

Each selection points to a tab-separated counts file produced by DELT-Hit demultiplexing. Required columns:

| Column | Description |
|--------|-------------|
| `code_1`, `code_2`, ... | Building block identifiers for each cycle |
| `count` | Read count for this compound identity |

Optional columns (e.g. `id`) are read silently and ignored by the z-score analysis. Any column whose name starts with `code_` is treated as a building-block identifier. The number of cycles is detected automatically from how many `code_*` columns are present — two-cycle, three-cycle, and four-cycle libraries are all supported without configuration.

---

## Processing pipeline

### Step 1: Data preparation (`Analyse.prepare`)

`Analyse.prepare()` is called before enrichment runs regardless of the method chosen:

1. Reads the YAML config.
2. Loads each selection's count file, attaches the selection `name`.
3. Concatenates all selections into one tall DataFrame.
4. Writes `<save_dir>/<experiment>/data.csv` and `<save_dir>/<experiment>/samples.csv`.

### Step 2: `zscore_analysis()` — orchestrator

#### Validation

| Check | Error raised if violated |
|-------|--------------------------|
| At least one `code_*` column in data | `ValueError: z-score analysis requires at least one code_* column` |
| `name` and `count` columns present in data | `ValueError: Count data is missing required columns` |
| `name` and `group` columns present in samples | `ValueError: Sample metadata is missing required columns` |
| All data selection names appear in samples | `ValueError: Count data contains selections absent from samples metadata` |
| No duplicate selection names in samples | `ValueError: Duplicate sample names in metadata` |
| No missing count values | `ValueError: Count data contains missing count values` |
| No negative count values | `ValueError: Count data contains negative count values` |
| No missing code identifiers | `ValueError: Count data contains missing code_* values` |
| `library_size` is an integer ≥ 2 | `ValueError: library_size must be...` |
| Configured `library_size` ≥ observed compound union | `ValueError: Configured library_size is smaller than...` |

#### Compound grid completion (zero-fill)

After validation, the code builds a full compound × selection grid:

1. Collect all distinct compound identities (all unique code-column combinations) observed across all selections.
2. Take the Cartesian product with the selection list.
3. Left-join the actual counts onto this grid; fill any missing combinations with `count = 0`.

This zero-fill is **intentional and important**. Without it, a compound seen in the protein selection but not in the control would be absent from the control's score table, making cross-selection comparison impossible. Zero-filling assigns the minimally enriched z-score to genuinely undetected compounds.

#### Scoring (`compute_zscore`)

For each selection separately:

1. Compute `total_reads` = sum of all counts in that selection.
2. Raise `ValueError` if `total_reads = 0`.
3. Compute `p_i`, `p_observed`, `z_n`, CIs, and `fold_enrichment` using the vectorized NumPy functions in `delt_hit.analyse.enrichment`.
4. Attach all score columns to the selection DataFrame.

#### Group statistics (`_aggregate_group_stats`)

1. Group by `(code_*, group)` and average the z-scores across replicates in the same group.
2. Pivot to wide format: one column per group name.
3. Ensure `protein` and `no_protein` columns exist (fill with 0.0 if a group is absent).
4. Compute `enrichment = protein - no_protein`.
5. Column order: `code_*`, `protein`, `no_protein`, other groups, `enrichment`.

#### Hit file

`hits.csv` is the top-100 rows of `stats.csv`, sorted descending by:
- `protein` z-score, if a `protein` group is present in the samples; or
- `enrichment` otherwise.

The `hits_top_n` parameter (default 100) controls the number of rows written.

---

## Output files

All outputs are written under `<save_dir>/<experiment_name>/zscore/`.

### Per-selection files: `<selection_name>.csv`

One file per selection. Columns:

| Column | Description |
|--------|-------------|
| `code_1`, `code_2`, ... | Building block identifiers |
| `count` | Raw read count (0-filled for absent compounds) |
| `observed_fraction` | `count / total_reads` |
| `expected_fraction` | `1 / library_size` (same for every row) |
| `expected_count` | `total_reads / library_size` |
| `fold_enrichment` | `count / expected_count` |
| `zscore` | Faver normalized z-score |
| `zscore_ci_lower` | Lower bound of Agresti-Coull CI (default 95%) |
| `zscore_ci_upper` | Upper bound of Agresti-Coull CI (default 95%) |

Sorted descending by `zscore` within each selection.

### Aggregate file: `stats.csv`

One row per compound (across all selections). Columns:

| Column | Description |
|--------|-------------|
| `code_1`, `code_2`, ... | Building block identifiers |
| `protein` | Mean z-score across all `protein`-group replicates |
| `no_protein` | Mean z-score across all `no_protein`-group replicates |
| *(other groups)* | Mean z-score for any additional group present |
| `enrichment` | `protein - no_protein` |

If the experiment has no `protein` group, the column is added as 0.0. Same for `no_protein`.

### Hit file: `hits.csv`

Top 100 rows of `stats.csv` (same columns), sorted by `protein` z-score if that group exists, otherwise by `enrichment`.

---

## Code structure

```
src/delt_hit/
  analyse/
    __init__.py              package marker
    enrichment.py            math primitives (calculate_zscore, calculate_fold_enrichment)
  cli/
    analyse/
      zscore.py              CLI orchestration (zscore_analysis, compute_zscore, helpers)
      api.py                 Analyse class — wires 'zscore' case into enrichment()

tests/
  analyse/
    __init__.py
    test_zscore.py           unit tests for enrichment.py math
    test_cli_zscore_analysis.py  integration tests for zscore_analysis and compute_zscore
  cli/
    __init__.py
    test_analyse_zscore_cli.py   subprocess CLI smoke test
```

### Why a separate `enrichment.py` module?

The formula is placed in `delt_hit.analyse.enrichment` so it can be tested and imported independently of the CLI. The CLI wrapper in `zscore.py` handles data-loading concerns (validation, zero-fill, grouping, file I/O) while `enrichment.py` is kept to pure math. This also makes the formula reusable for future multi-level synthon scoring without depending on any CLI infrastructure.

---

## Design decisions

### `library_size` vs observed diversity

Using the theoretical library size as the denominator is the scientifically correct choice: it accounts for compounds that were designed but never appeared in any selection. Deriving it from observed data underestimates the diversity denominator for large sparse DEL libraries, which inflates z-scores for genuinely random noise.

Resolution order is:

1. explicit `library_size` in `analysis.yaml`;
2. exact whitelist product from `library_config_path` / `library_config` / `demultiplex_config_path` / `delt_hit_config_path`;
3. exact whitelist product from an auto-discovered parent DELT-Hit `config.yaml` next to the count files;
4. approximate product of observed per-cycle code values, with a warning.

A configured or inferred `library_size` smaller than the observed compound count raises `ValueError`, because that is impossible for a full combinatorial library and almost certainly means the wrong config was supplied.

### Zero-filling absent compounds

Compounds observed in one selection but absent from another must be zero-filled rather than excluded. The Faver z-score assigns a large negative value to `count = 0`, correctly reflecting that a compound was undetected. Excluding them would silently drop true depletions and make the hit list asymmetric.

### Agresti-Coull confidence intervals

Agresti-Coull is preferred over the simple Wald interval because Wald CI coverage fails at low counts (a common situation in DEL data, where most compounds appear 0–5 times). Agresti-Coull provides near-nominal coverage even at `count = 0`.

### Group naming convention

The `protein` and `no_protein` group names are conventions, not requirements. The implementation detects them by string matching and falls back gracefully if either is absent. This allows the CLI to handle experiments with non-standard group names (e.g. a single target group with no explicit control) without error.

### `hits.csv` sort order

When a `protein` group is present, `hits.csv` is sorted by `protein` z-score rather than `enrichment`. This prioritizes compounds that are strongly enriched in the target selection regardless of control behavior — important in experiments where the control is noisy or absent. When no `protein` group exists, `enrichment` (which defaults to `protein - no_protein = 0 - no_protein = -no_protein`) captures the next-best proxy.

### Double config read in `api.py`

The `zscore` case re-reads the YAML file to resolve `library_size` and optional library-config paths. This is a minor redundancy (the config was already read in `prepare()`), but `prepare()` does not currently return the raw config dict, so the local re-read avoids a larger refactor of the `Analyse` API. At typical config sizes the cost is negligible.

### `assert` statements in `prepare()` and `prepare_data()`

`Analyse.prepare()` and `prepare_data()` use `assert` statements for validation. These predate the z-score work and are not changed here. Production code should use explicit `ValueError` or `FileNotFoundError` raises — these `assert`s would silently pass under Python's `-O` optimized flag. This is a **pre-existing issue** not introduced by the z-score PR.

---

## Test coverage

### `tests/analyse/test_zscore.py` — math unit tests

| Test | What it checks |
|------|----------------|
| `test_zscore_canonical_value` | Exact numerical result against hand-calculated formula |
| `test_zscore_zero_when_count_equals_expected` | z = 0 at the expected count |
| `test_zscore_negative_when_below_expected` | z < 0 for depleted compound |
| `test_zscore_positive_when_above_expected` | z > 0 for enriched compound |
| `test_ci_ordering_positive_counts` | CI lower < z < CI upper for enriched compounds |
| `test_ci_ordering_depleted_counts` | CI ordering holds for zero/depleted counts |
| `test_ci_width_shrinks_with_more_reads` | CI narrows with more sequencing depth |
| `test_zscore_count_zero` | z is large-negative and matches formula for count=0 |
| `test_zscore_count_equals_total_reads` | Extreme enrichment case (all reads = one compound) |
| `test_zscore_diversity_one_returns_nan` | NaN returned when diversity=1 (denominator=0) |
| `test_zscore_array_shape_preserved` | Output shape matches input shape |
| `test_zscore_matches_scalar_formula_elementwise` | Vectorized output = scalar formula for each element |
| `test_zscore_raises_on_zero_total_reads` | `ValueError` on invalid total_reads |
| `test_zscore_raises_on_zero_diversity` | `ValueError` on invalid diversity |
| `test_fold_enrichment_canonical` | Expected = 10× enrichment at 10× count |
| `test_fold_enrichment_at_expected_is_one` | fold_enrichment=1 at the expected count |
| `test_fold_enrichment_zero_count` | fold_enrichment=0 for count=0 |

### `tests/analyse/test_cli_zscore_analysis.py` — integration tests

| Test | What it checks |
|------|----------------|
| `test_zscore_analysis_writes_protocol_outputs_and_zero_fills` | All output files created; zero-fill works; expected_fraction, expected_count, fold_enrichment correct; stats columns and hit ordering |
| `test_zscore_analysis_supports_three_cycle_code_columns` | Three `code_*` columns handled correctly |
| `test_compute_zscore_rejects_zero_total_selection` | Raises on selection with all-zero counts |
| `test_zscore_analysis_rejects_missing_code_columns` | Raises when no `code_*` column present |
| `test_zscore_analysis_rejects_missing_code_values` | Raises on NaN in code columns |
| `test_zscore_analysis_rejects_library_size_smaller_than_observed` | Raises when configured size < observed |
| `test_analyse_enrichment_zscore_wires_cli_method` | Full Analyse.enrichment() path writes expected outputs |
| `test_analyse_enrichment_rejects_unknown_method` | Raises `ValueError` on unknown method string |

### `tests/cli/test_analyse_zscore_cli.py` — subprocess smoke test

| Test | What it checks |
|------|----------------|
| `test_delt_hit_analyse_enrichment_zscore_cli` | Full CLI invocation via subprocess returns exit code 0 and writes stats.csv and hits.csv |

**All 26 tests pass** on Python 3.12 in the pixi test environment.

---

## Known limitations and future work

| Limitation | Notes |
|------------|-------|
| Compound-level only | Multi-level synthon z-score (per-building-block) is not wired into this CLI pass |
| No correlation PNGs | The `counts` and `edgeR` methods generate R/GGally correlation plots; the z-score method does not (it matches the manuscript protocol) |
| `library_size` must be integer | Float or None with auto-derive is supported; anything else raises |
| `hits_top_n` not exposed to CLI | Fixed at 100; adjustable via the Python API only |
| `alpha` not exposed to CLI | Fixed at 0.05 (95% CI); adjustable via the Python API only |
| Pre-existing `assert` in `prepare()` | `prepare()` and `prepare_data()` use `assert` for validation; these silently pass under `-O` |

---

## Real-data validation

The CLI was verified on Dimitar 3BB Thrombin SP data (122 247 distinct 3-cycle compounds, 6 selections — 3 protein replicates × 3 no-protein replicates):

- `stats.csv`: 122 247 rows, columns `code_1, code_2, code_3, protein, no_protein, enrichment`
- `hits.csv`: 100 rows, sorted by `protein` z-score
- Per-selection files: 122 247 rows each
- Known Thrombin hits (codes 23,44,8 and 23,44,14) ranked at positions 1 and 2 by protein z-score

This confirms that the implementation produces protocol-shaped outputs and recovers established positive controls.
