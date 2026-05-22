# Response to PR review — compound z-score enrichment

**PR:** #28 `feature/compound-zscore-enrichment`  
**Reviewer:** @adrianomartinelli  
**Date:** 2026-05-22

---

## Comment 1: `derive_library_size` is not the diversity of a library

> *Line R18, `zscore.py`: "as mentioned earlier this is not the diversity of a library."*

**You are correct, and this is the most important issue in the PR.**

### What the reviewer means

The Faver et al. normalized z-score formula is:

```
p_i = 1 / D
z_n = (p_observed - p_i) / sqrt(p_i * (1 - p_i))
```

where **D is the theoretical diversity** — the total number of distinct compounds that were *designed* in the library. This is the combinatorial product of all building blocks across all synthesis cycles:

```
D = n_BB1 × n_BB2 × n_BB3 × ...
```

What the old `derive_library_size()` computed was something different:

```python
# OLD — WRONG as a diversity fallback
def derive_library_size(data, code_columns):
    return int(data[code_columns].drop_duplicates().shape[0])
```

This counts the number of compound identities that actually appeared in the sequencing output — the **observed compound count**. In a real DEL experiment, this is always ≤ D, often substantially less:

| Library | Theoretical D | Typical observed |
|---------|--------------|-----------------|
| Small 2BB (50×50) | 2 500 | ~2 000–2 499 |
| Medium 3BB (50×50×50) | 125 000 | ~60 000–100 000 |
| Large 3BB (200×200×200) | 8 000 000 | ~1 000 000–4 000 000 |

Count files from DELT-Hit demultiplexing contain only compounds with count ≥ 1. Compounds that happened to be entirely undetected in a run are absent from the file entirely. They are still part of the library.

### Why this matters for the z-score

The null hypothesis underpinning the z-score is: *"each of the D designed compounds is equally likely to appear in any given decoded read."* Under this null, the expected frequency of any one compound is `p_i = 1/D`.

Using observed compound count instead of D:
- **Inflates p_i** (because n_observed < D)
- **Deflates z-scores** for enriched compounds (the expected baseline is set too high)
- **Produces non-comparable z-scores** across runs with different sequencing depths (deeper runs observe more compounds → lower fallback diversity → different p_i)

Concretely, for a compound with `count = 50`, `total_reads = 1 000 000`, in a library of `D = 125 000`:

| Diversity used | p_i | z-score |
|---|---|---|
| D = 125 000 (correct) | 8.0 × 10⁻⁶ | 14.86 |
| n_observed = 80 000 (wrong) | 1.25 × 10⁻⁵ | 10.61 |

The z-score is 29% lower when using the wrong denominator.

### The fix applied

**Renamed** `derive_library_size` → `count_observed_compounds` to make clear what it actually computes:

```python
def count_observed_compounds(data: pd.DataFrame, code_columns: list[str]) -> int:
    """Count distinct compound identities (code-column combinations) in a count table.

    This is NOT the library diversity. Do not use this as the z-score diversity
    denominator — use infer_diversity_from_cycle_bbs() or a configured
    library_size instead.
    """
    return int(data[code_columns].drop_duplicates().shape[0])
```

**Added exact config-based diversity resolution.** When `library_size` is not configured, the API now first tries to infer the designed diversity from the DELT-Hit demultiplex `config.yaml`:

```python
def infer_diversity_from_library_config(config: dict, code_columns: list[str] | None = None) -> int:
    building_blocks = config["library"]["building_blocks"]
    diversity = 1
    for building_block in building_blocks:
        diversity *= len(config["whitelists"][building_block])
    return diversity
```

This is exact for DELT-Hit output-style configs because the designed building-block families are represented in `library.building_blocks`, and each family's designed members are represented in `whitelists[B*]`.

The API resolves diversity in this order:

1. explicit `library_size` in the analysis experiment;
2. explicit library config path (`library_config_path`, `library_config`, `demultiplex_config_path`, or `delt_hit_config_path`);
3. auto-discovered parent `config.yaml` from the count-file paths;
4. observed per-cycle BB product as a warning-level fallback.

**Kept** `infer_diversity_from_cycle_bbs()` only as the final fallback when no DELT-Hit library config is available:

```python
def infer_diversity_from_cycle_bbs(data: pd.DataFrame, code_columns: list[str]) -> int:
    """Estimate theoretical library diversity from per-cycle building-block counts.

    Computes n_bb1 × n_bb2 × ... — the Cartesian-product size of the designed
    library and the correct null-hypothesis denominator for the Faver z-score.

    Exact when every designed BB appears in at least one sequenced compound.
    Underestimates when some BBs are entirely absent from the count file.
    Set library_size explicitly in the config whenever the designed diversity is known.
    """
    diversity = 1
    for col in code_columns:
        diversity *= int(data[col].nunique())
    return diversity
```

**Updated** the API path before calling `zscore_analysis()`:

```python
# OLD
observed_library_size = derive_library_size(data, code_columns)
if library_size is None:
    library_size = observed_library_size  # wrong

data = pd.read_csv(data_path)
library_size = resolve_library_size_from_experiment(
    exp=exp,
    analysis_config_path=config_path,
    code_columns=get_code_columns(data),
)
zscore_analysis(
    data=data,
    samples=pd.read_csv(samples_path),
    library_size=library_size,
    save_dir=save_dir / "zscore",
)
```

### Why the config-based whitelist product is the correct default

The `config.yaml` whitelist product is exact because it is the designed library definition, not a sequencing-derived observation. For a DELT-Hit config:

```
D = len(whitelists["B0"]) × len(whitelists["B1"]) × ...
```

The observed per-cycle product is still better than the observed compound count, but it remains an approximation. It is used only when no designed library config is available.

### When explicit `library_size` is still useful

Set `library_size` explicitly when the analysis count files do not live under the original DELT-Hit demultiplex output folder, when multiple library configs are nearby, or when the designed diversity differs from the raw whitelist product because of a deliberate downstream library filter:

```yaml
experiments:
  - name: my_experiment
    library_size: 125000   # 50 × 50 × 50 designed BBs
    save_dir: ...
    selections: ...
```

### Validation check preserved

The check `library_size >= n_observed` is preserved (now using `count_observed_compounds` for `n_observed`). This catches configurations where the user accidentally provides a `library_size` smaller than the number of distinct compounds in the data — an impossible situation that indicates a configuration error.

---

## Comment 2: Why compute fold enrichment?

> *Line R55, `zscore.py`: "why do we compute a fold enrichment here?"*

**Fold enrichment is included as a companion metric for hit picking and reporting.**

### What fold enrichment is

```
expected_count  = total_reads / library_size
fold_enrichment = count / expected_count
```

It answers: "how many times more often did this compound appear than expected by chance?" A compound at `fold_enrichment = 10` was sequenced 10× more often than a random library member.

### Why it is included alongside z-score

z-score and fold enrichment are monotonically related — ranking by one is identical to ranking by the other. They differ in what they express:

| | z-score | fold enrichment |
|---|---|---|
| Scale | standard deviations | multiplicative factor |
| Value at expectation | 0 | 1 |
| Interpretable to chemists | no — needs statistical vocabulary | yes — "10× enriched" |
| Accounts for sequencing depth | yes — bakes in sampling variance | no — pure ratio |
| Good for ranking | yes | yes (identical rank order) |
| Good for reporting | no | yes |
| Good for thresholding | no (z ≥ 1.0 is not a meaningful cutoff) | yes (FE ≥ 2 or 3 is interpretable) |

Including both in the per-selection CSV allows users to:
1. **Rank** by z-score (the statistically appropriate primary signal)
2. **Report** by fold enrichment in papers, presentations, and to collaborators ("compound X was 15-fold enriched over the no-protein control")
3. **Apply a magnitude gate** alongside the z-score: `zscore > 2 AND fold_enrichment >= 3` filters out statistical hits with negligible signal intensity

This two-metric convention is standard in the DEL literature and is explicitly described in the DELT-Hit manuscript (statistical methods section).

### Precedent in existing DELT-Hit output

The edgeR method outputs `logFC`, which is `log2(fold_enrichment)`. Fold enrichment is the natural-scale version of the same quantity. Including it here is consistent with the information already provided by the R-based methods.

### Could it be removed?

If the preference is to keep per-selection CSVs minimal, fold enrichment could be dropped and computed post-hoc by the user from `count` and `expected_count` (which are both present in the CSV). However, given that:
- It requires no additional computation (it is computed as an intermediate anyway)
- It is the single most interpretable metric for hit picking
- It matches the manuscript description

the recommendation is to keep it.

---

## Summary of code changes made

| File | Change |
|------|--------|
| `src/delt_hit/cli/analyse/zscore.py` | Replaced `derive_library_size` with `count_observed_compounds`; added exact config-based diversity inference from `library.building_blocks` and `whitelists`; retained observed per-cycle inference only as a final warning-level fallback; clarified validation errors |
| `src/delt_hit/cli/analyse/api.py` | Resolves `library_size` before z-score analysis by using explicit `library_size`, explicit library config path, auto-discovered parent `config.yaml`, then fallback inference |
| `tests/analyse/test_zscore_synthetic_libraries.py` | Added synthetic DELT-Hit-style YAML/count-file tests for uneven 1-, 2-, 3-, and 4-cycle libraries, sparse observed coverage, exact config diversity, explicit config path resolution, and API auto-discovery |

All 58 relevant analyse/CLI tests pass.

---

## Recommended additional action

Add `library_size` to example analysis configs when the count files are detached from the original demultiplex `config.yaml`; otherwise keep the auto-discovery behavior documented so standard DELT-Hit output folders work without duplicate metadata.
