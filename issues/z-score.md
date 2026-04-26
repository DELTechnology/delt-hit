# Normalized Z-Score Support for DELT

## Summary

This note documents the normalized z-score method proposed by Faver et al. for DNA-encoded library (DEL) enrichment analysis and outlines the smallest DELT change that would add it as a third enrichment mode for experiments without replicates.

Primary paper:

- John C. Faver et al., "Quantitative Comparison of Enrichment from DNA-Encoded Chemical Library Selections", *ACS Combinatorial Science* 2019, 21(2), 75-82.
- DOI: [10.1021/acscombsci.8b00116](https://doi.org/10.1021/acscombsci.8b00116)
- PubMed: [30672692](https://pubmed.ncbi.nlm.nih.gov/30672692/)

## What the Paper Proposes

Faver et al. propose a normalized z-score enrichment metric for DEL sequencing data using a binomial sampling model. Their motivation is that raw counts and fold-enrichment alone are hard to compare across selections because enrichment measurements are affected by sequencing depth, library diversity, expected barcode frequency, and stochastic sampling.

From the paper abstract and ACS article page:

- The method models DEL sequencing as sampling from a binomial distribution.
- It compares an observed population fraction against an expected population fraction.
- The normalized form is intended to be less sensitive to sequencing depth and library composition than naive enrichment metrics.
- The metric is intended to support quantitative comparison of enrichment across parallel DEL selections.

The paper's stated target use case is comparison of DEL features, especially `n`-synthons, across selections such as target versus control. This is importantly different from replicate-aware differential-expression style methods like `edgeR`.

## High-Level Formula

The ACS article page exposes the interpretation of the score, but the equation images themselves are not text-searchable in the page source. For this reason, DELT should document the standard binomial derivation explicitly and label it as an interpretation consistent with the paper's prose and later DEL literature.

At a high level:

- Let `p_o` be the observed fraction for a feature in a sequenced sample.
- Let `p_i` be the expected fraction for that feature in the input population.
- Let `q_i = 1 - p_i`.
- Let `n` be the total sequencing depth for that sample.

The standard binomial z-score is:

```text
z = (p_o - p_i) / sqrt(p_i * q_i / n)
```

The normalized z-score removes the explicit sampling-depth term:

```text
z_n = (p_o - p_i) / sqrt(p_i * q_i)
```

This is consistent with the ACS page summary that the final normalized z-score "showed a low sensitivity to sampling" and that, for a fixed fold-enrichment `p_o / p_i`, the score scales by a factor involving `sqrt(p_i / q_i)`.

## Interpretation and Scope in the Paper

The paper is not presenting a replicate-based significance test. It is presenting a per-sample enrichment metric that can still be useful when biological replicates are unavailable.

The intended scope is broader than just full compounds:

- monosynthons
- disynthons
- trisynthons
- more generally, `n`-synthons

This matters because the paper's comparative plots and thresholding guidance are about enriched DEL features, not specifically about compound-level rows in a demultiplex table.

The ACS article page also states the authors' practical rule of thumb:

- they view `z_n >= 1` for an `n`-synthon as an indicator of significant enrichment
- they note that different DEL compositions and selection protocols may require an adjusted threshold

That threshold guidance is useful as context, but DELT should treat it as ranking guidance rather than a hard universal decision rule.

## How This Compares to DELT Today

DELT currently exposes two enrichment modes in [`src/delt_hit/cli/analyse/api.py`](/Users/adrianomartinelli/projects/delt-hit/src/delt_hit/cli/analyse/api.py):

- `counts`: a simple count-based summary that averages counts across replicates within a group and computes `protein - no_protein`
- `edgeR`: a replicate-aware generalized linear model using TMM normalization and FDR-adjusted tests

Those two modes leave a gap:

- `edgeR` is the statistically stronger option when replicates exist
- `counts` is lightweight but does not model stochastic sampling noise
- users with no replicates currently do not have an intermediate method that is more principled than simple count differences

The normalized z-score is a good fit for that gap because it does not depend on biological replicates and still uses a probabilistic model of observed enrichment.

## Minimal DELT Changes

The smallest DELT change should add a new analysis mode without expanding the data model or analysis scope beyond what DELT already supports.

### Proposed CLI Addition

Extend:

```bash
delt-hit analyse enrichment --config_path <path/to/config.yaml> --name <experiment-name> --method <method>
```

to accept:

```text
--method=zscore
```

This would live beside the current `counts` and `edgeR` branches in [`src/delt_hit/cli/analyse/api.py`](/Users/adrianomartinelli/projects/delt-hit/src/delt_hit/cli/analyse/api.py).

### Proposed v1 Scope

Keep the first implementation deliberately narrow:

- score full compound rows only
- reuse the current `prepare_data()` output shape
- keep the current analysis config shape
- do not add monosynthon/disynthon/trisynthon aggregation in the same change

This is narrower than the paper's full scope, but it matches DELT's current enrichment pipeline and keeps the implementation small.

### Expected Probability for v1

For a minimal implementation, DELT should derive the expected probability from existing project information rather than introducing a new input file.

Preferred default:

- assume a uniform expected probability over the full enumerated library
- compute `p_i = 1 / library_size`
- derive `library_size` from existing DELT library/config metadata rather than from only the observed subset of sequenced rows

Observed probability would be computed per selection as:

```text
p_o = count / total_reads_in_selection
```

This is closer to the intended use in the paper than using only the subset of compounds that happened to be observed in a given count table.

### Proposed Output Shape

To stay parallel with existing analysis modes, write results under:

```text
<save_dir>/<experiment_name>/zscore/
```

Minimum useful outputs:

- per-selection stats table with:
  - `code_*`
  - raw `count`
  - `observed_fraction`
  - `expected_fraction`
  - `zscore`
- `hits.csv` ranked by descending z-score
- a merged summary table for side-by-side comparison of `protein` and `no_protein` selections

Optionally, the generator could also emit one CSV per selection, mirroring the current export style used by `counts` and `edgeR`.

### Suggested Implementation Shape

The smallest code change would look like this:

- add a `case 'zscore': ...` branch in `Analyse.enrichment(...)`
- add a `zscore_rscript(...)` generator adjacent to `counts_rscript(...)` and `edgeR_rscript(...)`
- keep using the existing `prepare_data(...)` output files (`data.csv`, `samples.csv`)
- write outputs into a new `zscore/` directory

This keeps the change localized to the current analysis entrypoint instead of introducing a new analysis subsystem.

## Explicit Non-Goals for v1

To keep the first implementation minimal, do not include the following in the same change:

- `n`-synthon aggregation
- confidence interval exports
- new dashboard support
- automatic threshold-based hit calling beyond rank-ordered exports
- replacement of `edgeR` for replicate-rich experiments

## Documentation Changes Once Implemented

When the feature exists, DELT should update:

- [`README.md`](/Users/adrianomartinelli/projects/delt-hit/README.md)
- [`documentation/cli.md`](/Users/adrianomartinelli/projects/delt-hit/documentation/cli.md)

The wording should describe `zscore` as the recommended analysis mode for selections without replicates, while still positioning `edgeR` as the stronger choice when replicate structure is available.

## Test Plan for a Future Implementation

### Unit-Level Checks

- higher observed counts at fixed library size and sequencing depth should increase z-score
- counts near the expected abundance should give z-scores near zero
- the normalized score should preserve comparable ranking when the same enrichment pattern is observed at different sequencing depths

### Integration Checks

- `Analyse.enrichment(..., method="zscore")` should create the expected `zscore/` output directory
- expected files should be written successfully from the generated script
- existing `counts` and `edgeR` behavior should remain unchanged

### Acceptance Scenario

A practical acceptance test should confirm that:

- one `protein` selection and one `no_protein` selection can be analyzed without replicates
- compounds are ranked by z-score
- the output is useful for identifying candidate enriched compounds even when replicate-based statistics are unavailable

## Assumptions

- This note is documentation-first only and does not implement the z-score analysis itself.
- The first DELT implementation should stay at compound-level analysis because that matches the current analysis pipeline.
- The paper-backed motivation is strongest for `n`-synthon enrichment, but compound-level scoring is the lowest-risk DELT entry point.
- For v1, expected population should come from full library size rather than from the observed subset in a count file.

## Sources

- Faver et al. abstract and metadata via PubMed: [https://pubmed.ncbi.nlm.nih.gov/30672692/](https://pubmed.ncbi.nlm.nih.gov/30672692/)
- ACS article landing page and abstract: [https://pubs.acs.org/doi/10.1021/acscombsci.8b00116](https://pubs.acs.org/doi/10.1021/acscombsci.8b00116)
- ACS abstract page summary with the normalized z-score discussion: [https://pubs.acs.org/doi/abs/10.1021/acscombsci.8b00116](https://pubs.acs.org/doi/abs/10.1021/acscombsci.8b00116)
- Later DEL application citing the same normalized z-score framing: [https://pmc.ncbi.nlm.nih.gov/articles/PMC10357458/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357458/)
