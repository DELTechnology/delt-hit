# Synthetic Enrichment Benchmark

This benchmark exists to compare ranking recovery of known positive controls under
configurable count noise. It starts from per-selection count tables on purpose:
the goal is to isolate enrichment-analysis behavior rather than test sequencing,
adapter trimming, or barcode decoding.

The generator writes:

- six DELT-Hit-compatible count tables across `no_protein_rep1..3` and `protein_rep1..3`
- `config.yaml` ready for `delt-hit analyse enrichment`
- `truth.tsv` with seeded positive tiers and target means
- `manifest.json` with the exact generation parameters

Noise in this benchmark means abundance heterogeneity plus replicate-to-replicate
count variation. It does not model sequencing corruption.

## Generate A Synthetic Enrichment Dataset

```bash
./.venv/bin/python benchmarks/enrichment/generate_synthetic_enrichment.py \
  --building-blocks-per-cycle 100 \
  --positive-low 50 \
  --positive-medium 25 \
  --positive-high 10 \
  --fc-low 1.5 \
  --fc-medium 3.0 \
  --fc-high 8.0 \
  --dispersion 0.20 \
  --library-size-cv 0.10 \
  --seed 13 \
  --output-dir target/synthetic_enrichment \
  --experiment-name synthetic_enrichment_default
```

The default model uses:

- a heterogeneous per-compound background mean
- explicit low / medium / high positive spike-ins
- negative-binomial replicate sampling
- three `protein` replicates and three `no_protein` replicates

`truth.tsv` includes metadata comment lines describing the distribution family and
the exact parameter values used for that run.

## Run The Comparison

```bash
module load r-light/4.4.1
./.venv/bin/python benchmarks/enrichment/compare_synthetic_enrichment.py \
  --dataset-dir target/synthetic_enrichment/synthetic_enrichment_default \
  --top-k 100
```

This runs both existing enrichment paths:

- `counts`
- `edgeR`

and writes a comparison bundle under `<dataset>/comparison/`:

- `comparison.tsv`: one row per compound with truth labels, scores, and ranks
- `comparison.json`: machine-readable metrics summary
- `top_hits.tsv`: side-by-side top-ranked compounds from each method
- `summary.md`: compact human-readable report

## Interpreting The Report

The report is relative rather than threshold-based. The main questions are:

- How many seeded positives appear in the top `K` results?
- Does `edgeR` recover more `positive_low` compounds than raw counts?
- Are negatives staying out of the top-ranked region?
- How do the median and mean positive ranks compare across methods?

## Analysis

Minimal run:

```bash
./.venv/bin/python benchmarks/enrichment/generate_synthetic_enrichment.py \
  --output-dir target/synthetic_enrichment \
  --experiment-name synthetic_enrichment_default

./.venv/bin/python benchmarks/enrichment/compare_synthetic_enrichment.py \
  --dataset-dir target/synthetic_enrichment/synthetic_enrichment_default \
  --top-k 100
```

Main output files:

- `target/synthetic_enrichment/synthetic_enrichment_default/comparison/summary.md`
- `target/synthetic_enrichment/synthetic_enrichment_default/comparison/comparison.tsv`
