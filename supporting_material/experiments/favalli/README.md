# Favalli et al. Re-Analysis

Re-analysis of the Favalli et al. dataset using DELT-Hit.  
Publication DOI: [10.1038/s41557-021-00660-y](https://doi.org/10.1038/s41557-021-00660-y)

## Files

- `lane-1.xlsx` — library definition for lane 1
- `analysis.yaml` — enrichment analysis configuration
- `compare_selections.py` — comparison script against published counts

## Reproduce

Download the raw FASTA file and published counts:

```bash
bash download.sh
```

Run DELT-Hit for lane 1:

```bash
bash run.sh
```

This produces per-selection counts under `lane-1/selections/`.

## Compare against published counts

```bash
pixi run python compare_selections.py
```

Comparison tables are written to `comparison/lane-1/`. Each selection produces a detailed CSV with an `identical` column plus a `report.csv` summary of observed-compound counts, total counts, and mismatch totals.
