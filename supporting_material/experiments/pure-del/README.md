# Pure-DEL Re-Analysis

Re-analysis of the Keller et al. dataset using DELT-Hit.  
Publication DOI: [10.1126/science.adn3412](https://doi.org/10.1126/science.adn3412)

## Files

- `lane-1.xlsx`, `lane-2.xlsx` — library definitions
- `analysis.yaml` — enrichment analysis configuration
- `compare_selections.py` — comparison script against published counts

## Reproduce

Download the raw FASTQ files and published selection tables:

```bash
bash download.sh
```

Run DELT-Hit for both lanes:

```bash
bash run.sh lane-1
bash run.sh lane-2
```

This produces per-selection counts under `lane-1/selections/` and `lane-2/selections/`.

## Compare against published counts

```bash
pixi run python compare_selections.py lane-1
pixi run python compare_selections.py lane-2
```

Comparison tables are written to `comparison/lane-1/` and `comparison/lane-2/`. Each selection produces a detailed CSV with an `identical` column plus a `report.csv` summary of observed-compound counts, total counts, and mismatch totals.
