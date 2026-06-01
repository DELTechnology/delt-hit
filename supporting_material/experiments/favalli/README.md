# Favalli et al. Re-Analysis

Re-analysis of the Favalli et al. dataset using DELT-Hit.  
Publication DOI: [10.1038/s41557-021-00660-y](https://doi.org/10.1038/s41557-021-00660-y)

## Files

- `lane-1.xlsx` — library definition for lane 1
- `analysis.yaml` — enrichment analysis configuration
- `compare_selections.py` — comparison script against published counts
- `enrichment_scores.py` — cross-method enrichment score comparison and hit export

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

## Enrichment analysis

Run the enrichment setup and execute the generated R scripts:

```bash
pixi run delt-hit analyse enrichment \
  --analysis_config=analysis.yaml \
  --save_dir=lane-1/analysis \
  --name=ca9_ds \
  --method=counts
pixi run Rscript --vanilla lane-1/analysis/counts/ca9_ds/enrichment_counts.R

pixi run delt-hit analyse enrichment \
  --analysis_config=analysis.yaml \
  --save_dir=lane-1/analysis \
  --name=ca9_ds \
  --method=edgeR
pixi run Rscript --vanilla lane-1/analysis/edgeR/ca9_ds/enrichment_edgeR.R

for selection in 4_2 5_2 6_2; do
  pixi run delt-hit analyse enrichment \
    --config_path=lane-1/config.yaml \
    --save_dir=lane-1/analysis \
    --counts=lane-1/selections/${selection}/counts.txt \
    --name=${selection} \
    --method=z_score
  pixi run Rscript --vanilla lane-1/analysis/z_score/${selection}/enrichment_z_score.R
done
```

Then compare enrichment scores across methods and export hit lists:

```bash
pixi run python enrichment_scores.py
```

Defaults to `--analysis-root lane-1/analysis`, `--analysis-config analysis.yaml`, and `--output-dir enrichment/${exp}`. The script produces:

- `hits_<N>.csv` — top-N hit lists for counts, edgeR, and z-score
- `top_<N>_<code>_counts.pdf/png` — bar plots of building-block frequency among top hits
- `controls.csv` — rank and score of query compounds across all three methods

## Compare against published counts

```bash
pixi run python compare_selections.py
```

Comparison tables are written to `comparison/lane-1/`. Each selection produces a detailed CSV with an `identical` column plus a `report.csv` summary of observed-compound counts, total counts, and mismatch totals.
