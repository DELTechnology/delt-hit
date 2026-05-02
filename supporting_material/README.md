# Supporting Material

This directory collects the reproducible example workflows and dataset re-analyses that accompany DELT-Hit.

- `supporting_material/experiments/example-single-stranded/` contains a full single-stranded DEL workflow example, including the Excel template, FASTQ download helper, analysis configuration, and `run.sh` workflow script.
- `supporting_material/experiments/example-dual-display/` contains a standalone dual-display example focused on library enumeration and visualization only.
- `supporting_material/experiments/favalli/` and `supporting_material/experiments/pure-del/` contain rerunnable re-analyses of previously published datasets, together with comparison scripts for validating DELT-Hit output against the published counts.
- `supporting_material/experiments/johanna/` contains an additional experiment-specific workflow example.

For Polybox-based datasets, download troubleshooting notes and validated direct-download URLs are collected in `supporting_material/polybox-download.md`.

## Full Example DELT-Hit Workflow

The single-stranded example in `supporting_material/experiments/example-single-stranded/` demonstrates the full DELT-Hit workflow from initialization through demultiplexing, counting, visualization, enrichment analysis, and library featurization.

The example uses:

- `example-single-stranded.xlsx`
- `campaign.fastq.gz`
- `analysis.yaml`

If `campaign.fastq.gz` is not present locally, download it from figshare first:

```bash
cd supporting_material/experiments/example-single-stranded
bash download.sh
```

Then run the workflow commands listed in `run.sh`:

```bash
cd supporting_material/experiments/example-single-stranded
bash run.sh
```

`run.sh` executes the full workflow in order:

- `delt-hit init`
- demultiplex preparation, report generation, QC, and processing
- enumeration and top-hit library extraction
- library visualization, properties, and dashboard generation
- enrichment setup plus the counts- and edgeR-based R analyses
- full library enumeration and molecular representations

Outputs are written under `supporting_material/experiments/example-single-stranded/campaign/`.

## Dual-Display Example

The dual-display example in `supporting_material/experiments/example-dual-display/` is a lightweight standalone example for enumeration only. It does not require a FASTQ file because the documented workflow only runs the library visualization and enumeration steps.

Run the commands in `run.sh` with:

```bash
cd supporting_material/experiments/example-dual-display
bash run.sh
```

This script runs:

- `delt-hit init --excel_path example-dual-display.xlsx`
- `delt-hit visualize enumerate --config_path=campaign/config.yaml`
- `delt-hit library enumerate --config_path=campaign/config.yaml`

Outputs are written under `supporting_material/experiments/example-dual-display/campaign/`.

## Re-Analysis of Previously Published Datasets

The publication re-analysis folders reproduce the DELT-Hit processing for the Favalli et al. and Pure-DEL datasets and provide comparison scripts that align DELT-Hit counts with the published selection counts.

### Favalli et al.

Source folder: `supporting_material/experiments/favalli/`  
Publication DOI: [10.1038/s41557-021-00660-y](https://doi.org/10.1038/s41557-021-00660-y)

Download the raw FASTQ files and published comparison table first:

```bash
cd supporting_material/experiments/favalli
bash download.sh
```

This script:

- downloads `20190812.A-1907_NF2GB2_s1_R1.fastq.gz`
- downloads `20190812.A-1907_NF2GB2_s2_R1.fastq.gz`
- downloads `20190812.A-1907_NF2GB2_s1_R1.fasta.gz`
- downloads `20190812.A-1907_NF2GB2_s2_R1.fasta.gz`
- downloads the published evaluation table into `supporting_material/experiments/favalli/published/`

Run DELT-Hit for the `lane-1-fasta` workbook prefix:

```bash
cd supporting_material/experiments/favalli
bash run.sh lane-1-fasta
```

`run.sh` performs workbook initialization, demultiplex preparation, QC, and DELT-Hit count generation for the requested workbook prefix.

To generate all files required by both the enrichment workflow and the comparison script for `lane-1-fasta`, run:

```bash
cd supporting_material/experiments/favalli
bash run.sh lane-1-fasta
```

This produces the per-selection DELT-Hit counts under `lane-1-fasta/selections/`, which are the inputs required by both the checked-in `analysis.yaml` enrichment workflow and the comparison script.

Run the Favalli enrichment analysis explicitly afterward:

```bash
cd supporting_material/experiments/favalli
for exp in ca9_ss ca9_ds;
do
    pixi run delt-hit analyse enrichment \
      --config_path=analysis.yaml \
      --name=${exp} \
      --method=counts
    Rscript --vanilla lane-1-fasta/analysis/${exp}/counts/enrichment_counts.R
    
    pixi run delt-hit analyse enrichment \
      --config_path=analysis.yaml \
      --name=${exp} \
      --method=edgeR
    Rscript --vanilla lane-1-fasta/analysis/${exp}/edgeR/enrichment_edgeR.R
done
```

This workflow produces:

- per-selection DELT-Hit counts under `lane-1-fasta/selections/`
- enrichment analysis inputs and generated R scripts under `lane-1-fasta/analysis/`
- counts- and edgeR-based enrichment outputs for the analyses defined in `analysis.yaml`

Other workbook prefixes remain available for count comparison only:

```bash
cd supporting_material/experiments/favalli
bash run.sh lane-1
bash run.sh lane-2
bash run.sh lane-2-fasta
```

Compare DELT-Hit counts against the published counts:

```bash
cd supporting_material/experiments/favalli
pixi run python compare_selections.py lane-1-fasta
```

Comparison tables are written to `supporting_material/experiments/favalli/comparison/<prefix>/` for each of the four prefixes above. Each output folder contains one detailed CSV per selection plus a `report.csv` summary with the observed-compound counts, total counts, and mismatch totals for the published and DELT-Hit counts. Inspect the detailed CSV `identical` column and the summary report to confirm whether the published and DELT-Hit counts match for each selection.

### Pure-DEL

Source folder: `supporting_material/experiments/pure-del/`  
Publication DOI: [10.1126/science.adn3412](https://doi.org/10.1126/science.adn3412)

Download the raw FASTQ files and published comparison tables first:

```bash
cd supporting_material/experiments/pure-del
bash download.sh
```

This script:

- downloads `328321_1-230906_MK_TSD501701_S1_R1_001_SP_fl1.fastq.gz`
- downloads `328321_2-230906_MK_TSD504704_S2_R1_001_UP_fl2.fastq.gz`
- downloads the Zenodo ZIP with published selection counts
- extracts the `selection_*_.txt` files into `supporting_material/experiments/pure-del/published/`

Run DELT-Hit for a given lane: `run.sh` performs lane initialization, demultiplex preparation, QC, and count generation for the requested lane.

To generate all per-selection counts required by both the enrichment workflow and the comparison script, run both lanes:

```bash
cd supporting_material/experiments/pure-del
bash run.sh lane-1
bash run.sh lane-2
```

The cross-lane enrichment analyses in `analysis.yaml` require those two selection directories, so run the enrichment steps explicitly afterward:

```bash
cd supporting_material/experiments/pure-del
pixi run delt-hit analyse enrichment --config_path=analysis.yaml --name=his_pure_up --method=counts
Rscript --vanilla lane-2/analysis/his_pure_up/counts/enrichment_counts.R

pixi run delt-hit analyse enrichment --config_path=analysis.yaml --name=his_pure_up --method=edgeR
Rscript --vanilla lane-2/analysis/his_pure_up/edgeR/enrichment_edgeR.R

pixi run python enrichment.py \
  --data-dir lane-2/analysis/his_pure_up \
  --output-dir enrichment/his_pure_up
```

Repeat the same `delt-hit analyse enrichment`, `Rscript`, and `pixi run python enrichment.py` pattern for `his_pure_sp`, `dyna_up`, and `dyna_sp`.

This workflow produces:

- per-selection DELT-Hit counts under `lane-1/selections/` and `lane-2/selections/`
- enrichment analysis inputs and R scripts under `lane-1/analysis/` and `lane-2/analysis/`
- counts- and edgeR-based enrichment outputs for `his_pure_up`, `his_pure_sp`, `dyna_up`, and `dyna_sp`
- summary plots and CSV exports under `supporting_material/experiments/pure-del/enrichment/`

Compare DELT-Hit counts against the published counts:

```bash
cd supporting_material/experiments/pure-del
pixi run python compare_selections.py lane-1
# or
pixi run python compare_selections.py lane-2
```

Comparison tables are written to `supporting_material/experiments/pure-del/comparison/lane-1/` or `supporting_material/experiments/pure-del/comparison/lane-2/`. Each output folder contains one detailed CSV per selection plus a `report.csv` summary with the observed-compound counts, total counts, and mismatch totals for the published and DELT-Hit counts. Inspect the detailed CSV `identical` column and the summary report to verify the agreement between the published selections and DELT-Hit output.

## Additional Example

### Johanna

Source folder: `supporting_material/experiments/johanna/`

Download the FASTQ file first:

```bash
cd supporting_material/experiments/johanna
bash download.sh
```

This downloads `407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz`, which is the filename expected by `JP-1/config.yaml`.

Then run:

```bash
cd supporting_material/experiments/johanna
bash run.sh
```
