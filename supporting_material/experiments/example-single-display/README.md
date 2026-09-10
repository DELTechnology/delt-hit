# Example Single-Display Workflow

Full end-to-end DELT-Hit workflow example using a synthetic single-display DEL.

## Files

- `example-single-display.xlsx` — library definition
- `campaign.fastq.gz` — example sequencing data (approximately 6.6 GB compressed; downloaded separately)
- `analysis.yaml` — enrichment analysis configuration

## Reproduce

Follow the [Conda installation guide](../../../README.md#2-installation-guide), including R packages. From this directory, activate the environment and run:

```bash
conda activate delt-hit
bash download.sh
bash run.sh
```

The download script skips an existing non-empty FASTQ file. `run.sh` lists and times all workflow commands, including execution of the generated R scripts. The full-library commands at the end are optional for selection analysis.

## Expected output and runtime

Outputs are written under `campaign/`:

- `config.yaml`: parsed experiment and library configuration.
- `demultiplex/`: Cutadapt input files, execution script, logs, and mapped reads.
- `qc/`: read-retention report and barcode-recovery plots.
- `selections/`: selection-level count tables.
- `library/`: enumerated compounds, chemistry visualizations, and molecular properties.
- `analysis/`: counts, edgeR, and normalized z-score results, including `stats.csv` and `hits.csv` after R execution.
- `representations/`: molecular representations from optional full-library processing.

Allow one to several hours on a desktop with 32 GB RAM and 16 CPU cores, excluding the dataset download. This follows the protocol's indicative timing guidance; it is not a separately measured demo benchmark. Runtime depends on hardware and which optional stages are run.

The current `bert` CLI option writes Morgan fingerprints to `bert.npz`, not transformer embeddings. See the [main README](../../../README.md#3-demo) for this implementation limitation and the full output tree.
