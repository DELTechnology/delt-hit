# Example Single-Display Workflow

Full end-to-end DELT-Hit workflow example using a synthetic single-display DEL.

## Files

- `example-single-display.xlsx` — library definition
- `campaign.fastq.gz` — example sequencing data (download if not present)
- `analysis.yaml` — enrichment analysis configuration

## Reproduce

Download the FASTQ file if not present:

```bash
bash download.sh
```

Run the full workflow:

```bash
bash run.sh
```

Outputs are written under `campaign/`.
