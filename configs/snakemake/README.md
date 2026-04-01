# Snakemake Execution Profiles

Each subdirectory is a Snakemake profile — a `config.yaml` that sets execution
defaults so you don't have to repeat flags on every invocation.

## Local

Runs all jobs on the current machine using up to 4 cores.

```bash
snakemake --profile configs/snakemake/local/
snakemake --profile configs/snakemake/local/ -n          # dry-run
snakemake --profile configs/snakemake/local/ --forceall  # re-run everything
```

Change `cores:` in [local/config.yaml](local/config.yaml) to match your machine.

## SLURM

Submits each rule as a separate `sbatch` job. Requires the SLURM executor plugin:

```bash
uv add snakemake-executor-plugin-slurm   # install once
snakemake --profile configs/snakemake/slurm/
```

**Before first use**, edit [slurm/config.yaml](slurm/config.yaml) and set
`slurm_partition` under `default-resources:` to your cluster's partition name
(e.g. `short`, `compute`, `gpu`).

Per-rule memory and walltime overrides are already configured for the
resource-heavy rules (`demultiplex`, `enumerate`, `represent`). Add more
entries under `set-resources:` as needed.

## Dry-run before submitting to SLURM

Always check the execution plan before submitting to a cluster:

```bash
snakemake --profile configs/snakemake/slurm/ -n
snakemake --profile configs/snakemake/slurm/ --dag | dot -Tpng > dag.png
```
