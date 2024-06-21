# Package `delt-core`
Core functionalities to work with DECL libraries

## Installation

Install the package for development purposes.
Navigate to the root folder of the package and run the following command:

```bash
pip install -e ".[dev,test]"
```

Create a `.env` file to store environment configurations, keys and secrets.
```bash
touch .env
```
These configurations can be accessed using the `python-dotenv` package.


## SMILES construction

Standard use:
```bash
delt-cli compute smiles library1.xlsx
```

Hybridization of two libraries (the order of the libraries must match the final sequence in the 5'-to-3' direction, see README):
```bash
delt-cli compute smiles library1.xlsx library2.xlsx
```

## Simulation
```bash
delt-cli simulate
```

## Demultiplexing
```bash
delt-cli demultiplex run config.yml
```

```bash
delt-cli demultiplex init config.yml
```

```bash
delt-cli demultiplex create-lists config.yml
```

```bash
delt-cli demultiplex create-cutadapt-input config.yml
```

```bash
delt-cli demultiplex compute-counts reads_with_adapters.gz
```

```bash
delt-cli demultiplex convert structure.txt
```

## Quality Control
```bash
delt-cli qc
```
