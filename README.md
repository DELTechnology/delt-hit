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


## Initialization

Initialize folder structure:
```bash
delt-cli init
```


## Simulation

Create configuration file for simulation:
```bash
delt-cli simulate init
```

Generate reads (with or without erros):
```bash
delt-cli simulate run experiments/default-*/config.yml
```


## Demultiplexing

Initialize folder structure for demultiplexing:
```bash
delt-cli demultiplex init
```

Run demultiplexing:
```bash
delt-cli demultiplex run experiments/default-*/config.yml
```

Create codon lists for the library specified in the configuration file:
```bash
delt-cli demultiplex create-lists experiments/default-*/config.yml
```

Create input files for Cutadapt:
```bash
delt-cli demultiplex create-cutadapt-input experiments/default-*/config.yml
```

Run bash script:
```bash
bash experiments/default-*/cutadapt_input_files/demultiplex.sh
```

Compute count tables for the final reads:
```bash
delt-cli demultiplex compute-counts experiments/default-*/config.yml experiments/default-*/cutadapt_output_files/reads_with_adapters.gz output_dir
```

Convert old structure file to new configuration file:
```bash
delt-cli demultiplex convert structure.txt
```


## Quality Control

Plot codon hits:
```bash
delt-cli qc plot experiments/default-*
```

Print report:
```bash
delt-cli qc report experiments/default-*
```


## Workflow

Initialize the folder structure and move the library, selection, and FASTQ files to the corresponding directories:
```bash
delt-cli init
mv /path/to/input.fastq.gz fastq_files
mv /path/to/library.xlsx libraries
mv /path/to/selection.xlsx selections
```

Switch to the ```smiles``` branch and compute the SMILES and some chemical properties of a library:
```bash
git checkout smiles
delt-cli compute smiles libraries/NF.xlsx
delt-cli compute properties libraries/smiles/NF_smiles.txt.gz
delt-cli compute plot libraries/properties/properties_L1.txt.gz
```

Switch back to the ```main``` branch, create the configuration file, demultiplex the FASTQ file, and store the counts according to the selections defined in the selection file:
```bash
git checkout main
delt-cli demultiplex init -f fastq_files/input.fastq.gz -l libraries/NF.xlsx -s selections/selection.xlsx
delt-cli demultiplex run experiments/default-*/config.yml
```

Switch to the ```quality_control``` branch, report and plot the results:
```bash
git checkout quality_control
delt-cli qc report experiments/default-*
delt-cli qc plot experiments/default-*
```

If there are no library, selection, or FASTQ files available, one can generate them using the following commands:
```bash
git checkout main
delt-cli simulate init
delt-cli simulate run config.yml
```

The default initialization of the simulation generates a new library and a selection template that contain random codons. Alternatively, one can pass existing files using the following commands:
```bash
delt-cli simulate init -l libraries/library.xlsx -s selections/selection.xlsx -f fastq_files/input.fastq.gz -o fastq_files/input.fastq.gz
delt-cli simulate run config.yml
```
