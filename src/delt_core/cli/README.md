# API

## Design (This code base might come from Alice)
delt-cli design init  # create default config file 
delt-cli design init  --min-intra-codon-dist 3 --avoid-prefix-postfix-overlaps # create default config file
delt-cli design run <PATH_TO_CONFIG_FILE>  # run design

## Compute
delt-cli compute smiles  # PR#3
delt-cli compute embeddings
delt-cli compute diversity

## Demultiplex
delt-cli demultiplex init  # create default config file
delt-cli demultiplex init --from-fastq <PATH_TO_FASTQ>
delt-cli demultiplex init --from-old-stucture <PATH_TO_OLD_STRUCT_FILE>
delt-cli demultiplex init --from-library <PATH_TO_EXCEL_FILE_OF_LIBRARY>
delt-cli demultiplex run <PATH_TO_CONFIG_FILE>  # run demultiplexing

## Simulation
delt-cli simulate init  # create default config file
delt-cli simulate run <PATH_TO_CONFIG_FILE>  # run simulation

## Quality Control
delt-cli qc 

## DB Management
delt-cli db register library 
delt-cli db register fastq