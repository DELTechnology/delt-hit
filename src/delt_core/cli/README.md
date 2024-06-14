# API

## Design (This code base might come from Alice)
delt-cli design init  # create default config file 
delt-cli design init  --min-intra-codon-dist 3 --avoid-prefix-postfix-overlaps # create default config file
delt-cli design run <PATH_TO_CONFIG_FILE>  # run design

## Demultiplex
delt-cli demultiplex convert <PATH_TO_OLD_STRUCT_FILE> <PATH_TO_NEW_STRUCT_FILE> -> codon lists -> cutadapt_input -> demultiplex.sh
delt-cli demultiplex init --from-library <PATH_TO_EXCEL_FILE_OF_LIBRARY> -> codon lists -> cutadapt_input -> demultiplex.sh

delt-cli demultiplex create-lists <PATH_TO_EXCEL_FILE_OF_LIBRARY>
delt-cli demultiplex create-cutadapt-input <PATH_TO_STRUCTUR_FILE> --input input.fastq.gz 
./cutadapt-input-files/demultiplex.sh 
delt-cli demultiplex compute-counts <PAHT_TO_FILE_WITH_ADAPTER> 
delt-cli demultiplex report

## Compute
delt-cli compute smiles  # PR#3
delt-cli compute embeddings
delt-cli compute features # hydrophilic, pH, ..., etc
delt-cli compute diversity

## Simulation
delt-cli simulate init  # create default config file
delt-cli simulate run <PATH_TO_CONFIG_FILE>  # run simulation

## Quality Control
delt-cli qc CMD

## DB Management
delt-cli db register library 
delt-cli db register fastq