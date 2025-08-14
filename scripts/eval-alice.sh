#!/bin/bash

cd /work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB || exit
export ROOT="/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB"

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-112-123 \
-f $ROOT/fastq_files/394251_1-2508_AL_linker_1_AGP_S1_R1_001.fastq.gz \
-l $ROOT/libraries/AL_AGP.xlsx \
-s $ROOT/selections/selections_250812.xlsx

EXPERIMENT_DIR=/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/experiments/alice-112-123-2025-08-14-09-17-03
delt-cli demultiplex run $EXPERIMENT_DIR/config.yml
delt-cli qc report $EXPERIMENT_DIR
delt-cli qc plot $EXPERIMENT_DIR

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-124-135 \
-f $ROOT/fastq_files/394251_2-2508_AL_linker_2_fc_S2_R1_001.fastq.gz \
-l $ROOT/libraries/AL_AGP_FC.xlsx \
-s $ROOT/selections/selections_250812.xlsx

EXPERIMENT_DIR=/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/experiments/alice-124-135-2025-08-14-09-34-40
delt-cli demultiplex run $EXPERIMENT_DIR/config.yml
delt-cli qc report $EXPERIMENT_DIR
delt-cli qc plot $EXPERIMENT_DIR

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-136-151 \
-f $ROOT/fastq_files/394251_3-2508_AL_linker_3_CA_S3_R1_001.fastq.gz \
-l $ROOT/libraries/AL_CA.xlsx \
-s $ROOT/selections/selections_250812.xlsx

EXPERIMENT_DIR=/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/experiments/alice-136-151-2025-08-14-09-35-45
delt-cli demultiplex run $EXPERIMENT_DIR/config.yml
delt-cli qc report $EXPERIMENT_DIR
delt-cli qc plot $EXPERIMENT_DIR