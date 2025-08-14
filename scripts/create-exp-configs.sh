#!/bin/bash

cd /work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB || exit
export ROOT="/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB"

# --- #

delt-cli demultiplex init \
-e eval-1-27 \
-f $ROOT/fastq_files/368061_1-241105_AG_BZ_NC_pool1_NF_S3_R1_001.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-28-54 \
-f $ROOT/fastq_files/368061_2-241105_AG_BZ_NC_pool2_GB_S1_R1_001.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-55-81 \
-f $ROOT/fastq_files/368061_3-241105_AG_BZ_NC_pool3_2p1_S2_R1_001.fastq.gz \
-l $ROOT/libraries/NF-KS_corrected_codes3RC.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-82-87 \
-f $ROOT/fastq_files/20211210.A-o26804_1_3-NF_GB_YO_screenAG1_R1.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-88-96 \
-f $ROOT/fastq_files/o306001_2-sample2_701_501_TRIM_BAF1_S58_R1_001.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-97-102 \
-f $ROOT/fastq_files/20211210.A-o26804_1_3-NF_GB_YO_screenAG1_R1.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- #

delt-cli demultiplex init \
-e eval-103-111 \
-f $ROOT/fastq_files/o306001_4-sample4_703_503_TRIM_BAF2_S57_R1_001.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-112-123 \
-f $ROOT/fastq_files/394251_1-2508_AL_linker_1_AGP_S1_R1_001.fastq.gz \
-l $ROOT/libraries/AL_AGP.xlsx \
-s $ROOT/selections/selections_250812.xlsx

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-124-135 \
-f $ROOT/fastq_files/394251_2-2508_AL_linker_2_fc_S2_R1_001.fastq.gz \
-l $ROOT/libraries/AL_AGP_FC.xlsx \
-s $ROOT/selections/selections_250812.xlsx

# --- ALICE --- #

delt-cli demultiplex init \
-e alice-136-151 \
-f $ROOT/fastq_files/394251_3-2508_AL_linker_3_CA_S3_R1_001.fastq.gz \
-l $ROOT/libraries/AL_CA.xlsx \
-s $ROOT/selections/selections_250812.xlsx