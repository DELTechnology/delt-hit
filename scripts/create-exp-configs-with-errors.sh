#!/bin/bash

cd /work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB || exit
ROOT="/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB"

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-1-27-e=1.9 \
-f $ROOT/fastq_files/368061_1-241105_AG_BZ_NC_pool1_NF_S3_R1_001.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-28-54-e=1.9 \
-f $ROOT/fastq_files/368061_2-241105_AG_BZ_NC_pool2_GB_S1_R1_001.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-55-81-e=1.9 \
-f $ROOT/fastq_files/368061_3-241105_AG_BZ_NC_pool3_2p1_S2_R1_001.fastq.gz \
-l $ROOT/libraries/NF-KS_corrected_codes3RC.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9 --errors C4 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-82-87-e=1.9 \
-f $ROOT/fastq_files/20211210.A-o26804_1_3-NF_GB_YO_screenAG1_R1.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-88-96-e=1.9 \
-f $ROOT/fastq_files/o306001_2-sample2_701_501_TRIM_BAF1_S58_R1_001.fastq.gz \
-l $ROOT/libraries/NF.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-97-102-e=1.9 \
-f $ROOT/fastq_files/20211210.A-o26804_1_3-NF_GB_YO_screenAG1_R1.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #

delt-cli demultiplex init \
-r "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB" \
-e eval-103-111-e=1.9 \
-f $ROOT/fastq_files/o306001_4-sample4_703_503_TRIM_BAF2_S57_R1_001.fastq.gz \
-l $ROOT/libraries/GB2.xlsx \
-s $ROOT/selections/selections_250603.xlsx \
--errors C1 1.9 --errors C2 1.9 --errors C3 1.9

# --- #