#!/bin/sh
UCSD_ROOT_PATH="/fs/home/jw865/peptide-catalysis"
DATA_FILE="${UCSD_ROOT_PATH}/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv"
UCSD_outname="type2_hit"
Gamma0=500
Gamma1=0.05
Prior=0.2
rm -R ${UCSD_ROOT_PATH}/temp/cv/*

N=$(($(wc -w < ${DATA_FILE}) - 1))
for i in $(seq 1 $N)
do
jsub "${UCSD_ROOT_PATH}/src/NB_Greedy_library/CV/crossValidation.nbs $i ${UCSD_ROOT_PATH} ${DATA_FILE} ${UCSD_outname} $Gamma0 $Gamma1 $Prior"
done
