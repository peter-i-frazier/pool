#!/bin/sh
UCSD_ROOT_PATH="/fs/home/jw865/peptide-catalysis"
UCSD_outname="sfp"
# UCSD_datafile="${UCSD_ROOT_PATH}/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv"
rm -R ${UCSD_ROOT_PATH}/temp/loocv_find_param/${UCSD_outname}
mkdir ${UCSD_ROOT_PATH}/temp/loocv_find_param/${UCSD_outname}

N=210
for i in $(seq 1 $N)
do
jsub "${UCSD_ROOT_PATH}/src/NB_Greedy_library/loocv_find_param/loocv_find_param.nbs ${UCSD_ROOT_PATH} $i ${UCSD_outname}" 
done

