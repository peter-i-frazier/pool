#!/bin/sh
UCSD_ROOT_PATH="/fs/home/jw865/peptide-catalysis"
N=2573
UCSD_datafile="${UCSD_ROOT_PATH}/data/all_TS_data.csv"
result_path="${UCSD_ROOT_PATH}/temp/get_prob"
rm -R ${result_path}
mkdir ${result_path}

for i in $(seq 1 $N)
do
    jsub "${UCSD_ROOT_PATH}/src/NB_Greedy_library/get_prob/get_prob_of_all_data.nbs ${UCSD_ROOT_PATH} ${UCSD_datafile} $i ${result_path}" -mfail -email charleswang304@gmail.com -xhost sys_pf
done


