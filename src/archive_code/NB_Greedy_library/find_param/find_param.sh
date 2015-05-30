#!/bin/sh
UCSD_ROOT_PATH="/fs/home/jw865/peptide-catalysis"
UCSD_outname="AcpS"
N=378
num_fold=260

UCSD_datafile="${UCSD_ROOT_PATH}/data/all_TS_data.csv"
result_path="${UCSD_ROOT_PATH}/temp/find_param/allTS_${num_fold}fold_${UCSD_outname}"
rm -R ${result_path}
mkdir ${result_path}

for i in $(seq 1 $N)
do
    jsub "${UCSD_ROOT_PATH}/src/NB_Greedy_library/find_param/find_param.nbs ${UCSD_ROOT_PATH} ${UCSD_datafile} $i ${UCSD_outname} ${result_path} ${num_fold}" -mfail -email charleswang304@gmail.com -xhost sys_pf -queue long
done

