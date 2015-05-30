#!/bin/sh
rm -R ${UCSD_ROOT_PATH}/temp/CV/*

N=$(($(wc -w < ${ROOT_PATH}/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv) - 1))
echo $N

for g in (seq 1 $UCSD_num_gamma)
    for i in $(seq 1 $N)
    do
        jsub "crossValidation.nbs $g $i"
    done
done




