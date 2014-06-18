#!/bin/sh
ROOT_PATH = "/fs/home/jw865/peptide-catalysis"
rm -R "${ROOT_PATH}/temp/*/"

N=$(($(wc -l < "${ROOT_PATH}/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv") - 1))

for i in $(seq 1 $N)
do
	jsub "crossValidation.nbs $i ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}" 
done




