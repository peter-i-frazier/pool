#!/bin/sh

rm -R /fs/home/py75/Documents/peptide/temp/*/

N=$(($(wc -l < /fs/home/py75/Documents/peptide/data/Data_10.csv) - 1))

for i in $(seq 1 $N)
do
	jsub "crossValidation.nbs $i ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}" 
done




