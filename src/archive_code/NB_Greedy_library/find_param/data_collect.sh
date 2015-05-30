ROOT_PATH="/fs/home/jw865/peptide-catalysis"
UCSD_outname="PfAcpH"
result_path="${ROOT_PATH}/temp/find_param/allTS_500fold_${UCSD_outname}"
auc_collection="${ROOT_PATH}/src/NB_Greedy_library/find_param/auc"

rm ${auc_collection}

numdirs=(${result_path}/*)
numdirs=${#numdirs[@]}

for n in $(seq 1 $numdirs)
do
        folder_name=param_$n
        f=${result_path}/${folder_name}
        cat $f/auc >> ${auc_collection}
done


