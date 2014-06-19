ROOT_PATH="/fs/home/jw865/peptide-catalysis"
UCSD_outname="AcpS_specific"

rm ${ROOT_PATH}/src/NB_Greedy_library/find_param/auc

numdirs=(${ROOT_PATH}/temp/find_param/${UCSD_outname}/*)
numdirs=${#numdirs[@]}

for n in $(seq 1 $numdirs)
do
        folder_name=param_$n
        f=${ROOT_PATH}/temp/find_param/${UCSD_outname}/${folder_name}
        cat $f/auc >> ${ROOT_PATH}/src/NB_Greedy_library/find_param/auc
done


