ROOT_PATH="/fs/home/jw865/peptide-catalysis"

rm ${ROOT_PATH}/output/auc

numdirs=(${ROOT_PATH}/temp/find_param/*/)
numdirs=${#numdirs[@]}

for n in $(seq 1 $numdirs)
do
        folder_name=param_$n
        f=${ROOT_PATH}/temp/find_param/${folder_name}
        cat $f/auc >> ${ROOT_PATH}/output/auc
done


