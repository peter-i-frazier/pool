ROOT_PATH="/fs/home/jw865/peptide-catalysis"

rm ${ROOT_PATH}/output/cv

numdirs=(${ROOT_PATH}/temp/cv/*/)
numdirs=${#numdirs[@]}

for n in $(seq 1 $numdirs)
do
        folder_name=cv_$n
        f=${ROOT_PATH}/temp/cv/${folder_name}
        cat $f/prob >> ${ROOT_PATH}/output/cv
done


