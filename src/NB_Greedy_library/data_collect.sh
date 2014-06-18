ROOT_PATH = "/fs/home/jw865/peptide-catalysis"

rm "${ROOT_PATH}/output/output"

numdirs=("${ROOT_PATH}/temp/*/")
numdirs=${#numdirs[@]}

for n in $(seq 1 $numdirs)
do
        folder_name=temp_$n
        f="${ROOT_PATH}/temp/${folder_name}"
        cat $f/output >> "${ROOT_PATH}/output/output"
done


