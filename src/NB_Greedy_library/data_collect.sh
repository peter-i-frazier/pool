rm /fs/home/py75/Documents/peptide/output/output

numdirs=(/fs/home/py75/Documents/peptide/temp/*/)
numdirs=${#numdirs[@]}

for n in $(seq 0 $(($numdirs -1)))
do
        folder_name=temp_$n
        f=/fs/home/py75/Documents/peptide/temp/$folder_name
        cat $f/output >> /fs/home/py75/Documents/peptide/output/output
done


