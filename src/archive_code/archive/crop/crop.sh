#!/bin/bash
# ./crop.sh dest_folder ori_folder filename num_row num_col num_row-1
convert $2/$3.png -crop 1x$4@ $1/temp_$3_%d.png
for i in `seq 0 $6`;
do
    ./crop_part.sh $1 $1 temp_$3_$i $5
done
