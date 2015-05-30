#!/bin/bash
# ./crop_part.sh dest_folder ori_folder filename num_col
convert $2/$3.png -crop $4x1@ $1/$3_%d.png

