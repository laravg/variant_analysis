#!/bin/bash

for dir_slash in */; do
    dir=${dir_slash::-1}
    mkdir ${dir}/summary_corrected_bc_wrong_strains
done

# dirs=$(ls -d -1 */ | awk '{print substr($0, 1, length($0)-1)}' | egrep -v "^genomes$" )

# echo "$dirs" | while read dir; do
#     folder=$(echo "$dir" | sed -e 's/.*.\.txt_//')
#     echo "Creating $folder"
#     echo "Moving $dir to $folder"
#     mkdir $folder 2> /dev/null
#     mv $dir $folder/
# done
