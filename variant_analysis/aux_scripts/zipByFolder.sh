#!/bin/bash

for dir in $(ls -d -1 */); do
	echo $dir
	#zip -r ${dir%?}_summary.zip $(find $dir -name '*_all.xlsx') ${dir}fasta
done
