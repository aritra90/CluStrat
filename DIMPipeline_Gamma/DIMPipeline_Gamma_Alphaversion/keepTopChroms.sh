#!/bin/bash

# First run this: >>  awk '{print $1}' alpha_PRK_top.assoc | sort -g | cat | uniq -c | sort -g | tac > alpha_PRK_top_chromsort.txt
# Second run this: >> ./keepTopChroms.sh 4 alpha_PRK_top_chromsort.txt alpha_PRK_top.assoc

# $1 is the number of top chromosomes to keep
# $2 is the descending sorted chromosome file 
# $3 is the original 'assoc' file prefix

awk 'NR==1' ${3} > "top${1}chroms.assoc"

# Keep top 4 chromosomes 
chromArr=($(awk -v var="$1" 'FNR < var+1 {print $2}' $2)) 

echo "${chromArr[@]}"

for chromosome in "${chromArr[@]}"
do
	awk -v var="$chromosome" '$1==var' ${3} >> "top${1}chroms.assoc"
done

