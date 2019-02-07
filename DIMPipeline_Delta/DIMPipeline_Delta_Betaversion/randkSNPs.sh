#!/bin/bash

# $1 needs to be the 'assoc' file prefix 

pvals=( 100 400 800 1200 1600 2000 )

# for different selections of markers...
for value in "${pvals[@]}"
do
	awk 'NR==1' HeaderAssoc.txt > ${1}_rand_${value}.assoc
	awk '{if(NR>1)print}' ${1}.assoc | shuf -n ${value} >> ${1}_rand_${value}.assoc
done


