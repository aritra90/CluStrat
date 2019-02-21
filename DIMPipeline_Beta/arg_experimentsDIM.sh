#!/bin/bash

#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=30:00:00

# set name of job
#PBS -N mysonDIMExperiments

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# #PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environments
#PBS -V

# start job from the directory it was submitted
cd "$PBS_O_WORKDIR" || exit $?

if [ -d "Experiments" ]; then
	rm -r "Experiments/"
	mkdir "Experiments/"
else
	mkdir "Experiments/"
fi

expVal=$1

./clean.sh

# $1 is the amount of experiments
counter=0
while [ $counter -lt $expVal ]
do
	let counter++
	bash arg_new_GWAS_ML_script.sh /depot/pdrineas/data/DIMs/Repo/PKHCGSRRS kfold "Ridge LDA RidgeSVM"
	# sleep for approximately how long the job should take (with some wiggle room)
	# sleep 25m
	# save the tables as exp 1/2/3
	for f in *.csv
	do
    		mv "$f" "Exp_${counter}_$f"
	done
	mv *.csv Experiments/
	./clean.sh
done

# Average the tables
python avgTables.py --experiments "${expVal}"
