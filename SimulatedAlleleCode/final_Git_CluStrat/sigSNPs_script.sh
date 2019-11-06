#!/bin/bash
#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

# set name of job
#PBS -N myson_PRK_CluStrat

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environment
#PBS -V

cd "$PBS_O_WORKDIR" || exit $?

source /depot/pdrineas/data/gwas_scripts/software_paths.conf

#echo '####################### STARTING QUALITY CONTROL ##############################'
#echo ' '
#/depot/pdrineas/data/gwas_software/plink2/plink --bfile /depot/pdrineas/data/DIMs/Repo/PKHCGSRRS --maf 0.01 --geno 0.1 --hwe 0.001 \
#                                                --allow-no-sex \
#						--indep-pairwise 1000 50 0.2 \
#                                                --make-bed --out PKHCGSRRS_qc

#/depot/pdrineas/data/gwas_software/plink2/plink --bfile PKHCGSRRS_qc --extract PKHCGSRRS_qc.prune.in --make-bed --out PKHCGSRRS_pruned

#echo '####################### STARTING EMMAX ##############################'
#echo ' '
#./emmax_script.sh PKHCGSRRS_pruned
 
#echo '####################### STARTING GEMMA ##############################'
#echo ' '
#./gemma_script.sh PKHCGSRRS_pruned

echo '####################### STARTING CLUSTRAT ##############################'
echo ' '
#python3 realdata_clustrat.py 1kgrs_chr22_subset
#python3 CluStrat_wrapper.py --dir ../PKHCGSRRS_pruned
python3 CluStrat_wrapper.py --sim 1
