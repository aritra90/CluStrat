#!/bin/bash

#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=30:00:00

# set name of job
#PBS -N GEMMA_Job

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# #PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environments
#PBS -V

# start job from the directory it was submitted
cd "$PBS_O_WORKDIR" || exit $?

source /depot/pdrineas/data/gwas_scripts/software_paths.conf

dataset_path=$1

# https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf

# 1) Basic QC
/depot/pdrineas/data/gwas_software/plink2/plink --bfile "${dataset_path}" --maf 0.01 --geno 0.1 --hwe 0.001 --make-bed --out dummy_BED_qc
echo ' '

# 2) Estimate relatedness matrix (gk = 1 is "center relatedness" matrix)
/depot/pdrineas/data/gwas_software/gemma/gemma-0.98.1 -bfile dummy_BED_qc -gk 1 -o dummy_cntr_relatedness
echo ' '

# 3) Eigen-Decomp of relatedness matrix (eigenvectors U and eigenvalues D) 
/depot/pdrineas/data/gwas_software/gemma/gemma-0.98.1 -bfile dummy_BED_qc -k output/dummy_cntr_relatedness.cXX.txt -eigen -o dummy_eigdecomp 
echo ' '

# 4) Assoc. tests w/ LMM 
/depot/pdrineas/data/gwas_software/gemma/gemma-0.98.1 -bfile dummy_BED_qc -k output/dummy_cntr_relatedness.cXX.txt -lmm 4 -o dummy_LMM
echo ' '


# Example run: ./gemma_script.sh simdata_TGP_1440957
# End of script
