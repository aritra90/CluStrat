#!/bin/bash

#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=30:00:00

# set name of job
#PBS -N EMMAX_Job

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

# https://genome.sph.umich.edu/wiki/EMMAX
# http://genetics.cs.ucla.edu/emmax/install.html

# 1) Convert to tped format 
/depot/pdrineas/data/gwas_software/plink2/plink --bfile "${dataset_path}" --recode12 --output-missing-genotype 0 --transpose --out dummy_TPED

# 1b) Basic QC
/depot/pdrineas/data/gwas_software/plink2/plink --tfile dummy_TPED --maf 0.01 --geno 0.1 --hwe 0.001 --make-bed --out dummy_BED_qc

# 1c) Convert back 
/depot/pdrineas/data/gwas_software/plink2/plink --bfile dummy_BED_qc --recode12 --output-missing-genotype 0 --transpose --out dummy_TPED_qc

# 2) Compute kinship/relatedness matrix 
/depot/pdrineas/data/gwas_software/emmax/emmax-kin -v -h -d 10 dummy_TPED_qc

# 3) Create phenotype file according to documentation
awk '{print $1, $2, $6}' dummy_TPED_qc.tfam > dummy_TPED_qc.pheno

# 4) Run association using mixed model
/depot/pdrineas/data/gwas_software/emmax/emmax -v -d 10 -t dummy_TPED_qc -p dummy_TPED_qc.pheno -k dummy_TPED_qc.hBN.kinf -o dummy_outfile


# Example run: ./emmax_script.sh simdata_TGP_1440957
# End of script
