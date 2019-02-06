#!/bin/bash
#PBS -q pdrineas
#PBS -l nodes=1:ppn=4,naccesspolicy=shared
#PBS -l walltime=05:00:00

# set name of job
#PBS -N mysonGWASJob

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environment
#PBS -V

cd /scratch/brown/mcburch/GWASPipeline/MLGWAS/DIMPipeline_Beta

module load texlive

dataset=$1
echo ${dataset}

echo '#######################STARTING PREPROCESSING##############################'
echo ' '
bash array_illumina_psycharray.sh "${dataset}" yes yes yes 2048 > preprocessing.out


echo '#######################STARTING QUALITY CONTROL##############################'
echo ' '
bash dataset_qc.sh --geno1 0.05 --geno2 0.02 --ibrc 0.2  "${dataset}_array/${dataset}" mixed 2048 > qc.out



echo '#######################STARTING IBD CHECK##############################'
echo ' '
bash ibd_check.sh "${dataset}_qc/${dataset}_qcind" emptyfiles/outliers 8 2048 --maf 0.05 --hwe 0.001 --snpmissing 0.02 --remainingSNPS 150000 --relatednessThresh 0.15 > ibd.out



echo '#######################STARTING PCA##############################'
echo ' '
bash /scratch/brown/mcburch/GWASPipeline/pca2.sh "${dataset}_qc/${dataset}"_qcind_qcsnp emptyfiles/outliers no 8 2048 --maf 0.05 --hwe 0.001 --snpmissing 0.02 --remainingSNPS 150000 > pca.out


echo '#######################STARTING FINAL ASSOCIATIONS##############################'
echo ' '
bash /scratch/brown/mcburch/GWASPipeline/final_associations.sh "${dataset}_qc/${dataset}_qcind_qcsnp" "${dataset}_qcind_qcsnp_pca_outliers/pca.covar" 1 simple "${dataset}_qcind_ibd_outliers/ibdout.ind" 8 2048 --pvalThresh 0.003 > finalassoc.out

#END
