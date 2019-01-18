#!/bin/bash

#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=33:00:00

# set name of job
#PBS -N MLGWASJob

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environments
#PBS -V

# start job from the directory it was submitted
cd /home/bose6/PopGen/DIM/PRK

source /depot/pdrineas/data/gwas_scripts/software_paths.conf

function usage() {
    echo ""
    echo "Usage: qsub -v dataset_prefix=dataset_prefix,extract=yes/no,type=LDA/SVM/Ridge...,cv=kfold/loocv,help=yes/no new_GWAS_ML_script.sh"
    echo "See MLGWAS.rmd for all classifier type options..."
    echo ""
}

# # Default values
# extract=""
# type="LDA"
# cv="kfold"
# dataset_prefix=""
# help=""

if [[ $help == "yes" ]];
then
    usage
    exit 1
fi

if [ -z "$extract" ]
then
      echo "No extraction specified so doing so as a precaution..."
      # usage
      extract="yes"
fi

if [ -z "$type" ]
then
      echo "No type specified so doing LDA by default..."
      # usage
      type="LDA"
fi

if [ -z "$cv" ]
then
      echo "No cross validation specified so doing k-fold by default..."
      # usage
      cv="kfold"
fi

echo $cv > "crossvalType.txt"

######---------------------------------------------------------
######    Input arguments
dataset=`readlink -e $dataset_prefix.bed | sed 's/\.bed//'`

if [ -z "$dataset" ]
then
      echo "Please enter a dataset prefix..."
      usage
      exit 1
fi

# load necessary modules
module load gcc
module load bioinfo
module load plink
module load texlive

# variable for different values of SNPs kept
value=0
filename=""
trfilename=""
tefilename=""
# start timer
start_time=`date +%s`
dt=$(date +"%Y-%m-%d_%T");
outdir=""
pvals=( 100 500 1000 1500 2000 )
assocfile=""
prefix=""
######---------------------------------------------------------
# echo ${dataset##*/}
if [[ ${dataset##*/} == "SZHCGSRRS" ]];
then
    assocfile="SCZ_top.assoc"
    prefix="SCZ"
elif [[ ${dataset##*/} == "PKHCGSRRS" ]];
then
    assocfile="PRK_top.assoc"
    prefix="PRK"
else
    echo "No associations file found..."
    exit 1
fi


if [[ $extract == "yes" ]];
then
    # for different selections of markers...
    for value in "${pvals[@]}"
    do
        filename="${prefix}_SNPs_top_${value}.txt"
        # extract top 'k' p values
        # the 'top.assoc' file is found from the 'final_associations' step in the GWAS pipeline
        awk -v var="$value" 'FNR < var+2 && FNR !=1 {print $2}' $assocfile > $filename

        echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
        ${PLINK2PATH} --bfile ${dataset##*/} \
                                        --extract $filename \
                                        --make-bed \
                                        --out ${dataset##*/}_extracted_$value
        echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

        # outdir="${dataset}_${type}_${value}_${dt}"
        # mkdir -p $outdir
        echo '(0) Splitting data into train and test sets...'
        filename="${dataset}_${type}_extracted_${value}_stats_${dt}.out"
        # split between train and test
        bash bose_splitdata.sh --subsample 250 "${dataset}_extracted_${value}"
        trfilename="${dataset}_extracted_${value}_trainset"
        tefilename="${dataset}_extracted_${value}_testset"
        # run the analytics
        python DIMdetect2.py -tr $trfilename -te $tefilename -pf $prefix -cl $type -cv $cv -pval $value | tee $filename
        # mv $filename $outdir
    done
else
    # for different selections of markers...
    for value in "${pvals[@]}"
    do
        # outdir="${dataset}_${type}_${value}_${dt}"
        # mkdir -p $outdir
        echo '(1) Splitting data into train and test sets...'
        filename="${dataset}_${type}_extracted_${value}_stats_${dt}.out"
        # split between train and test
        bash bose_splitdata.sh --subsample 250 "${dataset}_extracted_${value}"
        trfilename="${dataset}_extracted_${value}_trainset"
        tefilename="${dataset}_extracted_${value}_testset"
        # run the analytics
        python DIMdetect2.py -tr $trfilename -te $tefilename -pf $prefix  -cl $type -cv $cv -pval $value | tee $filename
        # mv $filename $outdir
    done
fi

# output the time taken
echo run time is $(expr `date +%s` - $start_time) s

#END
