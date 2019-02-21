#!/bin/bash

#PBS -q pdrineas
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=30:00:00

# set name of job
#PBS -N mysonMLGWASJob

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# #PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environments
#PBS -V

# start job from the directory it was submitted
#cd /scratch/brown/mcburch/GWASPipeline/MLGWAS/DIMPipeline_Alpha

cd "$PBS_O_WORKDIR" || exit $?

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
IFS=' ' read -r -a typelist <<< "$type"
echo ${typelist[@]}

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
      echo "No type specified so doing Ridge by default..."
      # usage
      # type="Ridge"
      typelist=( "Ridge" )
fi

if [[ $type == "all" ]];
then
    typelist=( "Ridge" "LDA" "SVR" "QDA" "PolyReg" "Lasso" "Elastic" "kSVM" "RidgeSVM" "RFESVM" "RandomForest" )
    # typelist=( "Ridge" "LDA" "SVR" "QDA" "PolyReg" "Lasso" "Elastic" "RidgeSVM" "RandomForest" )
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

echo $dataset
echo ${dataset##*/}

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
pvals=( 100 400 800 1200 1600 2000 )
# pvals=( 100 500 1000 )
# ratios=( "1:1" "1:1.25" "1:1.5" "1:2" "1:2.25" "1:2.5" "1:2.75" "1:3" )
ratios=( "1:1")
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
    # for different algorithm types...
    for type in "${typelist[@]}"
    do
        # for different selections of markers...
        for value in "${pvals[@]}"
        do
            if [ -d "${prefix}_${type}_${value}" ]; then
                rm -r "${prefix}_${type}_${value}/"
                mkdir "${prefix}_${type}_${value}/"
            else
                mkdir "${prefix}_${type}_${value}/"
            fi

            filename="${prefix}_SNPs_top_${value}.txt"
            # extract top 'k' p values
            # the 'top.assoc' file is found from the 'final_associations' step in the GWAS pipeline
            awk -v var="$value" 'FNR < var+2 && FNR !=1 {print $2}' "/depot/pdrineas/data/DIMs/Repo/${assocfile}" > $filename

            echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
            ${PLINK2PATH} --bfile "/depot/pdrineas/data/DIMs/Repo/${dataset##*/}" \
                                            --extract $filename \
                                            --make-bed \
                                            --out "${dataset##*/}_extracted_${value}"
            echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'


            echo 'Splitting data into train and test sets...'
            # filename="${dataset}_${type}_extracted_${value}_stats_${dt}.out"

            # split between train and test
            bash bose_splitdata.sh --cases 1500 --controls 1500 "${dataset##*/}_extracted_${value}"

            mv "${dataset##*/}_extracted_${value}"* "${filename}" "${prefix}_${type}_${value}/"

            dirPrefix=`readlink -e ${prefix}_${type}_${value}`

            trfilename="${dirPrefix}/${dataset##*/}_extracted_${value}_trainset"
            tefilename="${dirPrefix}/${dataset##*/}_extracted_${value}_testset"
            # run the analytics
            python DIMdetect2.py -tr $trfilename -te $tefilename -pf $prefix -cl $type -cv $cv -pval $value
            # mv $filename $outdir


        done
    done
else
    # for different algorithm types...
    for type in "${typelist[@]}"
    do
        # for different selections of markers...
        for value in "${pvals[@]}"
        do
            echo 'Splitting data into train and test sets...'
            # filename="${dataset}_${type}_extracted_${value}_stats_${dt}.out"

            # split between train and test
            bash bose_splitdata.sh --cases 1500 --controls 1500 "${dataset##*/}_extracted_${value}"

            mv "${dataset##*/}_extracted_${value}"* "${filename}" "${prefix}_${type}_${value}/"

            dirPrefix=`readlink -e ${prefix}_${type}_${value}`

            trfilename="${dirPrefix}/${dataset##*/}_extracted_${value}_trainset"
            tefilename="${dirPrefix}/${dataset##*/}_extracted_${value}_testset"
            # run the analytics
            python DIMdetect2.py -tr $trfilename -te $tefilename -pf $prefix -cl $type -cv $cv -pval $value
            # mv $filename $outdir

        done
    done
fi

python createTable.py --pvalList "${pvals[@]}" --typeList "${typelist[@]}" --data "${prefix}"

# output the time taken
echo run time is $(expr `date +%s` - $start_time) s

#END
