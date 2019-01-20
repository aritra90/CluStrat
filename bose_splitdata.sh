#!/bin/bash
#Author: Aritra Bose

source /depot/pdrineas/data/gwas_scripts/software_paths.conf

int_re='^[0-9]+$'

PARAMS=""
# Default values
subsample=""
num_cases="300"
num_controls="300"

while (( "$#" )); do
  case "$1" in
    --cases)
      num_cases=$2
      if ! [[ "$(($num_cases))" =~ $int_re ]] ; then
          echo "Number of cases needs to be a positive integer.";
          exit 1
      fi
      shift 2
      ;;
    --controls)
      num_controls=$2
      if ! [[ "$(($num_controls))" =~ $int_re ]] ; then
          echo "Number of controls needs to be a positive integer.";
          exit 1
      fi
      shift 2
      ;;
    # --subsample)
      # echo "$(($2))"
      # if ! [[ "$(($subsample))" =~ $int_re ]] ; then
      #     echo "Subsampling value needs to be a positive integer.";
      #     exit 1
      # fi

      # subsample=$2
      # arrIN=(${subsample//:/ })
      # num_controls=$(echo "300*${arrIN[1]}" | bc -l)
      # num_controls=$(echo "($num_controls+0.5)/1" | bc)
      # shift 2
      # ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

######---------------------------------------------------------
######    Input arguments
dataset=`readlink -e $1.bed | sed 's/\.bed//'`
######---------------------------------------------------------


##################
###split into Test and Train
##################

echo `awk '{if($6==2) print $1,$2}' ${dataset##*/}.fam > ${dataset##*/}_cases.txt`

echo `awk '{if($6==1) print $1,$2}' ${dataset##*/}.fam > ${dataset##*/}_controls.txt`
##################
#20 and 250 here is a placeholder, can use user input as well
##################

shuf -n 20 ${dataset##*/}_cases.txt > ${dataset##*/}_cases_test.txt
shuf -n 20 ${dataset##*/}_controls.txt > ${dataset##*/}_controls_test.txt

cat ${dataset##*/}_cases_test.txt > ${dataset##*/}_testsamples.txt
cat ${dataset##*/}_controls_test.txt >> ${dataset##*/}_testsamples.txt

# Storing the difference between the test cases/controls as the trainset
grep -Fvxf ${dataset##*/}_cases_test.txt ${dataset##*/}_cases.txt > ${dataset##*/}_cases_train.txt
grep -Fvxf ${dataset##*/}_controls_test.txt ${dataset##*/}_controls.txt > ${dataset##*/}_controls_train.txt

# subsampling the training set
if [ -z "$num_cases" ]
then
      cat ${dataset##*/}_cases_train.txt > ${dataset##*/}_trainsamples.txt
      cat ${dataset##*/}_controls_train.txt >> ${dataset##*/}_trainsamples.txt
else
      shuf -n "$(($num_cases))" ${dataset##*/}_cases_train.txt > ${dataset##*/}_subcases_train.txt
      shuf -n "$(($num_controls))" ${dataset##*/}_controls_train.txt > ${dataset##*/}_subcontrols_train.txt

      cat ${dataset##*/}_subcases_train.txt > ${dataset##*/}_trainsamples.txt
      cat ${dataset##*/}_subcontrols_train.txt >> ${dataset##*/}_trainsamples.txt
fi

echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset##*/} \
                                --allow-no-sex \
                                --keep  ${dataset##*/}_trainsamples.txt \
                                --make-bed \
                                --out ${dataset##*/}_trainset
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset##*/} \
                                --allow-no-sex \
                                --keep  ${dataset##*/}_testsamples.txt \
                                --make-bed \
                                --out ${dataset##*/}_testset
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
