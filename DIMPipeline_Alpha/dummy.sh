#!/bin/bash

echo $#
CURR_DIR=$(pwd)
#set -x

#PBS -q pdrineas
#PBS -l nodes=1:ppn=2,naccesspolicy=shared
#PBS -l walltime=05:00:00

# set name of job
#PBS -N mysonMLGWASDummyJob

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m bea

# send mail to the following address
#PBS -M mcburch@purdue.edu

# use submission environment
#PBS -V

# start job from the directory it was submitted
cd /scratch/brown/mcburch/GWASPipeline/MLGWAS

# declare an array called array and define 3 vales
array=( one two three )
for i in "${array[@]}"
do
	echo $i
done

function usage() {
    echo ""
    echo "Usage: $0 dataset_prefix [--extract yes/no] [--type LDA/SVM/Ridge...] [--cv kfold/loo]"
    echo "See MLGWAS.rmd for all classifier type options..."
    echo ""
}

PARAMS=""
# Default values
extract=""
type="LDA"
cv="kfold"
ratio="1:1.5"

while (( "$#" )); do
  case "$1" in
    --help)
          usage
          exit 1
          ;;
    --extract)
        #echo "$2"
        extract=$2
        if [ $extract == "yes" ];
        then
            echo "Extracting top SNPs..."
        # else
        #     usage
        #     echo "Not a proper flag!"
        #     exit 1
        fi
        shift 2
        ;;
    --type)
        #echo "$2"
        type=$2
        shift 2
        ;;
    --cv)
        #echo "$2"
        cv=$2
        shift 2
        ;;
    --ratio)
	echo "$2"
	ratio=$2
        arrIN=(${ratio//:/ })
        echo ${arrIN[0]}
        echo ${arrIN[1]}
        thisVar=$(echo "scale=5;300*${arrIN[1]}" | bc -l)
        thisVar=$(echo "($thisVar+0.5)/1" | bc)
        #echo $thisVar
        #echo "$((300*${arrIN[1]}))"
	shift 2
        ;;
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

thenum="200"
ratios=( "1:1" "1:1.25" "1:1.5" "1:2" "1:2.25" "1:2.5" "1:2.75" "1:3" )

for value in "${ratios[@]}"
    do
        arrIN=(${value//:/ })
        echo " "
        echo ${arrIN[0]}
        echo ${arrIN[1]}
        thisVar=$(echo "scale=5;$thenum*${arrIN[1]}" | bc -l)
        thisVar=$(echo "($thisVar+0.5)/1" | bc)
        echo $thisVar
        echo " "
    done

echo "${CURR_DIR}/zzzz"

if [ -d "zzzz" ]; then
    echo "HEY WHATS UP"
    rmdir "zzzz/"
fi


