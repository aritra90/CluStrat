#!/bin/bash

# start timer
start_time=`date +%s`

PARAMS=""
data="PRK"

while (( "$#" )); do
  case "$1" in
    --data)
       data=$2
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

echo $data

files=( "${data}_caseaccuracy.csv" "${data}_controlaccuracy.csv" "${data}_f1score.csv" "${data}_trainingaccuracy.csv" "${data}_testingaccuracy.csv" )

for thefile in "${files[@]}"
do
	echo $thefile
	python spagplots.py --file "${thefile}"
done

# output the time taken
echo run time is $(expr `date +%s` - $start_time) s
