#!/usr/bin/env bash
INPUT_FILE=$1
dt=`date +"%Y_%m_%d_%T"`
cp nextflow.nohup.log ./nextflow.nohup_$dt.log2 || echo 'first time running'
# activate Nextflow conda env

# clean up previous run files
rm -f *.log
rm -f nextflow.nohup.PID.txt 

# start Nextflow in background:

CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
# export RUN_ID="${parentdir##*/}"
export RUN_ID="${PWD##*/}"
mkdir $PWD/work || echo 'exists'
mkdir $PWD/work/tmp || echo 'exists'
echo $RUN_ID | nextflow run /software/hgi/pipelines/yascp -profile sanger,test --nf_ci_loc $PWD -resume > nextflow.nohup.log 2>&1 & 

# get process PID 
sleep 1 && export PID=$(pgrep -f "\\-\\-nf_ci_loc $RUN_DIR")
echo $PID > nextflow.nohup.PID.txt
echo "Nextflow PID is $PID (saved in ./nextflow.nohup.PID.txt)" 
echo kill with \"kill $PID\"
echo "check logs files nextflow.nohup.log and .nextflow.log"
