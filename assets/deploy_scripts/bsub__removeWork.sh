#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"
sample="$RUN_ID.yascp"
echo "Cleaning the work directory (https://github.com/wtsi-hgi/yascp) with input file $INPUT_FILE by using '-entry WORK_DIR_REMOVAL --remove_work_dir' "
bsub -R'select[mem>4000] rusage[mem=4000]' -J $sample -n 1 -M 4000 -o $sample.o -e $sample.e -q long bash /software/hgi/pipelines/yascp/assets/deploy_scripts/nohup_start_nextflow_lsf__removeWork.sh $INPUT_FILE
echo "Submitted job can be killed with: bkill -J $sample"