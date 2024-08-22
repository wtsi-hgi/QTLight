#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE="$@"
export RUN_ID="${PWD##*/}"
export NXF_OPTS="-Xms5G -Xmx5G"
export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

sample="$RUN_ID"
echo -e "\n Submitting eqtl (https://github.com/wtsi-hgi/eqtl) with input file $INPUT_FILE"
bsub -R'select[mem>3000] rusage[mem=3000]' -J $sample -n 1 -M 3000 -o $sample.o -e $sample.e -q $QUEUE bash $SCRIPT_DIR/../../assets/deploy_scripts/nohup_start_nextflow_lsf.sh $INPUT_FILE
echo "Submitted job can be killed with: bkill -J $sample"