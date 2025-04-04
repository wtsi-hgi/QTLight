#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE="$@"
export RUN_ID="${PWD##*/}"
export NXF_OPTS="-Xms35G -Xmx35G"

export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

sample="$RUN_ID"
bsub -R'select[mem>75000] rusage[mem=75000]' -J $sample -n 1 -M 75000 -o $sample.o -e $sample.e -q $QUEUE bash $SCRIPT_DIR/../../assets/deploy_scripts/nohup_start_nextflow_lsf.sh $INPUT_FILE
echo "Submitted job can be killed with: bkill -J $sample"