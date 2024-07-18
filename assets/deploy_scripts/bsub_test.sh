#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE="$@"
export RUN_ID="${PWD##*/}"

# export SINGULARITY_CACHEDIR='/software/hgi/containers/QTLight'
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export NXF_OPTS="-Xms5G -Xmx5G"
export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp

sample="$RUN_ID.QTLight"
echo -e "\nSubmitting QTLight (https://github.com/wtsi-hgi/QTLight) in test mode withsample OneK1k dataset"
bsub -R'select[mem>4000] rusage[mem=4000]' -J QTLight_test -n 1 -M 4000 -o QTLight_test.o -e QTLight_test.e -q normal bash $SCRIPT_DIR/../../assets/deploy_scripts/nohup_start_nextflow_lsf_test.sh $INPUT_FILE
echo "Submitted job can be killed with: bkill -J QTLight_test"