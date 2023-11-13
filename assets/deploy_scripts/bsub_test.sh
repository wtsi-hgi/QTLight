#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"

# export SINGULARITY_CACHEDIR='/software/hgi/containers/yascp'

export NXF_OPTS="-Xms5G -Xmx5G"
export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp

echo press ENTER to not fetch containers, otherwise provide writable path:
read varname

if "$varname"==''
    then
        export NXF_SINGULARITY_CACHEDIR='/software/hgi/containers/yascp'
        export SINGULARITY_DISABLE_CACHE=0
    else
        echo eqtl Will fetch the containers and place them in $varname
        export NXF_SINGULARITY_CACHEDIR=$varname
fi

sample="$RUN_ID.yascp"
echo -e "\nSubmitting yascp (https://github.com/wtsi-hgi/yascp) in test mode withsample OneK1k dataset"
bsub -R'select[mem>4000] rusage[mem=4000]' -J yascp_test -n 1 -M 4000 -o yascp_test.o -e yascp_test.e -q normal bash /software/hgi/pipelines/yascp/assets/deploy_scripts/nohup_start_nextflow_lsf_test.sh
echo "Submitted job can be killed with: bkill -J yascp_test"