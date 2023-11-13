
export TEMP=$PWD/work
export TMP_DIR=$PWD/work
nohup nextflow run yascp -profile test_full -c /scratch/cellfunc/mo246/configs/conf/leicester.config --nf_ci_loc $PWD -resume > nextflow.nohup.log
