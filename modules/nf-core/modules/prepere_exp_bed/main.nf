
process PREPERE_EXP_BED {
  label 'process_low'
  tag {condition}
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/eqtl.img"
      
  } else {
      container "quay.io/biocontainers/multiqc:1.10.1--py_0"
  }


  input:
    tuple(val(condition),path(mapping_file),path(expression_file),path(phenotype_pcs) )
    each path(annotation_file)
    each path(genotype_pcs)

  output:
    tuple(val(condition),path("Expression_Data.bed.gz"),path('Covariates.tsv'), emit: exp_bed)

  script:

    if(params.sample_covariates==''){
      sample_covar =''
    }else{
      sample_covar ="--sample_covariates ${params.sample_covariates}"
    }
    """
      echo ${condition}
      prepere_bed.py --annotation_file ${annotation_file} --mapping_file ${mapping_file} --expression_file ${expression_file}
      prepere_covariates_file.py --genotype_pcs ${genotype_pcs} --phenotype_pcs ${phenotype_pcs} ${sample_covar} --sample_mapping ${mapping_file}
    """
}