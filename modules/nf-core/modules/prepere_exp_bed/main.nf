
process PREPERE_EXP_BED {
  label 'process_low'
  tag {condition}
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/eqtl.img"
      
  } else {
      container "quay.io/biocontainers/multiqc:1.10.1--py_0"
  }


  input:
    tuple(val(condition),path(mapping_file),path(expression_file))
    path(annotation_file)

  output:
    tuple(val(condition),path("Expression_Data.bed.gz"), emit: exp_bed)

  script:

    """
      echo ${condition}
      prepere_bed.py --annotation_file ${annotation_file} --mapping_file ${mapping_file} --expression_file ${expression_file}
    """
}