process SPLIT_AGGREGATION_ADATA {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

  input:
    path(adata) // lists input files per donor
    val(agg_columns)

  output:
    path("*__split.h5ad", emit:split_phenotypes)

  """
    split_adata_per_condition.py --agg_columns '${agg_columns}' -h5ad ${adata}
  """
}


process AGGREGATE_UMI_COUNTS {
  publishDir  path: "${params.outdir}/aggregated_counts/${sanitized_columns}",mode: "${params.copy_mode}",
              overwrite: "true"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

  input:
    path(adata) // lists input files per donor
    val(agg_columns)
    val(gt_id_column)
    val(sample_column)
    val(n_cells_min)
    val(n_donors_min)

  output:
    path("*phenotype_file.tsv", emit:phenotype_file)
    path("*genotype_phenotype_mapping.tsv", emit:genotype_phenotype_mapping)
    tuple val(sanitized_columns), path("*phenotype_file.tsv"),  path("*genotype_phenotype_mapping.tsv"), emit:phenotype_genotype_file
    path('*.tsv')

  script:
    sanitized_columns = adata.getName().replaceAll(/[^a-zA-Z0-9]/, '_').replaceAll(/\.h5ad$/, '')
    """
      echo ${sanitized_columns}
      aggregate_sc_data.py --agg_columns '${agg_columns}' --gt_id_column '${gt_id_column}' --sample_column '${sample_column}' --n_cells ${n_cells_min} -n_individ ${n_donors_min} -h5ad ${adata} --method ${params.aggregation_method}
    """
}