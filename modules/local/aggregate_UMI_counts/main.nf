process SPLIT_AGGREGATION_ADATA {
    // label 'process_medium'
    memory { 
            sizeInGB = adata.size() / 1e9 * 0.5 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
      
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

  publishDir  path: "${params.outdir}/per_condition_adatas/${params.aggregation_subentry}",mode: "${params.copy_mode}",
              overwrite: "true"

  input:
    path(adata) // lists input files per donor
    val(agg_columns)

  output:
    path("*__split.h5ad", emit:split_phenotypes)
  script:
    if ("${params.aggregation_subentry}"==''){
        cond1 = ""
    }else{
        cond1 = " --condition '${params.aggregation_subentry}' "
    }
  """
    split_adata_per_condition.py --agg_columns '${agg_columns}' -h5ad ${adata} ${cond1}
  """
}


process CHUNK_ADATA_FOR_AGGREGATION {
    // label 'process_medium'
    memory { 
            sizeInGB = adata.size() / 1e9 * 0.5 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
      
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }


  input:
    path(adata) // lists input files per donor
    val(agg_columns)
    val(number_of_cunks)

  output:
    path("*__split.h5ad--*", emit:split_phenotypes)
  script:

  """
    split_adata_in_cunks.py --agg_columns '${agg_columns}' -h5ad ${adata} --number_of_cunks ${number_of_cunks}
  """
}


process ORGANISE_AGGREGATED_FILES{

  publishDir  path: "${params.outdir}/aggregated_counts/${sanitized_columns}",mode: "${params.copy_mode}",
              overwrite: "true"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }
    memory { 
            sizeInGB = adata.size() / 1e9 * 0.3 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
  input:
    tuple val(sanitized_columns), path(phenotype_files),  path(genotype_phenotype_mapping)

  output:
    path("clean_table.tsv", emit:phenotype_files_tsv) optional true


  script:
    
    """
      organise_channels.py -files "${phenotype_files}" -files2 "${genotype_phenotype_mapping}"
    """


}

process COMBINE_AGGREGATES {

    
    publishDir  path: "${params.outdir}/aggregated_counts/${sanitized_columns}",mode: "${params.copy_mode}",
              overwrite: "true"
    
    
    label 'process_medium'
    tag "$sanitized_columns"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

  input:
    tuple val(sanitized_columns), path(phenotype_file),  path(genotype_phenotype_mapping)
    val(n_min_individ)

  output:    
    tuple val(sanitized_columns), path("[!0-9]*___phenotype_file.tsv"),  path("[!0-9]*___genotype_phenotype_mapping.tsv"), emit:phenotype_genotype_file optional true
    path("[!0-9]*___phenotype_file.tsv", emit:phenotype_file) optional true
    path("[!0-9]*___genotype_phenotype_mapping.tsv", emit:genotype_phenotype_mapping) optional true    
    path("[!0-9]*___sample_covariates.tsv", emit:sample_covariates) optional true
    
  script:

    """
      echo ${sanitized_columns}
      echo ${phenotype_file}
      echo ${genotype_phenotype_mapping}
      combine_chunked_aggregates.py --n_min_individ ${n_min_individ}
    """
}

process AGGREGATE_UMI_COUNTS {

    
    publishDir  path: "${params.outdir}/aggregated_counts/${sanitized_columns}",mode: "${params.copy_mode}",
              overwrite: "true",
              saveAs: { fn ->
                  def name = fn instanceof Path ? fn.getName() : fn.toString()
                  // Skip chunked files like "1--foo___phenotype_file.tsv"
                  if (name ==~ /^\d+--.*/) return null
                  // otherwise publish
                  return name
              }
    
    
    label 'process_medium'
    tag "$sanitized_columns"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }
    memory { 
            sizeInGB = adata.size() / 1e9 * 0.3 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
  input:
    path(adata) // lists input files per donor
    val(agg_columns)
    val(gt_id_column)
    val(sample_column)
    val(n_cells_min)
    val(n_donors_min)
    val(mode)

  output:
    path("*phenotype_file.tsv", emit:phenotype_file) optional true
    path("*genotype_phenotype_mapping.tsv", emit:genotype_phenotype_mapping) optional true
    path("*___sample_covariates.tsv", emit:sample_covariates) optional true
    tuple val(sanitized_columns), path("*phenotype_file.tsv"),  path("*genotype_phenotype_mapping.tsv"), emit:phenotype_genotype_file optional true
    tuple val(sanitized_columns), path("*phenotype_file.tsv"),  path("*genotype_phenotype_mapping.tsv"), emit:phenotype_genotype_file_with_sample_covs optional true
    path('*.tsv') optional true
    path('*log.txt') optional true
  script:

    adata_name = adata.getName().split('--')[0]
    
    sanitized_columns = adata_name.replaceAll(/[^a-zA-Z0-9]/, '_').replaceAll(/\.h5ad$/, '')
    """
      echo ${sanitized_columns}
      aggregate_sc_data.py --agg_columns '${agg_columns}' --gt_id_column '${gt_id_column}' --sample_column '${sample_column}' --n_cells ${n_cells_min} -n_individ ${n_donors_min} -h5ad ${adata} --method ${params.aggregation_method} --cell_percentage_threshold ${params.cell_percentage_threshold} --covariates '${params.covariates.adata_obs_covariate}' --n_workers  ${task.cpus} 
    """
}