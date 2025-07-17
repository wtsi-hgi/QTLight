
process PREPERE_EXP_BED {
  label 'process_medium'
  tag "$condition"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
      // container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-singularity_images-eqtl_19_09_2023.img.img"
      
  } else {
      container "${params.eqtl_docker}"
  }


  input:
    tuple(path(phenotype_pcs),val(condition),path(mapping_file),path(expression_file))
    each path(annotation_file)

  output:
    tuple(val("${condition}__${phenotype_pcs}"),path("Expression_Data.bed.gz"), emit: exp_bed)

  script:
    if(params.covariates.extra_covariates_file==''){
      sample_covar =''
    }else{
      sample_covar ="--sample_covariates ${params.covariates.extra_covariates_file}"
    }


    if ("${params.chromosomes_to_test}"!=''){
        chromosomes_as_string = params.chromosomes_to_test.join(',')
        cond2 = " --chr ${chromosomes_as_string}"
    }else{
        cond2 = " "
    }

    """
      echo ${condition}
      prepere_bed.py --annotation_file ${annotation_file} --mapping_file ${mapping_file} --expression_file ${expression_file} --position ${params.position} --gtf_gene_identifier ${params.gtf_gene_identifier} --gtf_type ${params.gtf_type} ${cond2} 
    """
}

process PREPERE_COVARIATES {
  label 'process_medium'
  tag "$condition, $nr_phenotype_pcs"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
  } else {
      container "${params.eqtl_docker}"
  }

  input:
    tuple(path(phenotype_pcs),val(condition),path(mapping_file),path(expression_file))
    each path(genotype_pcs)

  output:
    tuple(val("${condition}__${phenotype_pcs}"),path('Covariates.tsv'), val(nr_phenotype_pcs), emit: exp_bed)

  script:
    nr_phenotype_pcs = phenotype_pcs.getSimpleName()
    if(params.covariates.extra_covariates_file==''){
      sample_covar =''
    }else{
      sample_covar ="--sample_covariates ${params.covariates.extra_covariates_file}"
    }

    """
      echo ${condition}
      prepere_covariates_file.py --genotype_pcs ${genotype_pcs} --phenotype_pcs ${phenotype_pcs} ${sample_covar} --sample_mapping ${mapping_file} --nr_gPCs ${params.covariates.nr_genotype_pcs}
    """
}

process PREP_SAIGE_COVS {
  label 'process_medium'
  tag "$condition"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
  } else {
      container "${params.eqtl_docker}"
  }

  input:
    path(genotype_pcs)
    path(extra_covariates_file)
  output:
    path('G_E_Covariates.tsv')

  script:

    if(params.covariates.extra_covariates_file==''){
      sample_covar =''
    }else{
      sample_covar ="--sample_covariates ${extra_covariates_file}"
    }
    // subset to number of gPCs required and also append the extra covariates that may be provided. 
    """
      prepere_covariates_file_SAIGE.py --genotype_pcs ${genotype_pcs} ${sample_covar} --nr_gPCs ${params.covariates.nr_genotype_pcs}
    """
}