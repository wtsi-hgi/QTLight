
process TENSORQTL {
    label 'gpu'
    tag {condition}
    // label 'process_high_memory'
    
  // /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    // container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img"
    container "${params.eqtl_container}"
  } else {
    container "${params.eqtl_docker}"
  }


  input:
    tuple(val(condition),path(aggrnorm_counts_bed),path(genotype_pcs_tsv))
    // each path(genotype_pcs_tsv)
    each path(plink_files_prefix)

  output:
    // path("mapqtl_${oufnprfx}.cis_eqtl.tsv.gz"), emit: qtl_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_dropped.tsv.gz"), emit: dropped_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_qval.tsv.gz"), emit: qval_tsv
    path("Cis_eqtls.tsv"), emit: qtl_bin, optional: true
    path("Cis_eqtls_qval.tsv"), emit: q_qtl_bin, optional: true
    path('nom_output')

  script:
  // If a file with interaction terms is provided, use the interaction script otherwise use the standard script
  tensor_qtl_script = params.TensorQTL.interaction ? "tensorqtl_analyse_interaction.py --inter ${params.TensorQTL.interaction_file}" : "tensorqtl_analyse.py"   
    
    """
      ${tensor_qtl_script} --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed ${aggrnorm_counts_bed} --covariates_file ${genotype_pcs_tsv} -window ${params.windowSize} -nperm ${params.numberOfPermutations}
    """
}


workflow TENSORQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        
    main:
  
      TENSORQTL(
          condition_bed,
          plink_genotype
      )
}
