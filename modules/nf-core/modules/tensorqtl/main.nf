
process TENSORQTL {
    label 'gpu'
    // tag {condition}
    // label 'process_high_memory'
    
  // /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    // container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img"
    container "/software/hgi/containers/eqtl.img"
  } else {
    container "wtsihgi/nf_cellbender_container:3cc9983"
  }


  input:
    tuple(val(condition),path(aggrnorm_counts_bed))
    each path(genotype_pcs_tsv)
    each path(plink_files_prefix)

  output:
    // path("mapqtl_${oufnprfx}.cis_eqtl.tsv.gz"), emit: qtl_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_dropped.tsv.gz"), emit: dropped_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_qval.tsv.gz"), emit: qval_tsv
    path("Cis_eqtls.tsv"), emit: qtl_bin
    path("Cis_eqtls_qval.tsv"), emit: q_qtl_bin

  script:

    
    """
      tensorqtl_analyse.py --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed ${aggrnorm_counts_bed} --covariates_file ${genotype_pcs_tsv} -window ${params.window} -nperm ${params.nperm}
    """
}


workflow TENSORQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        genotype_pcs_tsv
    main:
  
      TENSORQTL(
          condition_bed,
          genotype_pcs_tsv,
          plink_genotype
      )
}