tensor_label = params.TensorQTL.utilise_gpu ? 'gpu' : "process_low"   

process TENSORQTL {  
    label "${tensor_label}"
    tag "$condition, $nr_phenotype_pcs"
    
  // /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/${nr_phenotype_pcs}",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    // container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_tensorqtl_1.2.img"
    container "${params.eqtl_container}"
  } else {
    container "${params.eqtl_docker}"
  }


  input:
    tuple(val(condition),path(aggrnorm_counts_bed),path(covariates_tsv),val(nr_phenotype_pcs))
    each path(plink_files_prefix)

  output:
    // path("mapqtl_${oufnprfx}.cis_eqtl.tsv.gz"), emit: qtl_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_dropped.tsv.gz"), emit: dropped_tsv
    // path("mapqtl_${oufnprfx}.cis_eqtl_qval.tsv.gz"), emit: qval_tsv
    path("Cis_eqtls.tsv"), emit: qtl_bin, optional: true
    path("Cis_eqtls_qval.tsv"), emit: q_qtl_bin, optional: true
    path(outpath)

  script:
  // If a file with interaction terms is provided, use the interaction script otherwise use the standard script   
  if (params.TensorQTL.interaction_file?.trim()) {
    tensor_qtl_script = "tensorqtl_analyse_interaction.py --inter ${params.TensorQTL.interaction_file}"
    outpath = 'inter_output'
  } else {
    tensor_qtl_script = "tensorqtl_analyse.py -nperm ${params.numberOfPermutations}"
    outpath = 'nom_output'
  }
  if (params.TensorQTL.use_gt_dosage) {
    dosage = "--dosage"
  }else{
    dosage = ""
  }
    """
      /software/team152/bedtools sort -i ${aggrnorm_counts_bed} -header > Expression_Data.sorted.bed
      ${tensor_qtl_script} --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed Expression_Data.sorted.bed --covariates_file ${covariates_tsv} -window ${params.windowSize} ${dosage}
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
