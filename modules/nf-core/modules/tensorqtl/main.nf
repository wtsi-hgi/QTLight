tensor_label = params.utilise_gpu ? 'gpu' : "process_medium"   

process TENSORQTL {  
    label "${tensor_label}"
    tag "$condition, $nr_phenotype_pcs, optim: $skip_nominal"
    // cache false
    
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.eqtl_container}"
  } else {
    container "${params.eqtl_docker}"
  }
  

  input:
    tuple(val(condition),path(aggrnorm_counts_bed),path(covariates_tsv),val(nr_phenotype_pcs))
    each path(plink_files_prefix)
    path(interaction_file)
    val(skip_nominal)

  output:
    tuple val(condition), path("${outpath}"), emit: pc_qtls_path

  script:
  // If a file with interaction terms is provided, use the interaction script otherwise use the standard script   
  if ("${interaction_file}" != 'fake_file.fq') {
    tensor_qtl_script = "tensorqtl_analyse_interaction.py -inter ${interaction_file} --interaction_maf ${params.TensorQTL.interaction_maf}"
    inter_name = file(interaction_file).baseName
    outpath = "${nr_phenotype_pcs}/interaction_output/${inter_name}"
  } else {
    tensor_qtl_script = "tensorqtl_analyse.py -nperm ${params.numberOfPermutations}"
    outpath = "${nr_phenotype_pcs}/base_output/base"
  }

  if (skip_nominal) {
    map_nominal_flag = ""
  } else {
    map_nominal_flag = "--map_nominal"
  }

  if (params.genotypes.use_gt_dosage) {
    dosage = "--dosage"
  }else{
    dosage = ""
  }
    """ 
      bedtools sort -i ${aggrnorm_counts_bed} -header > Expression_Data.sorted.bed
      sed -i 's/^chr//' Expression_Data.sorted.bed
      ${tensor_qtl_script} --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed Expression_Data.sorted.bed --covariates_file ${covariates_tsv} -window ${params.windowSize} ${dosage} --maf ${params.maf} --outdir ${outpath} ${write_nominal_flag}
      cd ${outpath} && ln ../../../${covariates_tsv} ./ && ln ../../../Expression_Data.sorted.bed
    """
}

// PREP_OPTIMISE_PCS process to create symlinks
process PREP_OPTIMISE_PCS {
    label 'process_low'
    tag { condition }
    input:
    tuple val(condition), val(paths)

    output:
    tuple val(condition), path("${condition}_symlink")

    script:
    paths_str = paths.join(" ")
    """
    mkdir ${condition}_symlink
    cd ${condition}_symlink
    for path in ${paths_str}; do
         unique_name=\$(echo \$path | awk -F/ '{print \$(NF-2)"__"\$(NF-1)"__"\$NF}')
        ln -s \$path \$unique_name || echo 'already liked'
    done
    """
}

process OPTIMISE_PCS{
     
    // Choose the best eQTL results based on most eGenes found over a number of PCs
    // ------------------------------------------------------------------------
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    errorStrategy 'ignore'


    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(eqtl_dir))
        path(interaction_file)
        
    output:
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.pdf"), emit: optimise_nPCs_plot
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.txt"), emit: optimise_nPCs
        path("${outpath}/Cis_eqtls.tsv"), emit: optim_qtl_bin, optional: true
        path("${outpath}/Cis_eqtls_qval.tsv"), emit: optim_q_qtl_bin, optional: true
        path("${outpath}/Cis_eqtls_independent.tsv"), emit: optim_independent_qtl_bin, optional: true
        path("${outpath}/cis_inter1.cis_qtl_top_assoc.txt.gz "), emit: optim_int_qtl_bin, optional: true
        tuple val(condition), path("${outpath}/Covariates.tsv"), path("${outpath}/Expression_Data.sorted.bed"), path("${outpath}/Cis_eqtls_qval.tsv"), emit: combined_input, optional: true
        path(outpath)
        

    script:
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
      if ("${interaction_file}" != 'fake_file.fq') {
          inter_name = file(interaction_file).baseName
          outpath_end = "interaction_output__${inter_name}"
      } else {
          inter_name = "NA"
          outpath_end = "base_output__base"
          }
        alpha = "0.05"
        alpha_text = alpha.replaceAll("\\.", "pt")
        outpath = "./OPTIM_pcs/${outpath_end}"
        """  
          mkdir -p ${outpath}
          tensorqtl_optimise_pcs.R ./ ${alpha} ${inter_name} ${condition} ${outpath}
          var=\$(grep TRUE ${outpath}/optimise_nPCs-FDR${alpha_text}.txt | cut -f 1) && cp -r ${condition}_symlink/"\$var"pcs__${outpath_end}/* ${outpath}
          echo \${var} >> ${outpath}/optim_pcs.txt
        """
}

process TRANS_BY_CIS {
    label "${tensor_label}"
    tag "$condition"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
    } else {
      container "${params.eqtl_docker}"
    }

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:

        tuple(
          val(condition),
          path(covariates),
          path(phenotype_file),
          path(cis_eqtls_qval)
        )
        each path(plink_files_prefix) 

    output:
        //path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        //path("trans-by-cis_all.tsv.gz", emit: trans_res_all, optional:true)

    script:
      // Use dosage?
      if (params.genotypes.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }
      if (params.TensorQTL.trans_by_cis_variant_list !='') {
        command_subset_var_list = "grep -P 'variant_id\tcondition_name|${condition}' ${params.TensorQTL.trans_by_cis_variant_list} > var_list.tsv"
        variant_list = "--variant_list var_list.tsv"
      }else{
		command_subset_var_list = ""
        variant_list = ""
      }

      """
	  ${command_subset_var_list}
      tensor_analyse_trans_by_cis.py \
        --covariates_file ${covariates} \
        --phenotype_file ${phenotype_file} \
        --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
        --outdir "./" \
        ${dosage} \
        --maf ${params.maf} \
        --cis_qval_results ${cis_eqtls_qval} \
        --alpha ${params.TensorQTL.alpha} \
        --window ${params.windowSize} \
        ${variant_list}
        
      """

      
}

process TRANS_OF_CIS {
    label "process_high_memory"
    tag "$condition"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
    } else {
      container "${params.eqtl_docker}"
    }

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
          val(condition),
          path(covariates),
          path(phenotype_file),
          path(cis_eqtls_qval)
        )
        each path(plink_files_prefix) 

    output:
        //path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-of-cis_all.tsv", emit: trans_res, optional: true)
        path("trans-of-cis_filt.tsv", emit: trans_res_filt, optional: true)

    script:
      // Use dosage?
      if (params.genotypes.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }

      alpha = "0.05"

      """
      tensor_analyse_trans_of_cis.py \
        --covariates_file ${covariates} \
        --phenotype_file ${phenotype_file} \
        --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
        --outdir "./" \
        --dosage ${dosage} \
        --maf ${params.maf}  \
        --cis_qval_results ${cis_eqtls_qval} \
        --alpha ${alpha} \
        --window ${params.windowSize} \
        --pval_threshold ${params.TensorQTL.trans_by_cis_pval_threshold}
      """
      //cp trans-by-cis_bonf_fdr.tsv ${outpath}
      //"""
}

workflow TENSORQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        
    main:
  
      if(params.TensorQTL.interaction_file !=''){
          int_file = params.TensorQTL.interaction_file
      }else{
          int_file = "$projectDir/assets/fake_file.fq"
      }

      TENSORQTL(
          condition_bed,
          plink_genotype,
          int_file,
          params.TensorQTL.optimise_pcs
      )

      if (params.TensorQTL.optimise_pcs){
          // TENSORQTL.out.pc_qtls_path.view()
          // Make sure all input files are available before running the optimisation
          
          TENSORQTL.out.pc_qtls_path.collect().subscribe { println "TENSORQTL dist: $it" }
          // Fix the format of the output from TENSORQTL
          prep_optim_pc_channel = TENSORQTL.out.pc_qtls_path.groupTuple().map { key, values -> [key, values.flatten()] }
          // Create symlinks to the output files
          PREP_OPTIMISE_PCS(prep_optim_pc_channel)
          // Run the optimisation to get the eQTL output with the most eGenes
          OPTIMISE_PCS(PREP_OPTIMISE_PCS.out,int_file)

          TENSORQTL(
            OPTIMISE_PCS.out.combined_input,
            plink_genotype,
            int_file,
            false
          )
          
          if(params.TensorQTL.trans_by_cis){
            log.info 'Running trans-by-cis analysis on optimum nPCs'
            TRANS_BY_CIS(
              OPTIMISE_PCS.out.combined_input,
              plink_genotype
            )
          }
          
          if(params.TensorQTL.trans_of_cis){
            log.info 'Running trans-of-cis analysis on optimum nPCs'
            TRANS_OF_CIS(
              OPTIMISE_PCS.out.combined_input,
              plink_genotype
            )
          }
  }
}
