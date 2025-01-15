tensor_label = params.utilise_gpu ? 'gpu' : "process_medium"   

include { TENSORQTL as TENSORQTL } from './processes.nf'
include { TENSORQTL as TENSORQTL_OPTIM } from './processes.nf'

// PREP_OPTIMISE_PCS process to create symlinks
process PREP_OPTIMISE_PCS {
    label 'process_low'
    tag "$condition, $interaction"
    input:
    tuple val(condition), val(interaction), val(paths)
    output:
    tuple val(condition), val(interaction), path("${condition}_${interaction}_symlink")

    script:
    paths_str = paths.join(" ")
    """
    mkdir ${condition}_${interaction}_symlink
    cd ${condition}_${interaction}_symlink
    for path in ${paths_str}; do
         unique_name=\$(echo \$path | awk -F/ '{print \$(NF-2)"__"\$(NF-1)"__"\$NF}')
        ln -s \$path \$unique_name || echo 'already linked'
    done
    """
}

process OPTIMISE_PCS{
     
    // Choose the best eQTL results based on most eGenes found over a number of PCs
    // ------------------------------------------------------------------------
    tag "$condition, $interaction"
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
        tuple(val(condition),val(interaction),path(eqtl_dir))
    output:
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.pdf"), emit: optimise_nPCs_plot
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.txt"), emit: optimise_nPCs
        tuple val(condition), path("${tensor_input_path}/Expression_Data.sorted.bed"), path("${tensor_input_path}/Covariates.tsv"), val('OPTIM_pcs'), path("${tensor_input_path}/${interaction_file}"), emit: cis_input, optional: true
        tuple val(condition), path("${tensor_input_path}/Expression_Data.sorted.bed"), path("${tensor_input_path}/Covariates.tsv"), path("${tensor_input_path}/Cis_eqtls_qval.tsv"), emit: trans_input, optional: true
        path(outpath)
        

    script:
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
      if ("${interaction}" != 'base') {
          interaction_file = "${interaction}.tsv"
          outpath_end = "interaction_output__${interaction}"
      } else {
          interaction_file = "fake_file.fq"
          outpath_end = "base_output__base"
          }
        alpha = "${params.TensorQTL.alpha}"
        alpha_text = alpha.replaceAll("\\.", "pt")
        outpath = "./OPTIM_pcs/${outpath_end}"
        tensor_input_path = "./OPTIM_input/${outpath_end}"
        """  
          mkdir -p ${outpath}
          mkdir -p ${tensor_input_path}
          tensorqtl_optimise_pcs.R ./ ${alpha} ${interaction} ${condition} ${outpath}
          var=\$(grep TRUE ${outpath}/optimise_nPCs-FDR${alpha_text}.txt | cut -f 1) && cp -r ${condition}_${interaction}_symlink/"\$var"pcs__${outpath_end}/* ${tensor_input_path}
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

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/${pcs}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:

        tuple(
          val(condition),
          path(phenotype_file),
          path(covariates),
          path(cis_eqtls_qval)
        )
        each path(plink_files_prefix) 

    output:
        //path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-by-cis_all.tsv.gz", emit: trans_res_all, optional:true)

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

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/${pcs}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
          val(condition),
          path(phenotype_file),
          path(covariates),
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

      alpha = "${params.TensorQTL.alpha}"

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

process SPLIT_INTERACTIONS {
    label 'process_small'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
    } else {
      container "${params.eqtl_docker}"
    }

    input:
        path multi_interaction_file
    output:
        path("interaction_files/*.tsv", emit: interactions_files)
    script:
        """
        mkdir -p interaction_files
        num_columns=\$(head -1 ${multi_interaction_file} | awk -F '\\t' '{print NF}')
        if [ "\$num_columns" -eq 2 ]; then
            # Two columns, copy the file as is, retaining its original name
            cp ${multi_interaction_file} interaction_files/
        elif [ "\$num_columns" -ge 3 ]; then
            # Three or more columns, run the Python script to split
            split_interactions.py -i ${multi_interaction_file} -o interaction_files
        else
            echo "Error: Input file must have at least two columns."
            exit 1
        fi
        """
}


process RUN_GSEA {
    // Run fGSEA for each interaction result
    // ------------------------------------------------------------------------
    
    label "process_medium"
    tag "$condition, $interaction"
    // cache false

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/OPTIM_pcs/interaction_output/${interaction}/fgsea",
                overwrite: "true"

    input:
        tuple(
          val(condition), 
          val(interaction), 
          path(sumstats_dir)
        )

    output:
        tuple(
            val(condition),
            val(interaction),
            path("${outfile}-gsea_results.tsv.gz"),
            val(outdir),
            emit: gsea_results
        )

    script:
        outdir = "./fgsea"
        sumstats_path = "${sumstats_dir}/cis_inter1.cis_qtl_top_assoc.txt.gz"
        outfile = "${condition}__${interaction}"
        """
        fgsea_ieQTLs.R \
            --iegenes ${sumstats_path} \
            --ranking_var b_gi \
            --eps 0 \
            --unsigned_ranking \
            --gsets_gene_matrix ${projectDir}/assets/data/gene_set_genes.tsv.gz \
            --gsets_info_file ${projectDir}/assets/data/gene_set_info.tsv.gz \
            --database c2.cp.reactome \
            --n_cores ${task.cpus} \
            --output_file ${outfile} \
            --verbose 
        """
}


workflow TENSORQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        
    main:
      if (params.TensorQTL.interaction_file != '') {
          SPLIT_INTERACTIONS(params.TensorQTL.interaction_file)
          interaction_files = SPLIT_INTERACTIONS.out.interactions_files.flatten()

          condition_bed = condition_bed.combine(interaction_files)
          
      } else {
          interaction_files = Channel.of("$projectDir/assets/fake_file.fq")
          condition_bed = condition_bed.map { cond_bed_item ->
              cond_bed_item + ["$projectDir/assets/fake_file.fq"]
          }
      }


      // if (params.TensorQTL.aggregation_subentry != '') {
      //     log.info("------- Analysing ${params.SAIGE.aggregation_subentry} celltypes ------- ")
      //     // Split the aggregation_subentry parameter into a list of patterns
      //     valid_files = phenotype_file.filter { file ->
      //         params.SAIGE.aggregation_subentry.split(',').any { pattern -> "${file}".contains("__${pattern}__") }
      //     }
      // } else {
      //     log.info('------- Analysing all celltypes ------- ')
      //     valid_files = phenotype_file
      // }

      condition_bed.view { "condition_bed item: $it" }
      // dMean__CD4_Naive_all, /lustre/scratch127/humgen/teams/hgi/mo11/tmp_projects127/sle_project/8.eqtl_analysis/work/6b/2a35dfd61dcea16e2a43d45eaffd0d/Expression_Data.bed.gz, /lustre/scratch127/humgen/teams/hgi/mo11/tmp_projects127/sle_project/8.eqtl_analysis/work/6b/2a35dfd61dcea16e2a43d45eaffd0d/Covariates.tsv, 20pcs, /software/hgi/pipelines/QTLight/QTLight_v1.41/assets/fake_file.fq
      // plink_genotype.subscribe { println "TENSORQTL dist: $it" }
      TENSORQTL(
          condition_bed,
          plink_genotype,
          params.TensorQTL.optimise_pcs,
          true
      )

      if (params.TensorQTL.optimise_pcs){
          // TENSORQTL.out.pc_qtls_path.view()
          // Make sure all input files are available before running the optimisation
          
          
          // Fix the format of the output from TENSORQTL
          prep_optim_pc_channel = TENSORQTL.out.pc_qtls_path.groupTuple(by: [0,1]).map { key1, key2, values -> [key1, key2, values.flatten()]}
          // Create symlinks to the output files
          PREP_OPTIMISE_PCS(prep_optim_pc_channel)
          // Run the optimisation to get the eQTL output with the most eGenes
          OPTIMISE_PCS(PREP_OPTIMISE_PCS.out)

          TENSORQTL_OPTIM(
            OPTIMISE_PCS.out.cis_input,
            plink_genotype,
            false,
            false
          )

        if (params.TensorQTL.interaction_file != '' && params.TensorQTL.run_gsea)
         {
          RUN_GSEA(TENSORQTL_OPTIM.out.pc_qtls_path)
        }

          if(params.TensorQTL.trans_by_cis){
            log.info 'Running trans-by-cis analysis on optimum nPCs'
            TRANS_BY_CIS(
              OPTIMISE_PCS.out.trans_input,
              plink_genotype
            )
          }
          
          if(params.TensorQTL.trans_of_cis){
            log.info 'Running trans-of-cis analysis on optimum nPCs'
            TRANS_OF_CIS(
              OPTIMISE_PCS.out.trans_input,
              plink_genotype
            )
          }
  }
}