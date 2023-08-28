process OPTIMISE_PCS{
     
    // Choose the best eQTL results based on most eGenes found over a number of PCs
    // ------------------------------------------------------------------------
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'


    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/OPTIM_pcs",
                mode: "copy",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(eqtl_dir))
        
    output:
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.pdf"), emit: optimise_nPCs_plot
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.txt"), emit: optimise_nPCs
        path("${outpath}/Cis_eqtls.tsv"), emit: optim_qtl_bin, optional: true
        path("${outpath}/Cis_eqtls_qval.tsv"), emit: optim_q_qtl_bin, optional: true
        path("${outpath}/cis_inter1.cis_qtl_top_assoc.txt.gz "), emit: optim_int_qtl_bin, optional: true
        path(outpath)
    script:
        alpha = 0.05
        alpha_text = alpha.replaceAll("\\.", "pt")
      if (params.TensorQTL.interaction_file?.trim()) {
        inter_name = file(params.TensorQTL.interaction_file).baseName
        outpath = "./interaction_output/${inter_name}"
    } else {
        inter_name = "NA"
        outpath = "./base_output/base"
        """  
            Rscript tensorqtl_optimise.R .. ${alpha} ${inter_name} ${condition} ${outpath}
            var=$(grep 'TRUE' ${outpath}/optimise_nPCs-FDR${alpha_text}.txt | cut -f 1) && cp -r ../"$var"pcs/${outpath}/* ${outpath}
        """
    
}