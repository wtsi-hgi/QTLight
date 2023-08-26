process OPTIMISE_PCS{
     
    // Choose the best eQTL results based on most eGenes found over a number of PCs
    // ------------------------------------------------------------------------
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'


    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/optim_pcs",
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
        path("optimise_nPCs.pdf"), emit: optimise_nPCs_plot
        path("optimise_nPCs.txt"), emit: optimise_nPCs
        path("Cis_eqtls.tsv"), emit: qtl_bin, optional: true
        path("Cis_eqtls_qval.tsv"), emit: q_qtl_bin, optional: true
        path("cis_inter1.cis_qtl_top_assoc.txt.gz "), emit: int_qtl_bin, optional: true
        path(outpath)
    script:
        """  
            tensorqtl_optimise.R 0.05
            var=$(grep 'TRUE' ./*-optimise_nPCs.txt | cut -f 1); cp -r ../"$var"pcs/* .
        """
    
}