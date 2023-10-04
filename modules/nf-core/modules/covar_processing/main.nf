process SUBSET_PCS{

     
    // Normalise expression data and perform PCA
    // ------------------------------------------------------------------------
    tag { condition }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'

    publishDir  path: "${params.outdir}/norm_data/${condition}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(mappings_handeling_repeats),path(normalised_phenotype),path(all__pcs),val(pc1))

    output:
        tuple(val(condition),path(mappings_handeling_repeats),path(normalised_phenotype),path("${pc1}pcs.tsv"), emit: for_bed)

    script:
    // Here we split the pcs based on the input file. This is removed from normalise part to avoid normalisation multiple times if new pcs asre added.
        """  
            PC_subset.R ${all__pcs} ${pc1}
        """
    


}