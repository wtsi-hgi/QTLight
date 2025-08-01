process MERGE_COVARIATES {

    tag { condition }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_low'


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),val(condition2),path(mappings_handeling_repeats),path(normalised_phenotype),path(pc1),path(extra_cove))

    output:
        tuple val(condition2),path(mappings_handeling_repeats),path(normalised_phenotype),path("proc_*pcs.tsv"), emit: for_bed_covs optional true

    script:

        """
            echo "Phenotype: ${condition}"
            echo "Mapping: ${condition2}"
            echo "PCs: ${mappings_handeling_repeats}"
            echo "Sample covariates: ${normalised_phenotype}"
            echo "Sample covariates: ${extra_cove}"

            merge_pc_and_meta.py --meta ${extra_cove} --pc ${pc1} --output proc_${pc1}
        """
}

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
        tuple val(condition),path(mappings_handeling_repeats),path(normalised_phenotype),path("${pc1}pcs.tsv"), emit: for_bed optional true

    script:
    // Here we split the pcs based on the input file. This is removed from normalise part to avoid normalisation multiple times if new pcs asre added.
        """  
            PC_subset.R ${all__pcs} ${pc1}
        """
}