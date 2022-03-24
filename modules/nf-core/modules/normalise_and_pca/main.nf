process NORMALISE_and_PCA_PHENOTYPE{
     
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    tag { condition }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'

    publishDir  path: "${params.outdir}/norm_data/${condition}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/eqtl.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        tuple(val(condition),path(phenotype_file))
        path(grouping_file)
        
    output:
        tuple(val(condition),path("normalised_phenotype.tsv"), path("pcs.tsv") , emit: filtered_phenotype)
        tuple(val(condition),path(grouping_file),path("normalised_phenotype.tsv"),path("pcs.tsv"), emit: for_bed)
        val(condition), emit: cond1
        path("*.pdf")
        path(phenotype_file)
    script:
    
        """ 
            normalise_and_pca.R ${phenotype_file} ${grouping_file} ${params.filter_method}
        """
    
}