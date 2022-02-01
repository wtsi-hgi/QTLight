process SPLIT_PHENOTYPE_DATA{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run

    tag{condition}
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/eqtl.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        path(annotation_file)
        path(phenotype_file)
        val(condition)
    output:
        tuple val(condition),path("*_phenotype.tsv"), emit: phenotye_file
        
    script:
        
        """
            split_phenotype_for_condition.py --condition '${condition}' --genome_phenotype ${annotation_file} --phenotype ${phenotype_file}
        """
    
}