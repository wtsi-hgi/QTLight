process PREPROCESS_SAMPLE_MAPPING{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/eqtl.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        path(genotype_phenotype)
    output:
        path("genotype_phenotype.tsv") , emit: genotype_phenotype
    script:
        """
        
          genotype_phenotype_preprocess.py --genotype_phenotype ${genotype_phenotype}
        """
    
}