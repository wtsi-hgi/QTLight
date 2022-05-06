process CHUNK_GENOME{
    
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
        each path(genome_annotation)
        tuple(val(condition),path(phenotype_file), path(pcs))
    output:
        // path("Chunging_file.tsv"), emit: chunking_file_path
        // path('annotation_file_processed.tsv'), emit: annotation_file_processed
        tuple val(condition),path("Chunging_file.tsv"),path('annotation_file_processed.tsv'),path(phenotype_file), path(pcs) , emit: filtered_chunking_file
        path('limix_chunking.tsv'), emit: limix_condition_chunking
    script:
        """
            generate_chunking_file.py --genome_annotation ${genome_annotation} --chunk_size ${params.chunkSize} --phenotype_file ${phenotype_file} --covar_file ${pcs} --condition ${condition}
        """
    
}