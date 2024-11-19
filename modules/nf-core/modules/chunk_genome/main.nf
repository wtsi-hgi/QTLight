process CHUNK_GENOME{
    tag "$condition, $nr_phenotype_pcs"
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        each path(genome_annotation)
        tuple path(phenotype_pcs),val(condition),path(mapping_file),path(phenotype_file)

    output:
        // path("Chunging_file.tsv"), emit: chunking_file_path
        // path('annotation_file_processed.tsv'), emit: annotation_file_processed
        tuple val(condition),path("Chunging_file.tsv"),path('annotation_file_processed.tsv'),path(phenotype_file), path(phenotype_pcs) , emit: filtered_chunking_file
        path('limix_chunking.tsv'), emit: limix_condition_chunking
    script:
        nr_phenotype_pcs = phenotype_pcs.getSimpleName()
        """
            generate_chunking_file.py --genome_annotation ${genome_annotation} --chunk_size ${params.chunkSize} --phenotype_file ${phenotype_file} --covar_file ${phenotype_pcs} --condition ${condition}  --genotype_phenotype_file ${mapping_file}
        """
    
}