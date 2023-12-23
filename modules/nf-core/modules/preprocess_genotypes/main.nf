process PREPROCESS_GENOTYPES{
    
    // Subsets genotypes
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(file__vcf)
    output:
        path("filtered_vcf.vcf.gz") , emit: filtered_vcf
    script:
        """
            bcftools view --known ${params.bcftools_filters} --min-af ${params.maf}:minor ${file__vcf} -Oz -o filtered_vcf.vcf.gz
            
        """
    
}