process PREPROCESS_GENOTYPES{
    
    // Calulates bbknn neighbors and saves UMAPS of these
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
            bcftools index ${file__vcf}
            bcftools view --known ${params.bcftools_filters} ${file__vcf} -Oz -o filtered_vcf.vcf.gz
            #bcftools sort filtered_vcf.vcf.gz -Oz -o filtered_vcf2.vcf.gz
            #rm filtered_vcf.vcf.gz
        """
    
}