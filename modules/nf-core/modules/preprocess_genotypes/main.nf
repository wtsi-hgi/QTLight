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

    publishDir  path: "${params.outdir}/genotypes",
                mode: "${params.copy_mode}",
                overwrite: "true"       
                
    input:
        path(file__vcf)
    output:
        path("filtered_vcf.vcf.gz") , emit: filtered_vcf
    script:
        """
            bcftools index ${file__vcf} || echo 'exists'
            bcftools view ${params.bcftools_filters} ${file__vcf} -Oz -o filtered_vcf.vcf.gz
            #bcftools sort filtered_vcf.vcf.gz -Oz -o filtered_vcf2.vcf.gz
            #rm filtered_vcf.vcf.gz
        """
    
}