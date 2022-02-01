process PREPROCESS_GENOTYPES{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/eqtl.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        path(file__vcf)
    output:
        path("filtered_vcf.vcf.gz") , emit: filtered_vcf
    script:
        """
            bcftools view --known --max-alleles 2 -m2 -M2 -v snps  -q 0.1:minor ${file__vcf} -O z -o filtered_vcf.vcf.gz
            
        """
    
}