process PLINK_CONVERT{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    // --max-alleles 2
    // --snps-only 
    // --threads ${task.cpus}
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
        path("plink_genotypes"), emit: plink_path
    script:
        """
            mkdir plink_genotypes
            plink2 --make-bed ${params.plink2_filters} --hwe ${params.hwe} --vcf ${file__vcf} --out plink_genotypes/plink_genotypes
        """
    
}