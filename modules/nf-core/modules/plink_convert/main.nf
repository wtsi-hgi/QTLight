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
            plink2 --make-bed --allow-extra-chr 0 --chr 1-22 XY --output-chr chrM --hwe 1e-6 --indep-pairwise 250 50 0.2 --vcf ${file__vcf} --snps-only --rm-dup exclude-all --out plink_genotypes/plink_genotypes
        """
    
}