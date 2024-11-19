
process SUBSET_GENOTYPE {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    // publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }


    input:
        path(donor_vcf)
        val(file__reduced_dims)


    output:
    
        path("${samplename}.subset.vcf.gz"), emit: samplename_subsetvcf

    script:
        file__reduced_dims2 = file__reduced_dims.unique().join(",")
        samplename='subset'
        // sample_subset_file = donor_vcf.getSimpleName()
        // sample_names = file__reduced_dims.toString()
        """ 
            bcftools view ${donor_vcf} -s ${file__reduced_dims2} --force-samples -Oz -o ${samplename}.subset.vcf.gz --threads ${task.cpus}
        """
}
