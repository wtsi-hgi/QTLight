
process GENOTYPE_PC_CALCULATION {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"
    
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        log.info 'change the docker container - this is not the right one'
        container "${params.eqtl_docker}"
    }


    input:
        path(plink_bed)

    output:
    
        path("gtpca_plink.eigenvec"), emit: gtpca_plink

    script:

        """
            plink2 --freq counts --bfile ${plink_bed}/plink_genotypes --out tmp_gt_plink_freq
            plink2 --pca ${params.covariates.nr_genotype_pcs} --read-freq tmp_gt_plink_freq.acount  --bfile ${plink_bed}/plink_genotypes --out gtpca_plink
        """
}
