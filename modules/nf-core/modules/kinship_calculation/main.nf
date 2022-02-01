
process KINSHIP_CALCULATION {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"
    
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/eqtl.img"
    } else {
        log.info 'change the docker container - this is not the right one'
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }


    input:
        path(plink_path)

    output:
    
        path("kinship_matrix.tsv"), emit: kinship_matrix

    script:
        """
            plink2 --freq counts --bfile ${plink_path}/plink_genotypes --out tmp_gt_plink_freq
            plink2 --make-rel square --read-freq tmp_gt_plink_freq.acount --maf 0.05 --bfile ${plink_path}/plink_genotypes
            generate_kinship.py
        """
}
