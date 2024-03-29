
process KINSHIP_CALCULATION {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"
    
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }


    input:
        path(plink_path)

    output:
    
        path("kinship_matrix.tsv"), emit: kinship_matrix

    script:
        """
            plink2 --freq counts --bfile ${plink_path}/plink_genotypes --out tmp_gt_plink_freq
            plink2 --make-rel square --read-freq tmp_gt_plink_freq.acount --bfile ${plink_path}/plink_genotypes
            generate_kinship.py
        """
}
