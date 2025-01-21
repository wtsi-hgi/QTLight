
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
        if (params.genotypes.use_gt_dosage) {
            pgen_or_bed = "--pfile"
        }else{
            pgen_or_bed = "--bfile"
        }
        pgen_or_bed = "--bfile" //ATM only bed is used in LIMIX
        """
            plink2 --freq counts ${pgen_or_bed} ${plink_path}/plink_genotypes --out tmp_gt_plink_freq
            plink2 --make-rel square --read-freq tmp_gt_plink_freq.acount ${pgen_or_bed} ${plink_path}/plink_genotypes
            generate_kinship.py
        """
}
