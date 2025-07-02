
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
        path(plink_dir)

    output:
    
        path("gtpca_plink.eigenvec"), emit: gtpca_plink

    script:
        if (params.genotypes.use_gt_dosage) {
            pgen_or_bed = "--pfile"
        }else{
            pgen_or_bed = "--bfile"
        }
        """
            plink_dir="${plink_dir}"
            base_name=""

            if ls "\$plink_dir"/*.psam 1> /dev/null 2>&1; then
                base_name=\$(basename \$(ls "\$plink_dir"/*.psam | head -n 1) .psam)
                pgen_or_bed="--pfile"
            elif ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
                base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
                pgen_or_bed="--bfile"
            else
                echo " No .psam or .bed file found in \$plink_dir"
                exit 1
            fi

            echo "Detected base name: \$base_name"
            echo "Using mode: \$pgen_or_bed"

            plink2 ${pgen_or_bed} "\$plink_dir/\$base_name" ${params.covariates.genotype_pc_filters} --out tmp_gt_plink_freq
            plink2 ${pgen_or_bed} "\$plink_dir/\$base_name" --extract tmp_gt_plink_freq.prune.in --freq --out tmp_gt_plink_freq
            plink2 --pca 10 --read-freq tmp_gt_plink_freq.afreq --extract tmp_gt_plink_freq.prune.in ${pgen_or_bed} "\$plink_dir/\$base_name" --out gtpca_plink
        """
}
