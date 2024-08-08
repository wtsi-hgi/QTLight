process PLINK_CONVERT{
    
    // Converts VCF to PLINK format, makes bed/bim/fam if use_gt_dosage param is false
    // otherwise makes pgen/psam/pvar with dosages

    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"

    publishDir  path: "${params.outdir}/genotypes",
                mode: "${params.copy_mode}",
                overwrite: "true"        
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(file__vcf)
    output:
        path("plink_genotypes_bgen"), emit: plink_path
        tuple path("plink_genotypes_bed/plink_genotypes.bim"),path("plink_genotypes_bed/plink_genotypes.bed"),path("plink_genotypes_bed/plink_genotypes.fam"), emit: bim_bed_fam

    script:
        if ("${file__vcf}".contains(".vcf")) {
            ext1 = "--vcf"
        } else {
            ext1 = "--bcf"
        }

        if(params.TensorQTL.use_gt_dosage==true && params.TensorQTL.run==true){
        pgen = "plink2 ${ext1} ${file__vcf} 'dosage=DS' --max-alleles 2 --make-pgen ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes_bgen/plink_genotypes"
        }else{
        pgen = ""
        }

        """
            mkdir plink_genotypes_bgen
            mkdir plink_genotypes_bed
            ${pgen}
            plink2 ${ext1} ${file__vcf} --make-bed ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes_bed/plink_genotypes
            # Sort the .bim file by chromosome and position
            sort -k1,1n -k4,4n plink_genotypes_bed/plink_genotypes.bim > plink_genotypes_bed/plink_genotypes_sorted.bim

            # Replace the original .bim file with the sorted one
            mv plink_genotypes_bed/plink_genotypes_sorted.bim plink_genotypes_bed/plink_genotypes.bim

        """
    
}