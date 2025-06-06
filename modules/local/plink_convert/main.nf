process BGEN_CONVERT{
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
        path("genotypes_bgen.bgen"), emit: plink_path

    script:
        if ("${file__vcf}".contains(".vcf")) {
            ext1 = "--vcf"
        } else {
            ext1 = "--bcf"
        }


        """
            plink2 ${ext1} ${file__vcf} --make-pgen --out temp
            plink2 --pfile temp --export bgen-1.2 --out genotypes_bgen
        """    
}

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
        path("plink_genotypes_bed"), emit: plink_path
        tuple path("plink_genotypes_bed/plink_genotypes.bim"),path("plink_genotypes_bed/plink_genotypes.bed"),path("plink_genotypes_bed/plink_genotypes.fam"), emit: bim_bed_fam

    script:
        if ("${file__vcf}".contains(".vcf")) {
            ext1 = "--vcf"
        } else {
            ext1 = "--bcf"
        }


        """
            mkdir plink_genotypes_bed
            plink2 ${ext1} ${file__vcf} --make-bed ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes_bed/plink_genotypes

        """
    
}



process PGEN_CONVERT{
    
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

    script:

        if ("${file__vcf}".contains(".vcf")) {
            ext1 = "--vcf ${file__vcf} 'dosage=DS' "
        }  else if ("${file__vcf}"=='plink_genotypes_bed') {
            ext1 = "--bfile ${file__vcf}/plink_genotypes"
        }else {
            ext1 = "--bcf ${file__vcf} 'dosage=DS'"
        }

        """
            mkdir plink_genotypes_bgen
            plink2 ${ext1} --max-alleles 2 --make-pgen ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes_bgen/plink_genotypes
        """
    
    
}