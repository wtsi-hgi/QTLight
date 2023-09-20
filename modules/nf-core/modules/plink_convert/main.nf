process PLINK_CONVERT{
    
    // Converts VCF to PLINK format, makes bed/bim/fam if use_gt_dosage param is false
    // otherwise makes pgen/psam/pvar with dosages

    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(file__vcf)
    output:
        path("plink_genotypes"), emit: plink_path
    script:
    if(params.TensorQTL.use_gt_dosage==true && params.TensorQTL.run==true){
      pgen_or_bed = "'dosage=DS' --make-pgen"
    }else{
      pgen_or_bed = "--make-bed"
    }
        """
            mkdir plink_genotypes
            plink2 --vcf ${file__vcf} ${pgen_or_bed} ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes/plink_genotypes
        """
    
}