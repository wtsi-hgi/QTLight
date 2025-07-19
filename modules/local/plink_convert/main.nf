process PGEN_TO_BED_CONVERT {

    // Converts PGEN to PLINK BED format (for SAIGE or other tools needing .bed/.bim/.fam)

    scratch false
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        publishDir path: "${params.outdir}/genotypes",
                   mode: "${params.copy_mode}",
                   overwrite: "true"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(pgen_prefix_dir)

    output:
        path("plink_genotypes_bed"), emit: plink_path_bed
        tuple path("plink_genotypes_bed/plink_genotypes.bim"),
              path("plink_genotypes_bed/plink_genotypes.bed"),
              path("plink_genotypes_bed/plink_genotypes.fam"), emit: bim_bed_fam

    script:
        """
        # Detect base name inside pgen_prefix_dir
        plink_dir="${pgen_prefix_dir}"
        base_name=""

        if ls "\$plink_dir"/*.psam 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.psam | head -n 1) .psam)
        elif ls "\$plink_dir"/*.pvar 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.pvar | head -n 1) .pvar)
        else
            echo "No .psam or .pvar file found in \$plink_dir"
            exit 1
        fi

        echo "Detected base name: \$base_name"

        mkdir plink_genotypes_bed

        plink2 --pfile "\$plink_dir/\$base_name" \\
               --make-bed \\
               --out plink_genotypes_bed/plink_genotypes
        """
}



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
        path("plink_genotypes_pgen"), emit: plink_path

    script:

        if ("${file__vcf}".contains(".vcf")) {
            ext1 = "--vcf ${file__vcf} 'dosage=DS' "
        }  else if ("${file__vcf}"=='plink_genotypes_bed') {
            ext1 = "--bfile ${file__vcf}/plink_genotypes"
        }else {
            ext1 = "--bcf ${file__vcf} 'dosage=DS'"
        }

        """
            mkdir plink_genotypes_pgen
            plink2 ${ext1} --max-alleles 2 --make-pgen ${params.plink2_filters} --hwe ${params.hwe} --out plink_genotypes_pgen/plink_genotypes
        """
    
    
}