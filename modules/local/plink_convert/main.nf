process PGEN_TO_BED_CONVERT_FOR_QTLS {

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
        mkdir plink_genotypes_tmp

        # Step 1: Convert pgen to bed with relaxed hard-call threshold
        plink2 --pfile "\$plink_dir/\$base_name" \\
            --make-bed \\
            --hard-call-threshold ${params.genotypes.hard_call_threshold} \\
            --out plink_genotypes_tmp/plink_genotypes

        # Step 2: Filter SNPs with missing
        plink2 --bfile plink_genotypes_tmp/plink_genotypes \\
            --geno ${params.genotypes.geno} \\
            --make-bed \\
            --out plink_genotypes_bed/plink_genotypes
        rm -r plink_genotypes_tmp
        """
}

process PGEN_TO_BED_CONVERT_FOR_GRM {

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
        # Detect base name
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
        mkdir -p plink_genotypes_bed

        # Step 1: Convert PGEN to BED with hard-call threshold
        plink2 --pfile "\$plink_dir/\$base_name" \\
               --hard-call-threshold 0.05 \\
               --make-bed \\
               --out plink_genotypes_bed/tmp_hardcall

        # Step 2: Filter missingness (approximate --geno 0.01)
        plink2 --bfile plink_genotypes_bed/tmp_hardcall \\
               --geno 0.01 \\
               --make-bed \\
               --out plink_genotypes_bed/tmp_filtered

        # Step 3: LD pruning
        plink2 --bfile plink_genotypes_bed/tmp_filtered \\
               ${params.covariates.genotype_pc_filters} \\
               --out plink_genotypes_bed/prune

        # Step 4: Extract pruned SNPs
        plink2 --bfile plink_genotypes_bed/tmp_filtered \\
               --extract plink_genotypes_bed/prune.prune.in \\
               --make-bed \\
               --out plink_genotypes_bed/plink_genotypes

        rm plink_genotypes_bed/tmp_hardcall.* plink_genotypes_bed/tmp_filtered.*
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
            ext1 = "mkdir -p temp && plink2 --vcf ${file__vcf} --make-pgen --out temp/temp"
        } else if ("${file__vcf}".contains(".bcf")){
            ext1 = "mkdir -p temp && plink2 --bcf ${file__vcf} --make-pgen --out temp/temp"
        }else{
            ext1 = "ln -s ${file__vcf} temp"
        }

        """
            
            ${ext1}
            
            plink_dir="temp"
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

            plink2 --pfile "\$plink_dir/\$base_name" --export bgen-1.2 ref-first --out genotypes_bgen
            bgenix -g genotypes_bgen.bgen -index
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
            mkdir plink_genotypes_tmp

            # Step 1: Convert pgen to bed with relaxed hard-call threshold
            plink2 ${ext1} ${file__vcf} \\
                --make-bed \\
                --hard-call-threshold ${params.genotypes.hard_call_threshold} --hwe ${params.hwe}  \\
                --out plink_genotypes_tmp/plink_genotypes

            # Step 2: Filter SNPs with missing
            plink2 --bfile plink_genotypes_tmp/plink_genotypes \\
                --geno ${params.genotypes.geno} \\
                --make-bed \\
                --out plink_genotypes_bed/plink_genotypes
            rm -r plink_genotypes_tmp

        """
    
}

process PLINK_CONVERT__GRM{
    
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
            mkdir -p plink_genotypes_bed

            # Step 1: Convert PGEN to BED with hard-call threshold
            plink2 ${ext1} ${file__vcf}  \\
                --hard-call-threshold 0.05 \\
                --make-bed \\
                --out plink_genotypes_bed/tmp_hardcall

            # Step 2: Filter missingness (approximate --geno 0.01)
            plink2 --bfile plink_genotypes_bed/tmp_hardcall \\
                --geno 0.01 \\
                --make-bed \\
                --out plink_genotypes_bed/tmp_filtered

            # Step 3: LD pruning
            plink2 --bfile plink_genotypes_bed/tmp_filtered \\
                ${params.covariates.genotype_pc_filters} \\
                --out plink_genotypes_bed/prune

            # Step 4: Extract pruned SNPs
            plink2 --bfile plink_genotypes_bed/tmp_filtered \\
                --extract plink_genotypes_bed/prune.prune.in \\
                --make-bed \\
                --out plink_genotypes_bed/plink_genotypes

            rm plink_genotypes_bed/tmp_hardcall.* plink_genotypes_bed/tmp_filtered.*

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