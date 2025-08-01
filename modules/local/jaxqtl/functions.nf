jax_label = params.JAXQTL.use_gpu ? 'gpu' : "process_high_memory"   

process JAXQTL {  
    tag "$condition, $nr_phenotype_pcs"
    label "${jax_label}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.jax_container}"
    } else {
        container "${params.jax_docker}"
    }

    publishDir path: { "${params.outdir}/JAX_eQTLS/${group0}__${group1}__${group2}/OPTIM_PCs/${mode1}" },
            saveAs: {filename ->
                if (filename.contains("cis_score.tsv.gz")) {
                    null
                } else {
                    filename
                }
            },
            mode: "${params.copy_mode}",
            overwrite: true

    input:
        tuple(
            val(condition),
            path(aggrnorm_counts_bed),
            path(covariates_tsv),
            val(nr_phenotype_pcs),
            path(genelist),
            path(plink_files_prefix)
        )
        val(mode1)

    output:
        tuple val("${condition}__${nr_phenotype_pcs}"), path('*cis_score.tsv.gz'), emit: qtl_data optional true
        tuple val("${condition}__${nr_phenotype_pcs}"), path('*cis_qtl_pairs.*.parquet'), emit: nominal_qtl_data optional true
        

    script:

        group0 = "${condition}".split('__')[0]
        group1 =  "${condition}".split('__')[1]
        group2 =  "${condition}".split('__')[2]

        if (params.JAXQTL.use_gpu){
            gpu_cpu='gpu'
        }else{
            gpu_cpu='cpu'
        }

        """

            export DISABLE_PANDERA_IMPORT_WARNING=True
            plink_dir="${plink_files_prefix}"
            base_name=""
            outname=\$(basename ${genelist})
            if ls "\$plink_dir"/*.psam 1> /dev/null 2>&1; then
                base_name=\$(basename \$(ls "\$plink_dir"/*.psam | head -n 1) .psam)
                pgen_or_bed="--pfile"
            elif ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
                base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
                pgen_or_bed="--bfile"
            else
                echo "No .psam or .bed file found in \$plink_dir"
                exit 1
            fi
            transpose_covs.py --infile ${covariates_tsv} --outfile Covariates.fixed.tsv
            jaxqtl \
            --geno "\$plink_dir/\$base_name"  \
            --covar Covariates.fixed.tsv \
            --pheno ${aggrnorm_counts_bed} --genelist ${genelist}  \
            --model NB \
            --mode ${mode1} \
            --window ${params.windowSize} \
            --test-method score \
            --nperm ${params.numberOfPermutations} \
            --addpc 0 \
            --standardize \
            -p ${gpu_cpu} \
            --out ${mode1}__\$outname


        """
}