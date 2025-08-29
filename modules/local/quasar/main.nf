process QUASAR{
    tag "$condition, $pcs"
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"    
    } else {
        container "${params.eqtl_docker}"
    }

    publishDir  path: "${params.outdir}/QUASAR/$pcs/",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple (
            val(condition), 
            path(phenotype_file), 
            path(phenotype_pcs),
            val(pcs)
        )
        path(annotation_file)
        tuple path(bim), path(bed), path(fam)
        val(mode)
        val(model)

    output:
        tuple val(condition), path('*-variant.txt'), emit: quasar_variant
        tuple val(condition), path('*-region.txt'), emit: quasar_region

    script:
    //TO DO GRM input
    def outname="quasar_${model}"
    def api = (params.model in ['nb_glm', 'nb_glmm']) ? '--use-api' : ''

    """
    zcat ${phenotype_file} | sed s'/gene_id/phenotype_id/' > phenotype.bed

    transpose_covs.py --infile ${phenotype_pcs} --outfile Covariates.fixed_tmp.tsv
    sed s'/ /_/g' Covariates.fixed_tmp.tsv > Covariates.fixed.tsv

    echo ${bim.baseName}

    quasar \
        --plink ${bim.baseName} \
        --bed phenotype.bed \
        --cov Covariates.fixed.tsv \
        ${api} \
        --mode ${mode} \
        --model ${model} \
        --out ${outname}
    """

}