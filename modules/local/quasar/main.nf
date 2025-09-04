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
        tuple val(condition), 
        path(phenotype_file), 
        path(phenotype_pcs), 
        val(pcs),
        path(annotation_file),
        path(bim), 
        path(bed), 
        path(fam),
        path(sparseGRM),
        path(sparseGRM_sample),
        val(model),
        val (mode)
        

    output:
        tuple val(condition), path('*-variant.txt'), emit: quasar_variant
        tuple val(condition), path('*-region.txt'), emit: quasar_region

    script:
    //TO DO make GRM input more efficient - doesn't need parsing to tsv if not using lmm, nb_glmm or p_glm
    def outname="quasar_${model}"
    def apl = (model in ['nb_glm','nb_glmm']) ? '--use-apl' : ''
    def grm = (model in ['lmm', 'nb_glmm','p_glmm']) ? '--grm grm.tsv' : ''

    """
    zcat ${phenotype_file} | sed s'/gene_id/phenotype_id/' | grep -E '^#chr|^(chr|[0-9])' > phenotype.bed

    transpose_covs.py --infile ${phenotype_pcs} --outfile Covariates.fixed_tmp.tsv
    sed s'/ /_/g' Covariates.fixed_tmp.tsv > Covariates.fixed.tsv

    convert_grm_to_tsv.py --grm ${sparseGRM} --samples ${sparseGRM_sample} --output grm.tsv

    echo ${bim.baseName}

    quasar \
        --plink ${bim.baseName} \
        --bed phenotype.bed \
        --cov Covariates.fixed.tsv \
        ${apl} \
        ${grm}  \
        --mode ${mode} \
        --model ${model} \
        --out ${outname}
    """

}