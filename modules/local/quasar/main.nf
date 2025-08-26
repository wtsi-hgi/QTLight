process QUASAR{
    tag "$condition, $nr_phenotype_pcs"
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"    
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        each path(genome_annotation)
        tuple val(condition),path(phenotype_pcs),path(plink_files_prefix),path(phenotype_file)
        val(model)
        val(mode)

    output:
        tuple val(condition), path('*-variant.txt'), emit: quasar_variant
        tuple val(condition), path('*-region.txt'), emit: quasar_region

    script:
    //TO DO for models lmm, p_glmm, nb_glmm need to add GRM input
    """
    outname="quasar_${model}_${mode}"

    quasar \
        --plink ${plink_files_prefix} \
        --bed ${phenotype_file} \
        --cov ${phenotype_pcs} \
        --mode ${mode} \
        --model ${model} \
        --use-apl \
        --out ${outname}
    """

}