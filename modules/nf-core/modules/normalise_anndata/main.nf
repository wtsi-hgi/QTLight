

process NORMALISE_ANNDATA {
    publishDir  path: "${outdir}/normalise_anndata",mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium_memory'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
      path(adata)

    output:
      path("normAnnData.h5ad", emit:adata)
    script:
      outdir = params.outdir
      """
        normalise_anndata.py -h5ad ${adata} --method ${params.dMean_norm_method}
      """
}


process REMAP_GENOTPE_ID{
    publishDir  path: "${params.outdir}/normalise_anndata",mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium_memory'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
      path(adata_emmited_file)
      path(mapping_file)
    output:
      path("remap_genotype_phenotype_mapping.tsv", emit:remap_genotype_phenotype_mapping)

    script:
      """
        replace_genotype_ids.py --mappings ${mapping_file} --genotype_phenotype_mapping ${adata_emmited_file}
      """


}