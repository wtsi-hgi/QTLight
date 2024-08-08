

process NORMALISE_ANNDATA {
    publishDir  path: "${params.outdir}/normalise_anndata/${sanitized_columns}",mode: "${params.copy_mode}",
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
      path("nAD_*.h5ad", emit:adata)
    script:
      sanitized_columns = adata.getName().replaceAll(/[^a-zA-Z0-9]/, '_').replaceAll(/\.h5ad$/, '')
      
      if ("${params.dMean_norm_method}"=='NONE'){
        comand="ln -s ${adata} nAD_${adata}"
      }else{
        comand="normalise_anndata.py -h5ad ${adata} --method ${params.dMean_norm_method}"
      }

      """
        ${comand}
      """
}


process REMAP_GENOTPE_ID{
    publishDir  path: "${params.outdir}/normalise_anndata",mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
      tuple val(sanitized_columns), path(phenotype_file),  path(adata_emmited_file)
      path(mapping_file)
    output:
      tuple val(sanitized_columns), path(phenotype_file),  path("remap_genotype_phenotype_mapping.tsv"), emit:remap_genotype_phenotype_mapping
      path("remap_genotype_phenotype_mapping.tsv"),emit: genotype_phenotype_mapping
    script:
      """
        replace_genotype_ids.py --mappings ${mapping_file} --genotype_phenotype_mapping ${adata_emmited_file}
      """


}