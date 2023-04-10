

process NORAMLISE_ANNDATA {
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
    path("normalised_anndata.h5ad", emit:adata)
  script:
  outdir = params.outdir
  """
    scTransform.py -h5ad ${adata}
  """
}