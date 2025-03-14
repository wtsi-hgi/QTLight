process NORMALISE_and_PCA_PHENOTYPE{
     
    // Normalise expression data and perform PCA
    // ------------------------------------------------------------------------
    tag { condition }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'

    publishDir  path: "${params.outdir}/norm_data/${condition}_${prefix}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(phenotype_file),path(grouping_file))
        

    output:
        tuple(val(outname),path("normalised_phenotype.tsv"), path("all__pcs.tsv") , emit: filtered_phenotype)
        tuple(val(outname),path('mappings_handeling_repeats.tsv'),path("normalised_phenotype.tsv"),path("all__pcs.tsv"), emit: for_bed)
        val(outname), emit: cond1
        path("*.pdf")
        path(phenotype_file)
        path('mappings_handeling_repeats.tsv'), emit: gen_phen_mapping
    script:
        matcher = (phenotype_file =~ /^([^_]+)___/)
        prefix = matcher ? matcher[0][1] : 'all'
        outname = "${condition}_${prefix}"
        """  
            echo ${prefix}
            echo ${outname}
            normalise_and_pca.R ${phenotype_file} ${grouping_file} ${params.filter_method} ${params.method} ${params.inverse_normal_transform} ${params.norm_method} ${params.percent_of_population_expressed} ${params.use_sample_pca}
        """
    
}