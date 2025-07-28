include {JAXQTL;JAXQTL as JAXQTL_NOMINAL} from './functions.nf'

process AGGREGATE_QTL_RESULTS{
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    

    publishDir  path: "${params.outdir}/JAX_eQTLS/${group0}__${group1}__${group2}/${group4}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(condition), path(full_output_list)
        
    output:
        tuple val("${group0}__${group1}__${group2}"), path("*merged_cis_scores.tsv.gz"), emit: jax_qtl_path

    script:
        all_qtl_results = "${condition}".replaceFirst(/.*\.tsv$/, '')

        group0 = "${condition}".split('__')[0]
        group4 = "${condition}".split('__')[4]
        group1 =  "${condition}".split('__')[1]
        group2 = "${condition}".split('__')[2].replaceFirst(/\.tsv/, '')

        """
            echo "${condition} ${group1} ${group2} ${group0}"
            merge_chunks.py -o ${group0}__${group1}__${group2}__${group4}__merged_cis_scores.tsv.gz
        """
}


process OPTIM_PCS{
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    

    publishDir  path: "${params.outdir}/JAX_eQTLS/${condition}/OPTIM_PCs",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(condition), path(full_combined_lists)
        
    output:
        tuple val(condition), path("results/optimise_nPCs-FDR*_optimal_PC.txt"), emit: optimal_pc_file
        path('results/*')

    script:

        """
            echo "${full_combined_lists}"
            plot_optimal_pcs.R ./ 0.05 "${condition}" ./results
        """
}







process CHUNK_BED_FILE{
    tag "$condition, $nr_phenotype_pcs"
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(
            val(condition),
            path(aggrnorm_counts_bed),
            path(covariates_tsv),
            val(nr_phenotype_pcs)
        )
        val(chunk_size)
        each path(plink_genotype) 

    output:
        tuple val(condition),path(aggrnorm_counts_bed), path(covariates_tsv),val(nr_phenotype_pcs),path("gene_chunks*.tsv"),path(plink_genotype), emit: chunked_bed_channel optional true

    when:
        "${condition}".contains("dSum") //Jax is supposed to work only of dSum

    script:
        """
            echo ${aggrnorm_counts_bed}
            generate_jax_chunking_file.py --chunk_size ${chunk_size} --bed_file ${aggrnorm_counts_bed} --output_prefix gene_chunks
        """
    
}


workflow JAXQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        
    main:

      CHUNK_BED_FILE(condition_bed,params.JAXQTL.number_of_genes_per_chunk,plink_genotype)
      chunking_channel=CHUNK_BED_FILE.out.chunked_bed_channel

      result = chunking_channel.flatMap { item ->
          def (condition,phenotype_file,phenotype_pcs,nr_phenotype_pcs,chunging_file,plink) = item
          if (!(chunging_file instanceof Collection)) {
              chunging_file = [chunging_file] // Wrap single value in a list
          }
          return chunging_file.collect { [condition,phenotype_file,phenotype_pcs,nr_phenotype_pcs,it,plink] }
      }

      JAXQTL(
          result,
          'cis'
      )

      inp_ch2 = JAXQTL.out.qtl_data
        .groupTuple(by: 0)
        .map { cond, files -> tuple(cond, files.unique { it.toString() }) }

      // Combine results and do Qval correction

      AGGREGATE_QTL_RESULTS(inp_ch2) // QTL results are then aggregated.
      all_basic_results = AGGREGATE_QTL_RESULTS.out.jax_qtl_path
        .groupTuple(by: 0)
      // Estimate the OptimPCs
      OPTIM_PCS(all_basic_results)
      optimal_pc_file = OPTIM_PCS.out.optimal_pc_file
      optimal_pc_file
        .map { condition, file ->
            def content = file.text.trim()
            if (content) {
            return ["${condition}__${content}pcs.tsv"]
            } else {
            return null
            }
        }
        .filter { it != null } // Skip if file was empty
        .set { optimal_pc_values }

      results_for_nominal = result.combine(optimal_pc_values,by:0)

      JAXQTL_NOMINAL(
          results_for_nominal,
          'nominal'
      )


      // Run the nominal QTLs.

}