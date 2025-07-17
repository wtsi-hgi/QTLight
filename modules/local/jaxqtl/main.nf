
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
        tuple val(), path("merged_cis_scores.tsv.gz"), emit: limix_qtl_path

    script:
        all_qtl_results = "${condition}".replaceFirst(/.*\.tsv$/, '')

        group0 = "${condition}".split('__')[0]
        group4 = "${condition}".split('__')[4]
        group1 =  "${condition}".split('__')[1]
        group2 = "${condition}".split('__')[2].replaceFirst(/\.tsv/, '')

        """
            echo "${condition} ${group1} ${group2} ${group0}"
            merge_chunks.py
        """
}



process JAXQTL {  
  tag "$condition, $nr_phenotype_pcs"
  label "process_high_memory"
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
        val(nr_phenotype_pcs),
        path(genelist)
    )
    each path(plink_files_prefix)

  output:
    tuple val("${condition}__${nr_phenotype_pcs}"), path('*cis_score.tsv.gz'), emit: qtl_data

  script:
    outpath = "${nr_phenotype_pcs}/base_output/base"
    """
        mkdir -p ${outpath}
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
        --mode cis \
        --window ${params.windowSize} \
        --test-method score \
        --nperm ${params.numberOfPermutations} \
        --addpc 0 \
        --standardize \
        -p cpu \
        --out \$outname

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

    output:
        tuple val(condition),path(aggrnorm_counts_bed), path(covariates_tsv),val(nr_phenotype_pcs),path("gene_chunks*.tsv"), emit: chunked_bed_channel optional true

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
      // 
    //   condition_bed.subscribe { println "condition_bed: $it" }

      CHUNK_BED_FILE(condition_bed,params.JAXQTL.number_of_genes_per_chunk)
      chunking_channel=CHUNK_BED_FILE.out.chunked_bed_channel
    //   chunking_channel.subscribe { println "chunking_channel: $it" }
      result = chunking_channel.flatMap { item ->
          def (condition,phenotype_file,phenotype_pcs,nr_phenotype_pcs,chunging_file) = item
          if (!(chunging_file instanceof Collection)) {
              chunging_file = [chunging_file] // Wrap single value in a list
          }
          return chunging_file.collect { [condition,phenotype_file,phenotype_pcs,nr_phenotype_pcs,it] }
      }

      JAXQTL(
          result,
          plink_genotype
      )

      inp_ch2 = JAXQTL.out.qtl_data
        .groupTuple(by: 0)
        .map { cond, files -> tuple(cond, files.unique { it.toString() }) }
    //   inp_ch2.subscribe { println "inp_ch2: $it" }
      // Combine results and do Qval correction

      AGGREGATE_QTL_RESULTS(inp_ch2) // QTL results are then aggregated.
      // MULTIPLE_TESTING_CORRECTION(AGGREGATE_QTL_RESULTS.out.limix_qtl_path) // Multiple testing is performed.

      // Estimate the OptimPCs

      // Run the nominal QTLs.

}