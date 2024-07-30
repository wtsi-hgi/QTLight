process CONDITIONAL_QTL {
    label 'process_tiny'

    // Finds top variant per gene and calculates up to 5 conditionally independent signals by including additional variant effects in the model

    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(for_conditioning)
    // output:
    //     tuple val(name),path("${output}/*_minimum_q.txt"),path("${output}/*_cis_genePval"), emit: q_out
        // echo step2_tests_qtl.R \
        //         --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
        //         --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
        //         --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
        //         --SAIGEOutputFile=${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round2.txt    \
        //         --chrom=${gene_chr}       \
        //         --minMAF=0.05 \
        //         --minMAC=20 \
        //         --LOCO=FALSE    \
        //         --is_imputed_data=TRUE \
        //         --GMMATmodelFile=${step1prefix}.rda     \
        //         --SPAcutoff=2 \
        //         --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        //         --rangestoIncludeFile=${step2prefix}_region_file.txt     \
        //         --markers_per_chunk=10000 \
        //         --condition=$topvariant
        // echo qvalue_correction.R -f ${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round2.txt -c "20" -n "c-${topvariant}.qvalues" -w "TRUE"
    script:
    """
        echo "~~~~~~~~~~~~~~~~PERFORMING THE SECOND ROUND OF EQTL ANALYSIS~~~~~~~~~~~~~~~~~"


        # Now perform an additional 3 rounds if we keep finding conditional 

    """
}


process SAIGE_S1 {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    


    input:
        tuple val(name),path(genes_list),path(pheno_file),path(cov),path(plink_bim), path(plink_bed), path(plink_fam)
        
    output:
        tuple val(name),path(genes_list),path("output"),emit:output

    // Define the Bash script to run for each array job
    script:
    """
        # Execute with the bash executable in an array (one job per gene within level)
        #// Genome wide for this we send a list of genes in chunks 
        mkdir output
        cat "${genes_list}" | while IFS= read -r i || [ -n "\$i" ]
        do
            step1_fitNULLGLMM_qtl.R \
                --useSparseGRMtoFitNULL=FALSE  \
                --useGRMtoFitNULL=FALSE \
                --phenoFile=${pheno_file}	\
                --phenoCol=\$i       \
                --covarColList=\$(head -n 1 ${cov})    \
                --sampleCovarColList=\$(sed -n '2p' ${cov})      \
                --sampleIDColinphenoFile=${params.gt_id_column} \
                --traitType=count \
                --outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\$i  \
                --skipVarianceRatioEstimation=FALSE  \
                --isRemoveZerosinPheno=FALSE \
                --isCovariateOffset=FALSE  \
                --isCovariateTransform=TRUE  \
                --skipModelFitting=FALSE  \
                --tol=0.00001   \
                --famFile ${plink_fam} \
                --bimFile ${plink_bim} \
                --bedFile ${plink_bed} \
                --IsOverwriteVarianceRatioFile=TRUE
        done

        cat "${genes_list}" | while IFS= read -r i || [ -n "\$i" ]
        do
            echo \$i >> ./output/gene_string.tsv
            step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\$i
            echo -e "\${i} \${step1prefix}.rda \${step1prefix}.varianceRatio.txt" >> ./output/step1_output_formultigenes.txt
        done
    """
}

process SAIGE_S2 {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output),path(plink_bim), path(plink_bed), path(plink_fam),val(chr)

    output:
        tuple val("${name}___${chr}"),path(genes_list),path(output),emit:output
        tuple val("${name}___${chr}"),path("${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_*"),emit:for_aggregation

    script:
    """
        
        step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
        step2prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis

        step2_tests_qtl.R       \
                --bedFile=${plink_bed}      \
                --bimFile=${plink_bim}      \
                --famFile=${plink_fam}      \
                --SAIGEOutputFile=\${step2prefix}     \
                --chrom=${chr}       \
                --minMAF=${params.SAIGE.minMAF} \
                --minMAC=${params.SAIGE.minMAC} \
                --LOCO=FALSE    \
                --GMMATmodel_varianceRatio_multiTraits_File=${output}/step1_output_formultigenes.txt     \
                --SPAcutoff=${params.SAIGE.SPAcutoff} \
                --markers_per_chunk=${params.SAIGE.markers_per_chunk}


        line_count=\$(wc -l < output/step1_output_formultigenes.txt)
        if [ "\$line_count" -eq 1 ]; then
            echo "File has exactly one line"
            var=\$(cat output/step1_output_formultigenes.txt | awk '{print \$1}')
            echo \$var
            mv output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\$var
            mv output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis.index output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\$var.index
        else
            echo "File does not have exactly one line"
        fi

    """
}

process SAIGE_QVAL_COR {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output)

    output:
        tuple val(name),path(genes_list),path(output),emit:output
        tuple val(name),path('for_conditioning.csv'), emit: for_conditioning optional true

    script:
    """
        step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
        step2prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
                
        cat "${output}/gene_string.tsv" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            qvalue_correction.R -f \${step2prefix}_\${gene} -c "13" -n "qvalues" -w "TRUE"
            mv \${step2prefix}_\${gene}_minimum_q.txt ${output}/cis_\${gene}_minimum_q.txt

            top_q=\$(awk -F'\t' 'NR==2 {print \$18}' ${output}/cis_\${gene}_minimum_q.txt)
            threshold=0.05
            if (( \$(echo "\$top_q <= \$threshold" | bc -l) )); then
                echo "Performing conditional analysis: q-value for first pass > \$threshold"
                echo \${gene},${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.rda,nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.varianceRatio.txt >> for_conditioning.csv
            fi            
        done
    """
}


process SAIGE_S3 {
    label 'process_tiny'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output)
    output:
        tuple val(name),path("${output}/*_minimum_q.txt"), emit: q_out
        tuple val(name),path("${output}/*_cis_genePval"), emit: q_out_2
    // Define the Bash script to run for each array job
    // nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:
    """
        cat "${output}/gene_string.tsv" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            step3_gene_pvalue_qtl.R \
            --assocFile=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${gene}        \
            --geneName=\${gene}       \
            --genePval_outputFile=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}_cis_genePval
        done
    """
}


process AGGREGATE_ACAT_RESULTS{
    tag { condition }
    scratch false      // use tmp directory
    label 'process_tiny'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${group}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(group),path(cis_genePval)
    output:
        tuple val(group),path("ACAT_all.tsv"), emit: acat_all
    script:
        
        """
            head -n 1 \$(find . -name "*_cis_genePval" | head -n 1)  >> ACAT_all.tsv
            find . -name "*_cis_genePval" -print0 | xargs -0 -I {} sh -c 'tail -n +2 "\$1"' sh {} >> ACAT_all.tsv
        """
}


process AGGREGATE_QTL_RESULTS{
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${group}",
                mode: "${params.copy_mode}",
                overwrite: "true"


    input:
        tuple val(group),path(qtl_q_results)
    output:
        tuple val(group),path("minimum_q_all_genes.tsv"), emit: all

    script:
        
        """
            prepend_gene.py --pattern 'cis_*_minimum_q.txt' --column 'gene' --outfile 'minimum_q_all_genes.tsv'
            qvalue_correction.R -f minimum_q_all_genes.tsv -c "18" -n "qvalues_across_genes" -w "FALSE"
        """
}

process AGGREGATE_QTL_ALLVARS{
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${group}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(group),path(qtl_q_results)
    output:
        tuple val(group),path("all_vars_genes.tsv"), emit: all

    script:
        
        """
            prepend_gene_large.py --pattern 'cis_*' --column 'gene' --outfile 'all_vars_genes.tsv'
        """
}


process H5AD_TO_SAIGE_FORMAT {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    


    input:
        each path(h5ad)  
        path(bridge)  
        val(aggregation_columns)
        path(genotype_pcs)

    output:
        tuple val(sanitized_columns), path("output_agg/*/*/saige_filt_expr_input.tsv"),path("output_agg/*/*/covariates.txt"),emit:output_pheno optional true
        tuple val(sanitized_columns),path("output_agg/*/*/test_genes.txt"),emit:gene_chunk optional true

    // Define the Bash script to run for each array job
    script:
    sanitized_columns = h5ad.getName().replaceAll(/[^a-zA-Z0-9]/, '_').replaceAll(/\.h5ad$/, '')
    """
        echo ${sanitized_columns}
        bridge='${bridge}'
        nperc=${params.percent_of_population_expressed}
        condition_col="NULL" #Specify 'NULL' if want to include all cells
        covariates="total_counts"
        scale_covariates=true
        expression_pca=${params.SAIGE.nr_expression_pcs}
        aggregate_on="${aggregation_columns}"

        mkdir output_agg
        prep_adata_saige.py \
            --phenotype__file ${h5ad} \
            --bridge \$bridge \
            --aggregate_on \$aggregate_on \
            --genotype_pc__file ${genotype_pcs} \
            --genotype_id ${params.gt_id_column} \
            --sample_id ${params.sample_column} \
            --general_file_dir ./output_agg \
            --nperc \$nperc \
            --min ${params.n_min_cells} \
            --condition_col \$condition_col \
            --condition \$condition_col \
            --covariates \$covariates \
            --scale_covariates \$scale_covariates \
            --expression_pca \$expression_pca
    """
}


process TEST {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(sanitized_columns), path(saige_filt_expr_input),path(test_genes) 
       
    // Define the Bash script to run for each array job
    // nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:
    """
        echo ${sanitized_columns}
        echo ${test_genes}
    """
}


process CHUNK_GENES {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(sanitized_columns),path(test_genes) 
        val(chunk_size)
    output:
        tuple val(sanitized_columns), path("chunk_${sanitized_columns}_*"),emit:output_genes
        
    // Define the Bash script to run for each array job
    // nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:
    """
        echo ${sanitized_columns}
        echo ${test_genes}
        split -l ${chunk_size} ${test_genes} chunk_${sanitized_columns}_
    """
}

workflow SAIGE_qtls{
    take:
        genotype_pcs
        phenotype_file
        bim_bed_fam

    main:
        log.info('------- Running SAIGE QTLs ------- ')

        H5AD_TO_SAIGE_FORMAT(phenotype_file.flatten(),params.genotype_phenotype_mapping_file,params.aggregation_columns,genotype_pcs)
        pheno = H5AD_TO_SAIGE_FORMAT.out.output_pheno.take(2)

        // We also need to define the chromosomes to test. 
        // If we go for a cis mode then we have to also define the ranges.
        // genes = Channel.of(['t1','/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/genes.txt'])
        // pheno = Channel.of(['t1','/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt','/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/covariates.txt'])
        // bim_bed_fam = Channel.of(["/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.bim","/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.bed","/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.fam"])

        genes = H5AD_TO_SAIGE_FORMAT.out.gene_chunk
        CHUNK_GENES(genes,params.chunkSize)
        result = CHUNK_GENES.out.output_genes.flatMap { item ->
            def (first, second) = item
            return second.collect { [first, it] }
        }

        result.combine(pheno, by: 0).set{pheno_chunk}
        SAIGE_S1(pheno_chunk.combine(bim_bed_fam))
        Channel.fromList(params.SAIGE.chromosomes_to_test)
                .set{chromosomes_to_test}
        SAIGE_S2(SAIGE_S1.out.output.combine(bim_bed_fam).combine(chromosomes_to_test))
        SAIGE_QVAL_COR(SAIGE_S2.out.output)
        SAIGE_S3(SAIGE_QVAL_COR.out.output)


        SAIGE_QVAL_COR.out.for_conditioning.subscribe { println "SAIGE_QVAL_COR dist: $it" }
        // ########## Collecting Chunk outputs.  ###############
        SAIGE_S2_for_aggregation = SAIGE_S2.out.for_aggregation.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { [first, it] }
        }

        SAIGE_S3_for_aggregation = SAIGE_S3.out.q_out.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { [first, it] }   
        }

        SAIGE_S3_for_aggregation_ACAT = SAIGE_S3.out.q_out_2.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { [first, it] }
        }

        // ########## Aggregating and emiting the results.  ###############
        AGGREGATE_QTL_RESULTS(SAIGE_S3_for_aggregation.groupTuple(by: 0))
        AGGREGATE_QTL_ALLVARS(SAIGE_S2_for_aggregation.groupTuple(by: 0))
        AGGREGATE_ACAT_RESULTS(SAIGE_S3_for_aggregation_ACAT.groupTuple(by: 0))
        CONDITIONAL_QTL(SAIGE_QVAL_COR.out.for_conditioning)

}