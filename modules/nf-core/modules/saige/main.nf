process RUNSAIGE {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    
    // Input files for the array job
    input:
        levelOption from PREPADATACATEGORY.out
        phenotype__file
        aggregate_on
        general_file_dir
        n_geno_pcs
        covariates
        covariates_cell
        genotype_id
        sample_id
        annotation__file
        condition_col
        condition
        cis_only
        cis_window
        knee
        gene from path("$general_file_dir/$aggregate_on/$levelOption/test_genes.txt")

    // Output files
    output:
        path("${params.outdir}/${general_file_dir}/$aggregate_on/$levelOption/chr*_nPC_$n_geno_pcs.txt")

    // Define the Bash script to run for each array job
    script:
    """
        # Execute with the bash executable in an array (one job per gene within level)
        bash bin/RUNSAIGE_1_2.sh 
            -c levelOption 
            -p $phenotype__file 
            -a $aggregate_on 
            -d $general_file_dir 
            -w $n_geno_pcs 
            -e $covariates 
            -k $covariates_cell 
            -i $genotype_id 
            -s $sample_id 
            -m $annotation__file 
            -o $condition_col 
            -t $condition 
            -x $cis_only 
            -y $cis_window 
            -k $knee 
            -g gene
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
        path(h5ad)  
        path(bridge)  
        val(aggregation_columns)
        path(genotype_pcs)
    // output:
    //     path("output"),emit:output

    // Define the Bash script to run for each array job
    script:
    """
        bridge='${bridge}'
        nperc=${params.percent_of_population_expressed}
        condition_col="NULL" #Specify 'NULL' if want to include all cells
        covariates="total_counts"
        scale_covariates=false
        expression_pca="true"
        aggregate_on="${aggregation_columns}"
        level="B"

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
            --expression_pca \$expression_pca \
            --level \$level
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


    output:
        path("output"),emit:output

    // Define the Bash script to run for each array job
    script:
    """
        # Execute with the bash executable in an array (one job per gene within level)
        mkdir output
        step1_fitNULLGLMM_qtl.R \
            --useSparseGRMtoFitNULL=FALSE  \
            --useGRMtoFitNULL=FALSE \
            --phenoFile=/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt	\
            --phenoCol=gene_1       \
            --covarColList=X1,X2,pf1,pf2    \
            --sampleCovarColList=X1,X2      \
            --sampleIDColinphenoFile=IND_ID \
            --traitType=count \
            --outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1 \
            --skipVarianceRatioEstimation=FALSE  \
            --isRemoveZerosinPheno=FALSE \
            --isCovariateOffset=FALSE  \
            --isCovariateTransform=TRUE  \
            --skipModelFitting=FALSE  \
            --tol=0.00001   \
            --plinkFile=/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1_01.step1       \
            --IsOverwriteVarianceRatioFile=TRUE
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
        path(output)   


    output:
        path("output"),emit:output

    // Define the Bash script to run for each array job
    script:
    """
        export regionFile=gene_1_cis_region.txt
        echo -e "2\t300001\t610001" > \${regionFile}
        step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1
        step2prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis

        step2_tests_qtl.R       \
                --bedFile=/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.bed      \
                --bimFile=/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.bim      \
                --famFile=/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/qtl/extdata/input/n.indep_100_n.cell_1.fam      \
                --SAIGEOutputFile=\${step2prefix}     \
                --chrom=2       \
                --minMAF=0 \
                --minMAC=20 \
                --LOCO=FALSE    \
                --GMMATmodelFile=\${step1prefix}.rda     \
                --SPAcutoff=2 \
                --varianceRatioFile=\${step1prefix}.varianceRatio.txt    \
                --rangestoIncludeFile=\${regionFile}     \
                --markers_per_chunk=10000
        
    """
}


process SAIGE_S3 {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        path(output)   

    // Define the Bash script to run for each array job
    script:
    """
        step3_gene_pvalue_qtl.R \
        --assocFile=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis        \
        --geneName=gene_1       \
        --genePval_outputFile=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval
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
    

    publishDir  path: "${params.outdir}/Limix_eQTLS",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path(all_qtl_results)
        
        
    output:
        path("${all_qtl_results}_qtls"), emit: limix_qtl_path
        // path("${all_qtl_results}_all_mpr/qtl_results_all.txt"), emit: qtl_results_all_all_mpr
        // path("${all_qtl_results}_all_mpr/top_qtl_results_all.txt"), emit: top_qtl_results_all_all_mpr
        // path("${all_qtl_results}_qtls")
    script:
        
        """
            export NUMBA_CACHE_DIR=/tmp
            export MPLCONFIGDIR=/tmp
            mkdir ${all_qtl_results}_all
            
            minimal_postprocess.py -id ${all_qtl_results} -od ${all_qtl_results}_all -sfo -tfb 
            minimal_postprocess.py -id ${all_qtl_results} -od ${all_qtl_results}_all -sfo -mrp 0.05 

            cp -Lr ${all_qtl_results}_all ${all_qtl_results}_qtls
        """
}

workflow SAIGE_qtls{
    take:
        genotype_pcs
    //     limix_condition_chunking
    //     plink_genotype
    //     kinship_file
    //     genotype_phenotype_mapping_file
    //     condition

    main:
        log.info('------- Running SAIGE QTLs ------- ')

        H5AD_TO_SAIGE_FORMAT(params.phenotype_file,params.bridge,param.aggregation_columns,genotype_pcs)
        SAIGE_S1()
        SAIGE_S2(SAIGE_S1.out.output)
        SAIGE_S3(SAIGE_S2.out.output)

}