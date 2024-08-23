process CONDITIONAL_QTL {
    label 'process_tiny'

    // Finds top variant per gene and calculates up to 5 conditionally independent signals by including additional variant effects in the model

    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(output),path(output_rda),path(output3),path(for_conditioning),path(genome_regions),path(plink_bim), path(plink_bed), path(plink_fam),val(chr)

    script:
        if (params.SAIGE.cis_trans_mode=='cis'){
            mode="--rangestoIncludeFile=\${step2prefix}_region_file.txt"
        }else{
            mode=""
        }
    """
        echo "~~~~~~~~~~~~~~~~PERFORMING THE SECOND ROUND OF EQTL ANALYSIS~~~~~~~~~~~~~~~~~"
        #

        # Extract the minimum p-value for which there is an FDR < 0.05 in the first round
        #pval_thresh=\$(awk -F'\t' 'NR==2 {print \$13}' ${output3}_minimum_q.txt)        

        # Now perform an additional 3 rounds if we keep finding conditional 

        run_step2_tests_qtl() {
            { 
                topvariant=\$(awk -F'\t' 'NR==2 {print \$3}' ${output3}/cis_\${gene}_${chr}_minimum_q.txt)
                step2_tests_qtl.R       \
                    --bedFile=${plink_bed}      \
                    --bimFile=${plink_bim}      \
                    --famFile=${plink_fam}      \
                    --SAIGEOutputFile=output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}    \
                    --chrom=${chr}       \
                    --minMAF=${params.SAIGE.minMAF} \
                    --minMAC=${params.SAIGE.minMAC} \
                    --LOCO=FALSE    \
                    --varianceRatioFile=\${step1prefix}_\${variable}.varianceRatio.txt    \
                    --GMMATmodelFile=\${step1prefix}_\${variable}.rda    \
                    --SPAcutoff=${params.SAIGE.SPAcutoff} \
                    --condition=\$topvariant \
                    --markers_per_chunk=${params.SAIGE.markers_per_chunk} ${mode} 

                if [ \$? -ne 0 ]; then
                    echo "step2_tests_qtl.R command failed" >&2
                    exit 1  # Exit the script with a non-zero status
                fi

                line_count=\$(wc -l < output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable})        
                if [ "\$line_count" -eq 1 ]; then
                    echo "File has exactly one line"
                    rm output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
                else
                    echo "\${variable}" >> genes_list2.tsv
                fi
                
            } || { 
                echo 'Failed since no markers present in range'
                rm output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
            }
        }


        if awk '\$2 == ${chr} {found=${chr}; exit} END {exit !found}' ${genome_regions}; then
            echo "The chromosome '${chr}' is found in the second column."

            step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
            step2prefix=output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            mkdir output_${name}___${chr}
            
            cat "${genome_regions}" | while IFS= read -r gene || [ -n "\$gene" ]
            do
                echo "\$gene" | cut -f2- >> regions_cis.tsv
                variable=\$(echo "\$gene" | cut -f1)
                echo \${variable}
                chr1=\$(echo "\$gene" | cut -f2)
                
                if [ "\$chr1" -eq ${chr} ]; then
                    run_step2_tests_qtl
                else
                    echo 'Not on the correct chromosome'
                    
                fi
                rm regions_cis.tsv
            done
        else
            echo "The chromosome '${chr}' is not found in the testing range, and hence ignored."
        fi

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
        path("*genes_droped_from_s1_due_to_error.tsv"), emit: dropped optional true

    // Define the Bash script to run for each array job
    script:
    """
        # Execute with the bash executable in an array (one job per gene within level)
        #// Genome wide for this we send a list of genes in chunks 
        mkdir output
        cat "${genes_list}" | while IFS= read -r i || [ -n "\$i" ]
        do
           {  # try
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
                --isCovariateOffset=TRUE  \
                --isCovariateTransform=TRUE  \
                --skipModelFitting=FALSE  \
                --tol=0.00001   \
                --famFile ${plink_fam} \
                --bimFile ${plink_bim} \
                --bedFile ${plink_bed} \
                --IsOverwriteVarianceRatioFile=TRUE 
            } || {
                # catch
                sed -i '/\$i/d' ${genes_list}
                echo \$i >> \${i}_genes_droped_from_s1_due_to_error.tsv
            }
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
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output),path(plink_bim), path(plink_bed), path(plink_fam),val(chr)

    output:
        tuple val("${name}___${chr}"),path(genes_list),path("output_${name}___${chr}"),path(output),path("regions_${genes_list}"),path(plink_bim), path(plink_bed), path(plink_fam),val(chr),emit:output
        tuple val("${name}___${chr}"),path("output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_*"),emit:for_aggregation

    script:
        if (params.SAIGE.cis_trans_mode=='cis'){
            mode="--rangestoIncludeFile=\${step2prefix}_region_file.txt"
        }else{
            mode=""
        }
    
    """
        cp ${genes_list} regions_${genes_list}
        step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
        step2prefix=output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
        mkdir output_${name}___${chr}

        run_step2_tests_qtl() {
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
                    --markers_per_chunk=${params.SAIGE.markers_per_chunk} ${mode}

            if [ \$? -ne 0 ]; then
                echo "step2_tests_qtl.R command failed" >&2
                exit 1  # Exit the script with a non-zero status
            fi

            line_count=\$(wc -l < output/step1_output_formultigenes.txt)
            if [ "\$line_count" -eq 1 ]; then
                echo "File has exactly one line"
                var=\$(cat output/step1_output_formultigenes.txt | awk '{print \$1}')
                echo \$var
                mv output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\$var
                mv output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis.index output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\$var.index
            else
                echo "File does not have exactly one line"
            fi
        }

        run_step2_tests_qtl
        
        cat "${genes_list}" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            f1="output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\$gene"
            first_line=\$(head -n 1 "\$f1")
            # Check if the first line contains both V1 and V2
            if [[ \$first_line == *"V1"* && \$first_line == *"V2"* ]]; then
                echo "The first line contains both V1 and V2."
                run_step2_tests_qtl
            else
                d="\$first_line : The first line does not contain both V1 and V2."
            fi
            line_count=\$(wc -l < \$f1)
            if [ "\$line_count" -eq 1 ]; then
                run_step2_tests_qtl
                echo "The file has more <= 1 line."
            elif [ "\$line_count" -lt 1 ]; then
                run_step2_tests_qtl
                echo "The file has more <= 1 line."
            else
                d="The file has more than 1 line."
            fi
        done


    """
}

process SAIGE_S2_CIS {
    label 'process_low'

    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output),path(genome_regions),path(plink_bim), path(plink_bed), path(plink_fam),val(chr)

    output:
        tuple val("${name}___${chr}"),path("genes_list2.tsv"),path("output_${name}___${chr}"),path(output),path(genome_regions),path(plink_bim), path(plink_bed), path(plink_fam),val(chr),emit:output optional true
        tuple val("${name}___${chr}"),path("output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_*"),emit:for_aggregation  optional true

    script:
        if (params.SAIGE.cis_trans_mode=='cis'){
            mode="--rangestoIncludeFile=regions_cis.tsv"
        }else{
            mode=""
        }
    
    """
        run_step2_tests_qtl() {
            { 
                step2_tests_qtl.R       \
                    --bedFile=${plink_bed}      \
                    --bimFile=${plink_bim}      \
                    --famFile=${plink_fam}      \
                    --SAIGEOutputFile=output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}    \
                    --chrom=${chr}       \
                    --minMAF=${params.SAIGE.minMAF} \
                    --minMAC=${params.SAIGE.minMAC} \
                    --LOCO=FALSE    \
                    --varianceRatioFile=\${step1prefix}_\${variable}.varianceRatio.txt    \
                    --GMMATmodelFile=\${step1prefix}_\${variable}.rda    \
                    --SPAcutoff=${params.SAIGE.SPAcutoff} \
                    --markers_per_chunk=${params.SAIGE.markers_per_chunk} ${mode} 

                if [ \$? -ne 0 ]; then
                    echo "step2_tests_qtl.R command failed" >&2
                    exit 1  # Exit the script with a non-zero status
                fi

                line_count=\$(wc -l < output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable})        
                if [ "\$line_count" -eq 1 ]; then
                    echo "File has exactly one line"
                    rm output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
                else
                    echo "\${variable}" >> genes_list2.tsv
                fi
                
            } || { 
                echo 'Failed since no markers present in range'
                rm output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
            }
        }


        if awk '\$2 == ${chr} {found=${chr}; exit} END {exit !found}' ${genome_regions}; then
            echo "The chromosome '${chr}' is found in the second column."

            step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
            step2prefix=output_${name}___${chr}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            mkdir output_${name}___${chr}
            
            cat "${genome_regions}" | while IFS= read -r gene || [ -n "\$gene" ]
            do
                echo "\$gene" | cut -f2- >> regions_cis.tsv
                variable=\$(echo "\$gene" | cut -f1)
                echo \${variable}
                chr1=\$(echo "\$gene" | cut -f2)
                
                if [ "\$chr1" -eq ${chr} ]; then
                    run_step2_tests_qtl
                else
                    echo 'Not on the correct chromosome'
                    #sed -i '/\$variable/d' ${genes_list}
                fi
                rm regions_cis.tsv
            done
        else
            echo "The chromosome '${chr}' is not found in the testing range, and hence ignored."
        fi

    """
}

process SAIGE_QVAL_COR {
    label 'process_low'

    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        tuple val(name),path(genes_list),path(output),path(output_rda),path(regions),path(plink_bim), path(plink_bed), path(plink_fam),val(chr)

    output:
        tuple val(name),path(genes_list),path(output),emit:output
        tuple val(name),path(output),path(output_rda),path('output3'),path('for_conditioning.csv'),path(regions),path(plink_bim), path(plink_bed), path(plink_fam),val(chr), emit: for_conditioning optional true
        tuple val(name),path("output3/*_minimum_q.txt"), emit: q_out
    script:

        parts = name.split('___')
        chr = parts[-1]
        exp = parts[0]
    """
        
        step2prefix=output_${name}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
        mkdir output3 
        cat "${genes_list}" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            qvalue_correction.py -f \${step2prefix}_\${gene} -c "13" -n "qvalues" -w "TRUE"
            mv \${step2prefix}_\${gene}_minimum_q.txt output3/cis_\${gene}_${chr}_minimum_q.txt

            top_q=\$(awk -F'	' 'NR==2 {print \$17}' output3/cis_\${gene}_${chr}_minimum_q.txt)
            threshold=${params.SAIGE.q_val_threshold_for_conditioning}
            echo "\$top_q"
            if awk -v tq="\$top_q" -v th="\$threshold" 'BEGIN {exit !(tq <= th)}'; then
                echo "Performing conditional analysis: q-value for first pass <= \$threshold"
                echo \${gene},${output_rda}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.rda,${output_rda}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.varianceRatio.txt >> for_conditioning.csv
            fi            
        done
    """
}


process SAIGE_S3 {
    label 'process_tiny'

    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output)
    output:
        tuple val(name),path("output_2/*_cis_genePval"), emit: q_out_2
    // Define the Bash script to run for each array job
    // nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:

    parts = name.split('___')
    chr = parts[-1]

    """
        mkdir output_2
        cat "${genes_list}" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            step3_gene_pvalue_qtl.R \
            --assocFile=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${gene}        \
            --geneName=\${gene}       \
            --genePval_outputFile=output_2/nindep_100_ncell_100_lambda_chr${chr}_tauIntraSample_0.5_\${gene}_cis_genePval
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
            qvalue_correction.py -f minimum_q_all_genes.tsv -c "18" -n "qvalues_across_genes" -w "FALSE"
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
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${exp}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(group),path(qtl_q_results)
    output:
        tuple val(group),path("chr${chr}__all_vars_genes.tsv"), emit: all

    script:
        parts = group.split('___')
        chr = parts[-1]
        exp = parts[0]

        """
            prepend_gene_large.py --pattern 'cis_*' --column 'gene' --outfile 'chr${chr}__all_vars_genes.tsv'
        """
}


process H5AD_TO_SAIGE_FORMAT {
    label 'process_medium'
    tag { sanitized_columns }
    memory { 
            sizeInGB = h5ad.size() / 1e9 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
    
    publishDir  path: "${params.outdir}/Saige_eQTLS",
        saveAs: { filename ->
            if (filename.contains("covariates.txt")) {
                return null
            } else if (filename.contains("saige_filt_expr_input.tsv")) {
                return null
            } else if (filename.contains("test_genes.txt")) {
                return null
            } else if (filename.contains("output_agg/")) {
                // Assuming `filename` contains the full path including directories
                // Remove 'output_agg/azimuth.celltyp.l0/' from the path
                def newFilename = filename.replaceAll(".*output_agg/[^/]+/", "")
                return newFilename
            } else {
                return null
            }
        },
        mode: "${params.copy_mode}",
        overwrite: "true"
    input:
        each path(h5ad)  
        path(bridge)  
        val(aggregation_columns)
        path(genotype_pcs)

    output:
        tuple val(sanitized_columns), path("output_agg/*/*/saige_filt_expr_input.tsv"),path("output_agg/*/*/covariates.txt"),emit:output_pheno optional true
        tuple val(sanitized_columns),path("output_agg/*/*/test_genes.txt"),emit:gene_chunk optional true
        path("output_agg/*"),emit:output_agg optional true

    // Define the Bash script to run for each array job
    script:
    sanitized_columns = h5ad.getName().replaceAll(/[^a-zA-Z0-9]/, '_').replaceAll(/\.h5ad$/, '')
    if("${params.SAIGE.covariate_obs_columns}"==''){
        cov_col = ""
    }else{
        cov_col = "--covariates ${params.SAIGE.covariate_obs_columns}"
    }
    sizeInGB = h5ad.size() / 1e9 * 3 + 5 * task.attempt

    if ("${params.aggregation_subentry}"==''){
        cond1 = " --condition_col 'NULL' --condition 'NULL' "
    }else{
        cond1 = " --condition_col '${aggregation_columns}' --condition '${params.aggregation_subentry}' "
    }

    """
        echo ${sizeInGB}
        bridge='${bridge}'
        nperc=${params.percent_of_population_expressed}
        condition_col="${aggregation_columns}" #Specify 'NULL' if want to include all cells
        condition="${aggregation_columns}" #Specify 'NULL' if want to include all cells
        
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
            --scale_covariates \$scale_covariates \
            --expression_pca \$expression_pca \
            ${cov_col} ${cond1}
    """
}


process TEST {
    label 'process_low'

    // Specify the number of forks (10k)
    // maxForks 1000

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
    tag { sanitized_columns }
    // Specify the number of forks (10k)
    // maxForks 1000

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

process DETERMINE_TSS_AND_TEST_REGIONS {
    label 'process_low'

    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    

    input:
        tuple val(name),path(genes_list),path(output)
        each path(annotation_file)

    output:
        tuple val(name),path(genes_list),path(output),path('gene_regions_to_test.tsv'),emit:output_genes
        
    // We want to test in a ceirtain window of the gene TSS +/-, to do this we need to propeare SAIGE regions files. 
    script:
    """
        echo 'Lets prepeare the regions file'
        prepeare_Saige_regions_for_cis.py --annotation_file ${annotation_file} --genes ${genes_list} --gtf_type ${params.gtf_type} --window ${params.windowSize}
    """
}

workflow SAIGE_qtls{
    take:
        genotype_pcs
        phenotype_file
        bim_bed_fam
        genome_annotation

    main:
        log.info('------- Running SAIGE QTLs ------- ')

        H5AD_TO_SAIGE_FORMAT(phenotype_file,params.genotype_phenotype_mapping_file,params.aggregation_columns,genotype_pcs)
        pheno = H5AD_TO_SAIGE_FORMAT.out.output_pheno

        H5AD_TO_SAIGE_FORMAT.out.gene_chunk.subscribe { println "H5AD_TO_SAIGE_FORMAT.out.gene_chunk dist: $it" }
        CHUNK_GENES(H5AD_TO_SAIGE_FORMAT.out.gene_chunk,params.chunkSize)
        result = CHUNK_GENES.out.output_genes.flatMap { item ->
            def (first, second) = item
            return second.collect { [first, it] }
        }

        result.combine(pheno, by: 0).set{pheno_chunk}

        
        Channel.fromList(params.SAIGE.chromosomes_to_test)
                .set{chromosomes_to_test}        

        SAIGE_S1(pheno_chunk.combine(bim_bed_fam))


        if(params.SAIGE.cis_trans_mode=='trans'){
            SAIGE_S2(SAIGE_S1.out.output.combine(bim_bed_fam).combine(chromosomes_to_test))
            output_s2 = SAIGE_S2.out.output
            agg_output = SAIGE_S2.out.for_aggregation

        }else if(params.SAIGE.cis_trans_mode=='cis'){
            DETERMINE_TSS_AND_TEST_REGIONS(SAIGE_S1.out.output,genome_annotation)
            for_cis_input = DETERMINE_TSS_AND_TEST_REGIONS.out.output_genes
            SAIGE_S2_CIS(for_cis_input.combine(bim_bed_fam).combine(chromosomes_to_test))
            output_s2 = SAIGE_S2_CIS.out.output
            agg_output = SAIGE_S2_CIS.out.for_aggregation
        }

        // HERE WE either run the cis or trans qtl mapping. For cis we loop through each of the chunks whereas in trans we can run all together.
         
        SAIGE_QVAL_COR(output_s2)
        SAIGE_S3(SAIGE_QVAL_COR.out.output)


        // SAIGE_QVAL_COR.out.for_conditioning.subscribe { println "SAIGE_QVAL_COR dist: $it" }
        // ########## Collecting Chunk outputs.  ###############
        SAIGE_S2_for_aggregation = agg_output.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { [first, it] }
        }

        SAIGE_S3_for_aggregation = SAIGE_QVAL_COR.out.q_out.flatMap { item ->
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
        SAIGE_S3_for_aggregation_ACAT = SAIGE_S3_for_aggregation_ACAT.map{row->tuple("${row[0]}".replaceFirst(/___.*/,""),
                                                    file(row[1])
                                                    )}  

        SAIGE_S3_for_aggregation = SAIGE_S3_for_aggregation.map{row->tuple("${row[0]}".replaceFirst(/___.*/,""),
                                                    file(row[1])
                                                    )}  
        // SAIGE_S3_for_aggregation_ACAT.subscribe { println "SAIGSAIGE_S3_for_aggregation_ACATE_QVAL_COR dist: $it" }
        // ########## Aggregating and emiting the results.  ###############
        AGGREGATE_QTL_RESULTS(SAIGE_S3_for_aggregation.groupTuple(by: 0))
        AGGREGATE_QTL_ALLVARS(SAIGE_S2_for_aggregation.groupTuple(by: 0))
        AGGREGATE_ACAT_RESULTS(SAIGE_S3_for_aggregation_ACAT.groupTuple(by: 0))
        // CONDITIONAL_QTL(SAIGE_QVAL_COR.out.for_conditioning)

}