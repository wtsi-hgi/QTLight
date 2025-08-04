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
                    --SAIGEOutputFile=output_${name}___${chr}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}    \
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

                line_count=\$(wc -l < output_${name}___${chr}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable})        
                if [ "\$line_count" -eq 1 ]; then
                    echo "File has exactly one line"
                    rm output_${name}___${chr}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
                else
                    echo "\${variable}" >> genes_list2.tsv
                fi
                
            } || { 
                echo 'Failed since no markers present in range'
                rm output_${name}___${chr}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
            }
        }


        if awk '\$2 == ${chr} {found=${chr}; exit} END {exit !found}' ${genome_regions}; then
            echo "The chromosome '${chr}' is found in the second column."

            step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
            step2prefix=output_${name}___${chr}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            mkdir -p output_${name}___${chr}
            
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

process CREATE_SPARSE_GRM {

    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_grm_container}"
    } else {
        container "${params.saige_grm_docker}"
    }    


    publishDir  path: "${params.outdir}/SPARSE_GRM",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path(plink)
        
    output:
        // tuple val(name),path(genes_list),path("output"),emit:output
        path("sparseGRM_*.mtx"), emit: sparseGRM
        path("sparseGRM_*.sampleIDs.txt"), emit: sparseGRM_sample

    // Define the Bash script to run for each array job
    script:
    """

        plink_dir="${plink}"
        base_name=""

        if ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
        elif ls "\$plink_dir"/*.bim 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bim | head -n 1) .bim)
        else
            echo "No .bed or .bim file found in \$plink_dir"
            exit 1
        fi

        echo "Detected base name: \$base_name"


        createSparseGRM.R \
            --nThreads=${task.cpus} \
            --outputPrefix=sparseGRM_output \
            --numRandomMarkerforSparseKin=${params.SAIGE.numRandomMarkerforSparseKin} \
            --relatednessCutoff ${params.SAIGE.relatednessCutoff} \
            --famFile  "\$plink_dir/\$base_name.fam" \
            --bimFile  "\$plink_dir/\$base_name.bim" \
            --bedFile  "\$plink_dir/\$base_name.bed" 
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
        tuple val(name),path(genes_list),path(pheno_file),path(cov),path(plink)
        each path(sparseGRM)
        each path(sparseGRM_samples)
    output:
        tuple val(name),path(genes_list),path("output"),emit:output
        path("*genes_droped_from_s1_due_to_error.tsv"), emit: dropped optional true

    // Define the Bash script to run for each array job
    script:
    """

        plink_dir="${plink}"
        base_name=""

        if ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
        elif ls "\$plink_dir"/*.bim 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bim | head -n 1) .bim)
        else
            echo "No .bed or .bim file found in \$plink_dir"
            exit 1
        fi

        echo "Detected base name: \$base_name"



        # Execute with the bash executable in an array (one job per gene within level)
        #// Genome wide for this we send a list of genes in chunks 
        mkdir -p output
        cat "${genes_list}" | while IFS= read -r i || [ -n "\$i" ]
        do
           {  # try
                step1_fitNULLGLMM_qtl.R \
                --useSparseGRMtoFitNULL=${params.SAIGE.useSparseGRMtoFitNULL} \
                --useSparseGRMforVarRatio=${params.SAIGE.useSparseGRMforVarRatio} \
                --isCateVarianceRatio=${params.SAIGE.isCateVarianceRatio} \
                --cateVarRatioMinMACVecExclude=${params.SAIGE.cateVarRatioMinMACVecExclude} \
                --cateVarRatioMaxMACVecInclude=${params.SAIGE.cateVarRatioMaxMACVecInclude} \
                --numRandomMarkerforVarianceRatio=${params.SAIGE.numRandomMarkerforVarianceRatio} \
                --relatednessCutoff=${params.SAIGE.relatednessCutoff} \
                --skipModelFitting=${params.SAIGE.skipModelFitting} \
                --skipVarianceRatioEstimation=${params.SAIGE.skipVarianceRatioEstimation} \
                --isCovariateTransform=${params.SAIGE.isCovariateTransform} \
                --isCovariateOffset=${params.SAIGE.isCovariateOffset} \
                --isRemoveZerosinPheno=${params.SAIGE.isRemoveZerosinPheno} \
                --tol=${params.SAIGE.tol} \
                --traceCVcutoff=${params.SAIGE.traceCVcutoff} \
                --nrun=${params.SAIGE.nrun} \
                --sparseGRMFile ${sparseGRM} \
                --sparseGRMSampleIDFile ${sparseGRM_samples} \
                --phenoFile=${pheno_file} \
                --phenoCol=\$i \
                --covarColList=\$(head -n 1 ${cov}) \
                --sampleCovarColList=\$(sed -n '2p' ${cov}) \
                --sampleIDColinphenoFile='${params.gt_id_column}' \
                --traitType=count \
                --outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\$i \
                --famFile "\$plink_dir/\$base_name.fam" \
                --bimFile "\$plink_dir/\$base_name.bim" \
                --bedFile "\$plink_dir/\$base_name.bed" \
                --IsOverwriteVarianceRatioFile=TRUE \
                --maxiterPCG=${params.SAIGE.maxiterPCG} --invNormalize=${params.SAIGE.invNormalize} --minCovariateCount=${params.SAIGE.minCovariateCount} \
                --minMAFforGRM=${params.SAIGE.minMAFforGRM} \
                --maxMissingRateforGRM=${params.SAIGE.maxMissingRateforGRM} \
                --useGRMtoFitNULL=${params.SAIGE.useGRMtoFitNULL} \
                ${params.SAIGE.step1_extra_flags}
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


process SAIGE_S2_CIS {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_container}"
    } else {
        container "${params.saige_docker}"
    }    

    input:
        tuple (val(name),path(genes_list),path(output),path(genome_regions),val(plink))
        each path(sparseGRM)
        each path(sparseGRM_samples)

    output:
        tuple val("${name}"),path("genes_list2.tsv"),path("output_${name}___*"),path(output),path(genome_regions),val(plink),emit:output optional true
        tuple val("${name}"),path("output_${name}___*/*___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_*"),emit:for_aggregation  optional true

    script:
        if (params.SAIGE.cis_trans_mode=='cis'){
            mode="--rangestoIncludeFile=regions_cis.tsv --chrom=\${chr1}"
        }else{
            if (params.SAIGE.trans_chr_to_test=='SAME'){
                mode="--chrom=\${chr1}"
            }else{
                mode="--chrom=${params.SAIGE.trans_chr_to_test}"
            }
            
        }
        plink_base = "${plink}".replaceFirst(/\.bgen$/, '')

        if ("${plink}".contains(".bgen")) {
            base_run = ""
            genotype_input = "--bgenFile ${plink} --bgenFileIndex ${plink}.bgi --sampleFile ${plink_base}.sample \\"
        }else{
            base_run =  """plink_dir="${plink}"
                            base_name=""
                            if ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
                                base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
                            elif ls "\$plink_dir"/*.bim 1> /dev/null 2>&1; then
                                base_name=\$(basename \$(ls "\$plink_dir"/*.bim | head -n 1) .bim)
                            else
                                echo "No .bed or .bim file found in \$plink_dir"
                                exit 1
                            fi
                            echo "Detected base name: \$base_name" """

            genotype_input = """--bedFile="\$plink_dir/\$base_name.bed" \\
                            --bimFile="\$plink_dir/\$base_name.bim" \\
                            --famFile="\$plink_dir/\$base_name.fam" \\"""
        }

    """

        ${base_run}

        run_step2_tests_qtl() {
            { 
                warning_output_file=\$(mktemp)
                step2_tests_qtl.R \
                ${genotype_input}
                --sparseGRMFile ${sparseGRM} --sparseGRMSampleIDFile ${sparseGRM_samples} --relatednessCutoff=${params.SAIGE.relatednessCutoff} \
                --SAIGEOutputFile=output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable} \
                --minMAF=${params.SAIGE.minMAF} \
                --minMAC=${params.SAIGE.minMAC} \
                --LOCO=${params.SAIGE.LOCO} \
                --varianceRatioFile=\${step1prefix}_\${variable}.varianceRatio.txt \
                --GMMATmodelFile=\${step1prefix}_\${variable}.rda    \
                --SPAcutoff=${params.SAIGE.SPAcutoff} \
                --markers_per_chunk=${params.SAIGE.markers_per_chunk} ${mode}  > "\$warning_output_file" 2>&1
                exit_code=\$?
                warning_output=\$(cat "\$warning_output_file")
                rm "\$warning_output_file"

                if [ \$exit_code -ne 0 ]; then
                    echo "step2_tests_qtl.R command failed" >&2
                    echo "\$warning_output"
                    return
                fi

                if [[ "\$warning_output" != *"Input/output error"* ]]; then
                    echo "proceed"
                else
                    echo "warning exists"
                    return 1
                fi

                line_count=\$(wc -l < output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable})        
                if [ "\$line_count" -eq 1 ]; then
                    echo "File has exactly one line"
                    rm output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
                elif [ "\$line_count" -eq 0 ]; then
                    echo "File has exactly zero lines"
                    return 1
                else
                    echo "\${variable}\t\${chr1}" >> genes_list2.tsv
                fi
                
            } || { 
                echo 'Failed since no markers present in range'
                rm output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_\${variable}
            }
        }


        step1prefix=${output}/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5           
        
        cat "${genome_regions}" | while IFS= read -r gene || [ -n "\$gene" ]
        do
            echo "\$gene" | cut -f2- >> regions_cis.tsv
            variable=\$(echo "\$gene" | cut -f1)
            echo \${variable}
            chr1=\$(echo "\$gene" | cut -f2)
            step2prefix=output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            
            mkdir -p output_${name}___\${chr1}
            run_step2_tests_qtl
            rm regions_cis.tsv
        done


    """
}

process SAIGE_QVAL_COR {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        tuple val(name),path(genes_list),path(output),path(output_rda),path(regions),path(plink)

    output:
        tuple val(name),path(genes_list),path(output),emit:output
        tuple val(name),path(output),path(output_rda),path('output3'),path('for_conditioning.csv'),path(regions),path(plink), emit: for_conditioning optional true
        tuple val(name),path("output3/*_minimum_q.txt"), emit: q_out
    script:

        parts = name.split('___')
        chr = parts[-1]
        exp = parts[0]
    """
        
        
        mkdir -p output3 
        cat "${genes_list}" | while IFS= read -r gene1 || [ -n "\$gene1" ]
        do
            chr1=\$(echo "\${gene1}" | cut -f2)
            gene=\$(echo "\${gene1}" | cut -f1)
            step2prefix=output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            qvalue_correction.py -f \${step2prefix}_\${gene} -c "13" -n "qvalues" -w "TRUE"
            mv \${step2prefix}_\${gene}_minimum_q.txt output3/cis_\${gene}_\${chr1}_minimum_q.txt

            top_q=\$(awk -F'	' 'NR==2 {print \$17}' output3/cis_\${gene}_\${chr1}_minimum_q.txt)
            threshold=${params.SAIGE.q_val_threshold_for_conditioning}
            echo "\$top_q"
            if awk -v tq="\$top_q" -v th="\$threshold" 'BEGIN {exit !(tq <= th)}'; then
                echo "Performing conditional analysis: q-value for first pass <= \$threshold"
                echo \${gene},${output_rda}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.rda,${output_rda}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_\${gene}.varianceRatio.txt >> for_conditioning.csv
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
    // \${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:

    parts = name.split('___')
    chr = parts[-1]

    """
        mkdir output_2
        cat "${genes_list}" | while IFS= read -r gene1 || [ -n "\$gene1" ]
        do
            chr1=\$(echo "\${gene1}" | cut -f2)
            gene=\$(echo "\${gene1}" | cut -f1)
            step2prefix=output_${name}___\${chr1}/\${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis
            step3_gene_pvalue_qtl.R \
            --assocFile=\${step2prefix}_\${gene}        \
            --geneName=\${gene}       \
            --genePval_outputFile=output_2/nindep_100_ncell_100_lambda_chr\${chr1}_tauIntraSample_0.5_\${gene}_cis_genePval
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
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${exp}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(group),path(cis_genePval)
    output:
        tuple val(group),path("ACAT_all.tsv"), emit: acat_all
    script:
        parts = group.split('__')
        chr = parts[-1]
        exp = parts[-2] 
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
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${exp}",
                mode: "${params.copy_mode}",
                overwrite: "true"


    input:
        tuple val(group),path(qtl_q_results)
    output:
        tuple val(group),path("minimum_q_all_genes.tsv"), emit: all

    script:
        parts = group.split('__')
        chr = parts[-1]
        exp = parts[-2] 
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
    

    publishDir  path: "${params.outdir}/Saige_eQTLS/${exp}/nominal",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(group),path(qtl_q_results)
    output:
        tuple val(group),path("${chr}__${exp}__all_vars_genes.tsv.gz"), emit: all  optional true
        tuple val(group),path("p_vals_lambd*"), emit: lambda optional true
        tuple val(group),path("*gene_lambdas.tsv"), emit: lambda_vals optional true
    script:
        parts = group.split('__')
        chr = parts[-1]
        exp = parts[-3]

        """
            prepend_gene_large.py --pattern 'cis_*' --column 'gene' --outfile '${chr}__${exp}__all_vars_genes.tsv'
            gzip ${chr}__${exp}__all_vars_genes.tsv
            mv cell_lambdas.tsv ${chr}__${exp}__gene_lambdas.tsv || echo 'not existant'
        """
}

process PHENOTYPE_PCs{
    // label 'process_medium'
    tag { sanitized_columns }

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }   

    publishDir  path: "${params.outdir}/Saige_eQTLS/PHENOTYPE_PCs/${sanitized_columns}",
        mode: "${params.copy_mode}",
        overwrite: "true"
    
    memory { 
        sizeInGB = saige_filt_expr_input.size() / 1e9 * 2 * task.attempt
        return (sizeInGB ).toString() + 'GB' 
    }   

    input:
        tuple val(sanitized_columns), path(saige_filt_expr_input),path(covariates)
        val(phenotype_pcs)

    output:
        tuple val(sanitized_columns), path("${sanitized_columns}_with_pheno_pcs.tsv"),path("covariates_new.txt"), emit: output_pheno optional true

    script:

    """
        export LD_LIBRARY_PATH=/opt/libstdc++-old:/usr/lib:/usr/local/lib:\$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH
        saige_phenotype_pcs_and_other_covs.py ${saige_filt_expr_input} ${sanitized_columns}_with_pheno_pcs.tsv ${phenotype_pcs} ${covariates}
    """
}

process H5AD_TO_SAIGE_FORMAT {
    label 'process_medium'
    tag { sanitized_columns }

    // Specify the number of forks (10k)
    // maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }   

    memory { 
        sizeInGB = adata.size() / 1e9 * 1.25 * task.attempt
        return (sizeInGB ).toString() + 'GB' 
    }   
    publishDir  path: "${params.outdir}/Saige_eQTLS/PHENOTYPE_PCs/${sanitized_columns}",
        mode: "${params.copy_mode}",
        overwrite: "true"

    input:
        each path(h5ad)  
        path(bridge)  
        val(aggregation_columns)
        path(genotype_pcs)
        path(genome_annotation)

    output:
        tuple val(sanitized_columns), path("output_agg/*/*/saige_filt_expr_input.tsv"),path("output_agg/*/*/covariates.txt"),emit:output_pheno
        tuple val(sanitized_columns),path("output_agg/*/*/test_genes.txt"),emit:gene_chunk
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

    if ("${params.analysis_subentry}"==''){
        cond1 = " --condition_col 'NULL' --condition 'NULL' "
    }else{
        cond1 = " --condition_col '${aggregation_columns}' --condition '${params.analysis_subentry}' "
    }

    if (params.chromosomes_to_test) {
        chromosomes_as_string = params.chromosomes_to_test.join(',')
        cond2 = " --chr ${chromosomes_as_string} --genome ${genome_annotation}"
    } else {
        cond2 = ""
    }

    """
        echo ${sizeInGB}
        bridge='${bridge}'
        nperc=${params.percent_of_population_expressed}
        condition_col="${aggregation_columns}" #Specify 'NULL' if want to include all cells
        condition="${aggregation_columns}" #Specify 'NULL' if want to include all cells
        
        scale_covariates=${params.SAIGE.scale_covariates}
        expression_pca=${params.SAIGE.nr_expression_pcs}
        aggregate_on="${aggregation_columns}"

        mkdir output_agg
        prep_adata_saige.py \
            --phenotype__file ${h5ad} \
            --bridge \$bridge \
            --aggregate_on \$aggregate_on \
            --genotype_pc__file ${genotype_pcs} \
            --genotype_id '${params.gt_id_column}' \
            --sample_id '${params.sample_column}' \
            --general_file_dir ./output_agg \
            --covariates '${params.covariates.adata_obs_covariate}' \
            --nperc \$nperc \
            --gtf_gene_identifier ${params.gtf_gene_identifier} \
            --min ${params.n_min_cells} \
            --scale_covariates \$scale_covariates \
            --expression_pca \$expression_pca --cell_percentage_threshold ${params.cell_percentage_threshold} \
            ${cov_col} ${cond1} ${cond2}
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
    // \${chr1}___nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_gene_1
    script:
    """
        echo ${sanitized_columns}
        echo ${test_genes}
        split -l ${chunk_size} ${test_genes} chunk_${sanitized_columns}_
    """
}

process DETERMINE_TSS_AND_TEST_REGIONS {
    label 'process_medium'

    // Specify the number of forks (10k)
    maxForks 200

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
        prepeare_Saige_regions_for_cis.py --annotation_file ${annotation_file} --genes ${genes_list} --gtf_type ${params.gtf_type} --window ${params.windowSize} --gtf_gene_identifier ${params.gtf_gene_identifier}
    """
}

workflow SAIGE_qtls{
    take:
        genotype_pcs
        phenotype_file
        plink_path
        plink_path_bed
        genotypes_saige
        genome_annotation
        genotype_phenotype_mapping_file

    main:
        log.info('------- Running SAIGE QTLs ------- ')
        pheno =  Channel.of()
        gene =  Channel.of()
        // Check if analysis_subentry is provided
        if (params.analysis_subentry != '') {
            log.info("------- Analysing ${params.analysis_subentry} celltypes ------- ")
            // Split the analysis_subentry parameter into a list of patterns
            valid_files = phenotype_file.filter { file ->
                params.analysis_subentry.split(',').any { pattern -> "${file}".contains("__${pattern}__") }
                
            }
        } else {
            log.info('------- Analysing all celltypes ------- ')
            valid_files = phenotype_file
        }

        if (params.existing_phenotype_pcs && file(params.existing_phenotype_pcs).exists()) {
            Channel
                .fromPath("${params.existing_phenotype_pcs}/*/*_with_pheno_pcs.tsv", checkIfExists: true)
                .map { pcs_file ->
                    def parent_dir = pcs_file.parent.name
                    def cov_file = pcs_file.parent.resolve("covariates_new.txt")
                    tuple(parent_dir, pcs_file, cov_file)
                }
                .set { pheno }

            Channel
                .fromPath("${params.existing_phenotype_pcs}/*/test_genes.txt", checkIfExists: true)
                .map { genes_file ->
                    def parent_dir = genes_file.parent.name
                    tuple(parent_dir, genes_file)
                }
                .set { gene }
        }else{

            H5AD_TO_SAIGE_FORMAT(
                valid_files,
                genotype_phenotype_mapping_file,
                params.aggregation_columns,
                genotype_pcs,
                genome_annotation
            )
            pheno = H5AD_TO_SAIGE_FORMAT.out.output_pheno
            gene = H5AD_TO_SAIGE_FORMAT.out.gene_chunk
            pheno.subscribe { println "pheno: $it" }
        
            PHENOTYPE_PCs(pheno,params.SAIGE.nr_expression_pcs)
            pheno = PHENOTYPE_PCs.out.output_pheno
        }


        CHUNK_GENES(gene,params.chunkSize)
        result = CHUNK_GENES.out.output_genes.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { [first, it] }
        }

        result.combine(pheno, by: 0).set{pheno_chunk}
 
        chromosomes_to_test = (params.chromosomes_to_test && params.chromosomes_to_test.size() > 0)
            ? Channel.of(params.chromosomes_to_test)
            : Channel.of((1..24).toList())   


        if (params.existing_sparse_grm){
            sparseGRM = Channel
                .fromPath(params.existing_sparse_grm + "/sparseGRM_*.mtx")
                .ifEmpty { error " No sparseGRM .mtx file found in ${params.existing_sparse_grm}" }

            sparseGRM_sample = Channel
                .fromPath(params.existing_sparse_grm + "/sparseGRM_*.sampleIDs.txt")
                .ifEmpty { error " No sparseGRM sample ID file found in ${params.existing_sparse_grm}" }

        }else{
            CREATE_SPARSE_GRM(plink_path)
            sparseGRM = CREATE_SPARSE_GRM.out.sparseGRM
            sparseGRM_sample = CREATE_SPARSE_GRM.out.sparseGRM_sample
        }



        SAIGE_S1(pheno_chunk.combine(plink_path_bed),sparseGRM,sparseGRM_sample)

        DETERMINE_TSS_AND_TEST_REGIONS(SAIGE_S1.out.output,genome_annotation)
        for_cis_input = DETERMINE_TSS_AND_TEST_REGIONS.out.output_genes
        SAIGE_S2_CIS(for_cis_input.combine(genotypes_saige),sparseGRM,sparseGRM_sample)
        output_s2 = SAIGE_S2_CIS.out.output
        // output_s2.subscribe { println "output_s2: $it" }
        
        agg_output = SAIGE_S2_CIS.out.for_aggregation

        // HERE WE either run the cis or trans qtl mapping. For cis we loop through each of the chunks whereas in trans we can run all together.
        SAIGE_QVAL_COR(output_s2)
        SAIGE_S3(SAIGE_QVAL_COR.out.output)

        // ########## Collecting Chunk outputs.  ###############
        SAIGE_S2_for_aggregation = agg_output.flatMap { item ->
            def (first, second) = item
            if (!(second instanceof Collection)) {
                second = [second] // Wrap single value in a list
            }
            return second.collect { 
                def fileName = "${it}".split('/').last() // Extract the filename by splitting on '/'
                def number = fileName.split('__')[0] // Extract the number before the first '__'
                ["${first}__${number}", it] 
            }
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
        AGGREGATE_QTL_RESULTS(SAIGE_S3_for_aggregation.groupTuple(by: 0))
        AGGREGATE_QTL_ALLVARS(SAIGE_S2_for_aggregation.groupTuple(by: 0))
        AGGREGATE_ACAT_RESULTS(SAIGE_S3_for_aggregation_ACAT.groupTuple(by: 0))
        // CONDITIONAL_QTL(SAIGE_QVAL_COR.out.for_conditioning)

}