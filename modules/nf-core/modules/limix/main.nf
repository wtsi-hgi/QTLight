

process COLLECT_RESULTS{
     tag { condition }
    
    input:
        val(full_output_list)
        // tuple(val(chunking_range),val(condition),path(phenotypeFile),path(covariateFile),path(annotationFile))
        tuple(val(condition),path(phenotypeFile),path(covariateFile))

    output:
        path("${condition}"), emit: condition_all_qtls
    script:
        full_output_list=full_output_list.join("\n")
        """

        echo -e "${full_output_list}">out.txt
        link_files.py --files_input out.txt --condition ${condition}
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

process TEST{
    tag { condition }
    input:
        each chunking_range
        tuple(val(condition),path(phenotypeFile),path(covariateFile))
        each path(genotypeFile)
        each path(annotationFile)

    output:
        path("out.txt"), emit: filtered_vcf
    script:
        
        """
        echo "${chunking_range}">out.txt
        multiCorrect.R ${}
        """
}


process MULTIPLE_TESTING_CORRECTION{
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }    
        
    input:
        path(limix_qtl_path)


    output:
        path("${limix_qtl_path}/qtl_results_all_FDR*"), emit: filtered_vcf
    script:
        
        """
            multiCorrect.R ${limix_qtl_path}
        """
}


process LIMIX{
    
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    tag { "${condition} ${chunking_range2}" }
    maxForks 100
    //cache false        // cache results from run
    // scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.limix_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(chunking_range),val(condition),path(phenotypeFile),path(covariateFile),path(annotationFile),path(genotypeFile),path(kinship_path),path(individual2sample_filename))
        
    output:
        tuple(path("${condition}_snp_metadata*"),path("${condition}_qtl_results_*"),path("${condition}_feature_metadata_*") , emit: qtl_data)
    script:
        numberOfPermutations = params.numberOfPermutations
        minorAlleleFrequency = params.maf
        hwe = params.hwe
        callRate = 0.95
        windowSize = params.windowSize 
        blockSize = 1500

        outputFolder='./'
        chunking_range="${chunking_range}"
        chunking_range=chunking_range.replaceAll('\\[', "")
        chunking_range=chunking_range.replaceAll('\\]', "")
        chunking_range2 = chunking_range.replaceAll(':', "_").replaceAll('-', "_")
        """
            export NUMBA_CACHE_DIR=/tmp
            export MPLCONFIGDIR=/tmp 
            run_limix_QTL_analysis.py --plink ${genotypeFile}/plink_genotypes -af ${annotationFile} -pf ${phenotypeFile} -cf ${covariateFile} -od ${outputFolder} -smf ${individual2sample_filename} -rf ${kinship_path} -gr ${chunking_range} -np ${numberOfPermutations} -maf ${minorAlleleFrequency} -hwe ${hwe} -cr ${callRate} -c -gm standardize -w ${windowSize} --block_size ${blockSize}
            ln -s snp_metadata_${chunking_range2}.txt ${condition}_snp_metadata_${chunking_range2}.txt
            ln -s qtl_results_${chunking_range2}.h5 ${condition}_qtl_results_${chunking_range2}.h5
            ln -s feature_metadata_${chunking_range2}.txt ${condition}_feature_metadata_${chunking_range2}.txt
        """
    
}

workflow LIMIX_eqtls{
    take:
        limix_condition_chunking
        plink_genotype
        kinship_file
        genotype_phenotype_mapping_file
        condition

    main:

        limix_condition_chunking.splitCsv(header: true, sep: "\t").map{row ->  tuple(row.Range,row.condition,row.phenotypeFile,row.covars,row.anotation_file)}.set{chunking_channel}
        chunking_channel =chunking_channel.combine(plink_genotype)
        chunking_channel =chunking_channel.combine(kinship_file)
        chunking_channel =chunking_channel.combine(genotype_phenotype_mapping_file)

        LIMIX(chunking_channel)

        COLLECT_RESULTS(LIMIX.out.qtl_data.collect(),condition)

        AGGREGATE_QTL_RESULTS(COLLECT_RESULTS.out.condition_all_qtls)
        
        MULTIPLE_TESTING_CORRECTION(AGGREGATE_QTL_RESULTS.out.limix_qtl_path)

}