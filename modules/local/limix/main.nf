
process CHUNK_GENOME{
    tag "$condition, $nr_phenotype_pcs"
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
        
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        each path(genome_annotation)
        tuple path(phenotype_pcs),val(condition),path(mapping_file),path(phenotype_file)
        val(chunkSize)

    output:
        tuple val(condition),path(phenotype_file), path(phenotype_pcs),path("Chunging_file*.tsv"),path(mapping_file) , emit: filtered_chunking_file optional true
    script:
        nr_phenotype_pcs = phenotype_pcs.getSimpleName()
        if ("${params.chromosomes_to_test}"!=''){
            chromosomes_as_string = params.chromosomes_to_test.join(',')
            cond2 = " --chr ${chromosomes_as_string}"
        }else{
            cond2 = " "
        }
        """
            generate_chunking_file.py --genome_annotation ${genome_annotation} --chunk_size ${chunkSize} --phenotype_file ${phenotype_file} --covar_file ${phenotype_pcs} --condition ${condition}  --genotype_phenotype_file ${mapping_file} ${cond2}
        """
    
}

process COLLECT_RESULTS{
    tag { condition }
    label 'process_tiny'
    input:
        tuple val(condition), path(metadata),path(full_output_list),path(feature_metadata)
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
    

    publishDir  path: "${params.outdir}/Limix_eQTLS/${group1}__${group0}/${group2}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple val(condition), path(metadata),path(full_output_list),path(feature_metadata)
        
        
    output:
        path("results_${group1}"), emit: limix_qtl_path
        // path("${all_qtl_results}_all_mpr/qtl_results_all.txt"), emit: qtl_results_all_all_mpr
        // path("${all_qtl_results}_all_mpr/top_qtl_results_all.txt"), emit: top_qtl_results_all_all_mpr
        // path("${all_qtl_results}_qtls")
    script:
        all_qtl_results = "${condition}".replaceFirst(/.*\.tsv$/, '')

        // matcher = ("${condition}" =~ /.*?__(.*?)__(.*?)\.tsv$/)
        group0 = "${condition}".split('__')[0]
        group1 =  "${condition}".split('__')[1]
        // group1 = matcher[0][1] // Extracts 'Mono_all'
        group2 = "${condition}".split('__')[2].replaceFirst(/\.tsv/, '')

        """
            mkdir results_${group1}
            echo "${condition}"
            minimal_postprocess.py -id ./ -od results_${group1} -sfo -tfb 
            minimal_postprocess.py -id ./ -od results_${group1} -sfo -mrp 0.05 
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

    tag { "${condition} ${annotationFile}" }
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.limix_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(phenotypeFile),path(covariateFile),path(annotationFile),path(individual2sample_filename),path(genotypeFile),path(kinship_path))
        
    output:
        tuple val("${condition}__${covariateFile}"), path("${condition}_snp_metadata*"),path("${condition}_qtl_results_*"),path("${condition}_feature_metadata_*") , emit: qtl_data optional true
    script:
        numberOfPermutations = params.LIMIX.numberOfPermutations
        minorAlleleFrequency = params.maf
        hwe = params.LIMIX.hwe
        callRate = params.LIMIX.callRate
        windowSize = params.windowSize 
        blockSize = params.LIMIX.blockSize

        outputFolder='./'
        chunk_number = "${annotationFile}".replaceFirst(/.*__(\d+)\.tsv$/, '$1')

        if (params.genotypes.use_gt_dosage) {
            genotypeFile2 = "--bgen ${genotypeFile}"
        }else{
            genotypeFile2 = "--plink ${genotypeFile}/plink_genotypes "
        }


        """
            export NUMBA_CACHE_DIR=\$PWD
            export MPLCONFIGDIR=\$PWD
            export HOME=\$PWD
            cut -f 1,2 -d \$'\\t' ${individual2sample_filename} > data_no_sample_category.txt
            run_limix_QTL_analysis.py ${genotypeFile2} -af ${annotationFile} -pf ${phenotypeFile} -cf ${covariateFile} -od ${outputFolder} -smf data_no_sample_category.txt -rf ${kinship_path} -np ${numberOfPermutations} -maf ${minorAlleleFrequency} -hwe ${hwe} -cr ${callRate} -c -gm standardize -w ${windowSize} --block_size ${blockSize}
            mv snp_metadata_all.txt ${condition}_snp_metadata_${chunk_number}.txt || echo 'not available'
            mv qtl_results_all.h5 ${condition}_qtl_results_${chunk_number}.h5 || echo 'not available'
            mv feature_metadata_all.txt ${condition}_feature_metadata_${chunk_number}.txt || echo 'not available'
        """
    
}

workflow LIMIX_eqtls{
    take:
        filtered_pheno_channel
        plink_genotype
        genome_annotation
        kinship_file
    main:


        CHUNK_GENOME(genome_annotation,filtered_pheno_channel,params.chunkSize)

        chunking_channel=CHUNK_GENOME.out.filtered_chunking_file
        

        result = chunking_channel.flatMap { item ->
            def (condition,phenotype_file,phenotype_pcs,chunging_file,mapping_file) = item
            if (!(chunging_file instanceof Collection)) {
                chunging_file = [chunging_file] // Wrap single value in a list
            }
            return chunging_file.collect { [condition,phenotype_file,phenotype_pcs,it,mapping_file] }
        }
        // condition.splitCsv(header: true, sep: "\t").map{row ->  tuple(row.Range,row.condition,row.phenotypeFile,row.covars,row.anotation_file,row.genotype_phenotype_file)}.set{chunking_channel}
        chunking_channel = result.combine(plink_genotype)
        chunking_channel = chunking_channel.combine(kinship_file)
        // chunking_channel.subscribe { println "chunking_channeloutput_s2 dist: $it" }
        LIMIX(chunking_channel) // These are then passed to limix model which is either run in cis or trans mode.
        inp_ch2 = LIMIX.out.qtl_data.groupTuple(by: 0)
        // COLLECT_RESULTS(inp_ch2) // Results are then collected.

        AGGREGATE_QTL_RESULTS(inp_ch2) // QTL results are then aggregated.
        
        MULTIPLE_TESTING_CORRECTION(AGGREGATE_QTL_RESULTS.out.limix_qtl_path) // Multiple testing is performed.

}