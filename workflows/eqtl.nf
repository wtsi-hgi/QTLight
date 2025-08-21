

include {PREPROCESS_GENOTYPES} from '../modules/local/preprocess_genotypes/main' 
include {PLINK_CONVERT;PLINK_CONVERT__GRM;PGEN_CONVERT;BGEN_CONVERT; PGEN_TO_BED_CONVERT_FOR_QTLS; PGEN_TO_BED_CONVERT_FOR_GRM} from '../modules/local/plink_convert/main' 
include {SUBSET_GENOTYPE} from '../modules/local/subset_genotype/main' 
include {GENOTYPE_PC_CALCULATION} from '../modules/local/genotype_pc_calculation/main' 
include {SPLIT_PHENOTYPE_DATA} from '../modules/local/split_phenotype_data/main' 
include {NORMALISE_and_PCA_PHENOTYPE} from '../modules/local/normalise_and_pca/main' 
include {LIMIX_eqtls} from '../modules/local/limix/main'
include {PREPROCESS_SAMPLE_MAPPING} from '../modules/local/preprocess_sample_mapping/main'
include {NORMALISE_ANNDATA; REMAP_GENOTPE_ID} from '../modules/local/normalise_anndata/main'
include {AGGREGATE_UMI_COUNTS; SPLIT_AGGREGATION_ADATA; ORGANISE_AGGREGATED_FILES} from '../modules/local/aggregate_UMI_counts/main'
include {PREPERE_EXP_BED;PREPERE_COVARIATES; PREP_SAIGE_COVS} from '../modules/local/prepere_exp_bed/main'
include {TENSORQTL_eqtls} from '../modules/local/tensorqtl/main'
include {JAXQTL_eqtls} from '../modules/local/jaxqtl/main'
include {SAIGE_qtls} from '../modules/local/saige/main'
include {SUBSET_PCS; MERGE_COVARIATES} from '../modules/local/covar_processing/main'
include {KINSHIP_CALCULATION} from "$projectDir/modules/local/kinship_calculation/main"

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow EQTL {

    log.info 'Lets run eQTL mapping'
    // if single cell data then have to prepere pseudo bulk dataset.
    if (params.method=='bulk'){
        log.info '------ Bulk analysis ------'

        if (params.genotype_phenotype_mapping_file!=''){
            log.info '------ Genotype - Phenotype file not used as the phenotype file already contains matching genotype IDs ------'
            genotype_phenotype_mapping_file=params.genotype_phenotype_mapping_file
            phenotype_file=params.phenotype_file
            input_channel = Channel.fromPath(genotype_phenotype_mapping_file)
            input_channel.splitCsv(header: true, sep: params.input_tables_column_delimiter)
                .map{row->tuple(row.Genotype)}.distinct()
                .set{channel_input_data_table}
            channel_input_data_table = channel_input_data_table.collect().flatten().distinct()
            input_channel.splitCsv(header: true, sep: params.input_tables_column_delimiter)
                .map{row->row.Sample_Category}.distinct().set{condition_channel}
            SPLIT_PHENOTYPE_DATA(genotype_phenotype_mapping_file,phenotype_file,condition_channel)
            phenotype_condition = SPLIT_PHENOTYPE_DATA.out.phenotye_file
        }else{
            phenotype_condition = Channel.from("foo").map { foo -> tuple("full",file(params.phenotype_file),file("$projectDir/assets/fake_file.fq")) }
        }

    }else if (params.method=='single_cell'){
        log.info '------ Scrna analysis ------'
        if (params.pre_aggregated_counts_folder==''){    
            if (params.normalise_before_or_after_aggregation=='before'){
                log.info '------ We are normalising the entire andata together prior spliting by condition defined by  params.aggregation_columns------'
                NORMALISE_ANNDATA(params.phenotype_file)
                pheno = NORMALISE_ANNDATA.out.adata
            }else{
                pheno = params.phenotype_file
            }
            
            if (params.split_aggregation_adata){
                if (params.existing_split_adata_dir && file(params.existing_split_adata_dir).exists()) {
                    log.info "Using user-provided split AnnData files from: ${params.existing_split_adata_dir}"
                    // Get all h5ad files from the provided directory
                    adata = Channel.fromPath("${params.existing_split_adata_dir}/*.h5ad")
                } else {
                    log.info "No precomputed split files provided, running SPLIT_AGGREGATION_ADATA"
                    // Run the module that generates splits
                    SPLIT_AGGREGATION_ADATA(pheno, params.aggregation_columns)
                    // Flatten outputs
                    adata = SPLIT_AGGREGATION_ADATA.out.split_phenotypes.flatten()
                }
            }else{
                adata = pheno
            }   
            
            if (params.normalise_before_or_after_aggregation=='after'){
                // here we normalise the adata per splits
                log.info '------ We are normalising the split andata defined by params.aggregation_columns------'
                NORMALISE_ANNDATA(adata)
                splits_h5ad = NORMALISE_ANNDATA.out.adata
            }else{
                splits_h5ad = adata
            }
        
            if (!params.SAIGE.run || params.TensorQTL.run || params.LIMIX.run || params.JAXQTL.run) {
                log.info "performing aggregation"
                if (params.analysis_subentry != '') {
                    log.info("------- Analysing ${params.analysis_subentry} celltypes ------- ")
                    // Split the analysis_subentry parameter into a list of patterns
                    splits_h5ad2 = splits_h5ad
                        .filter { file -> 
                            def matches = params.analysis_subentry.split(',').any { pattern -> "${file}".contains("__${pattern}__") }
                            if (matches) {
                                println "MATCH: File=${file}"
                            } 
                            return matches
                        }
                } else {
                    log.info('------- Analysing all celltypes ------- ')
                    valid_files =  umi_counts_phenotype_genotype_file
                }

                AGGREGATE_UMI_COUNTS(splits_h5ad2, params.aggregation_columns, params.gt_id_column, params.sample_column, params.n_min_cells, params.n_min_individ)
                covariates_by_name = AGGREGATE_UMI_COUNTS.out.sample_covariates.map { file ->
                    def fname = file.getBaseName().replaceAll(/___sample_covariates/, '')
                    return tuple(fname, file)
                }
                phenotype_genotype_file = AGGREGATE_UMI_COUNTS.out.phenotype_genotype_file
                genotype_phenotype_mappings = AGGREGATE_UMI_COUNTS.out.genotype_phenotype_mapping.flatten()
            } else {
                log.info "Skipping AGGREGATE_UMI_COUNTS — only SAIGE is enabled"
                phenotype_genotype_file = Channel.empty()
                genotype_phenotype_mappings = Channel.empty()
                covariates_by_name = Channel.empty()
            }
        }else{
            log.info "Looking for existing files '___phenotype_file.tsv' and '___genotype_phenotype_mapping.tsv' in ${params.pre_aggregated_counts_folder}/*/*phenotype_file.tsv"
            Channel
                .fromPath(params.pre_aggregated_counts_folder+'/*/*___phenotype_file.tsv').ifEmpty { error "No files found in data/ directory" }
                .map{file1 ->
                    def parts = "${file1}".split('___')
                    def name_pre = parts[ parts.size() - 2 ]
                    def name_parts = "${name_pre}".split('/')
                    def name = name_parts[ name_parts.size() - 1 ]
                    def replacedFullPath = "${file1}".replace('___phenotype_file.tsv', '___genotype_phenotype_mapping.tsv')
        
                    tuple( name, file(file1), file(replacedFullPath) )  }
                .set { phenotype_genotype_file }
            
            Channel
                .fromPath(params.pre_aggregated_counts_folder + '/*/*___sample_covariates.tsv')
                .ifEmpty { error "No sample_covariates files found in ${params.pre_aggregated_counts_folder}" }
                .map { file1 ->
                    def parts = "${file1}".split('___')
                    def name_pre = parts[ parts.size() - 2 ]
                    def name_parts = "${name_pre}".split('/')
                    def name = name_parts[ name_parts.size() - 1 ]
                    tuple(name, file(file1))
                }
                .set { covariates_by_name }
                
            Channel
                .fromPath(params.pre_aggregated_counts_folder+'/*/*___phenotype_file.tsv').ifEmpty { error "No FASTQ files found in data/ directory" }
                .map{file1 ->
                    def replacedFullPath = "${file1}".replace('___phenotype_file.tsv', '___genotype_phenotype_mapping.tsv')
                    file(replacedFullPath)   }
                .set { genotype_phenotype_mappings }
        }

        if (params.split_aggregation_adata==false){
            ORGANISE_AGGREGATED_FILES(phenotype_genotype_file)
            ORGANISE_AGGREGATED_FILES.out.phenotype_files_tsv.splitCsv(header: true, sep: params.input_tables_column_delimiter)
                .map{row->tuple(row.name, file(row.phen_file), file(row.gp_fi) )}
            .set{umi_counts_phenotype_genotype_file}
        }else{
            umi_counts_phenotype_genotype_file = phenotype_genotype_file
        }

        if (params.analysis_subentry != '') {
            log.info("------- Analysing ${params.analysis_subentry} celltypes ------- ")
            // Split the analysis_subentry parameter into a list of patterns
            valid_files = umi_counts_phenotype_genotype_file
                .filter { tuple -> 
                    def (sample, file, gp_mapping) = tuple
                    def matches = params.analysis_subentry.split(',').any { pattern -> "${file}".contains("__${pattern}__") }

                    if (matches) {
                        println "MATCH: Sample=${sample}, File=${file}, GO Mapping=${gp_mapping}"
                    } 
                    return matches
                }
        } else {
            log.info('------- Analysing all celltypes ------- ')
            valid_files =  umi_counts_phenotype_genotype_file
        }

        out2 = valid_files.map { data ->
                def (item, list1, list2) = data
                // Ensure list1 and list2 are processed as lists
                list1 = (list1 instanceof List) ? list1 : [list1]
                list2 = (list2 instanceof List) ? list2 : [list2]
                def result = []
                list1.eachWithIndex { val, idx ->
                    result << [val.toString().split('___')[0].split('/')[-1], val, list2[idx]]
                }
                return result
            }.flatMap { it }

        if(params.genotype_phenotype_mapping_file!=''){
            // Here user has provided a genotype phenotype file where the provided gt_id_column is contaiming a mapping file instead of actual genotype
            log.info("------- Genotype-Phenotype mapping file provided, will remap using this file - ${params.genotype_phenotype_mapping_file} ------- ")
            mapping_ch = Channel.fromPath(params.genotype_phenotype_mapping_file)
            REMAP_GENOTPE_ID(out2.combine(mapping_ch))
            phenotype_condition = REMAP_GENOTPE_ID.out.remap_genotype_phenotype_mapping
            genotype_phenotype_mapping_file = REMAP_GENOTPE_ID.out.genotype_phenotype_mapping
        }else{
            log.info("------- Genotype-Phenotype mapping file NOT provided, will asume that the ids are already correct ------- ")
            phenotype_condition = out2
            genotype_phenotype_mapping_file = genotype_phenotype_mappings.flatten()
        }

        genotype_phenotype_mapping_file.splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.Genotype)}.distinct()
            .set{channel_input_data_table2}
        channel_input_data_table=channel_input_data_table2.collect()
    }


    if (params.input_vcf){
        // VCF file processing
        log.info "---- VCF file provided - ${input_vcf} Lets preprocess genotypes and aply any bcftools filters ---"
        donorsvcf = Channel.from(params.input_vcf)
        if (params.genotypes.subset_genotypes_to_available){
            // Subset genotypes to available in expression data
            log.info "---- You have specified  ${params.genotypes.subset_genotypes_to_available} - will subset down the vcf file to only the donors available in the expression data---"
        
            SUBSET_GENOTYPE(donorsvcf,channel_input_data_table.collect())
            subset_genotypes = SUBSET_GENOTYPE.out.samplename_subsetvcf
        }else{
            subset_genotypes = donorsvcf
        }

        if (params.genotypes.apply_bcftools_filters){
            // preprocess vcf files to be in the right format for plink
            log.info "---- You have specified  ${params.genotypes.apply_bcftools_filters} - will apply these filters to vcf file before converting to other formats---"
            PREPROCESS_GENOTYPES(subset_genotypes)
            plink_convert_input=PREPROCESS_GENOTYPES.out.filtered_vcf
        }else{
            plink_convert_input=subset_genotypes  
        }

    }else{
        if (params.genotypes.preprocessed_pgen_file==''){
            plink_convert_input=Channel.of()
        }else{
            log.info "---- pgen file provided, will prioritise this over VCF---"
            plink_convert_input=Channel.from(params.genotypes.preprocessed_pgen_file)
        }
    }

    // This block is necesary for Saige as it needs BED for GRM construction and if we dont use dosages we also need this for the QTL matrix construction
    if (params.SAIGE.run || (params.genotypes.use_gt_dosage == false) || params.JAXQTL.run) {
        if (params.genotypes.preprocessed_bed_file==''){
            if (params.input_vcf){
                // BED file preparation
                log.info "---- Converting from vcf to plink to be used for GRMs and association testings if dosage is switched off ---"
                PLINK_CONVERT(plink_convert_input)
                PLINK_CONVERT__GRM(plink_convert_input)
                bim_bed_fam = PLINK_CONVERT.out.bim_bed_fam
                plink_path_bed = PLINK_CONVERT.out.plink_path
                bim_bed_fam__GRM = PLINK_CONVERT__GRM.out.bim_bed_fam
                plink_path_bed__GRM = PLINK_CONVERT__GRM.out.plink_path

            }else if (params.genotypes.preprocessed_pgen_file != '') {
                // BED file preparation from preprocessed PGEN
                log.info "--- PGEN provided and BED needed for GRMs — converting with PGEN_TO_BED_CONVERT ---"
                plink_path_pgen = Channel.from(params.genotypes.preprocessed_pgen_file)

                // Here we are creating two bed files - 1 that goes to QTL analysis for tools that only support BED formated genotypes and the other goes to GRM construction
                PGEN_TO_BED_CONVERT_FOR_QTLS(plink_path_pgen)
                bim_bed_fam    = PGEN_TO_BED_CONVERT_FOR_QTLS.out.bim_bed_fam
                plink_path_bed = PGEN_TO_BED_CONVERT_FOR_QTLS.out.plink_path_bed

                PGEN_TO_BED_CONVERT_FOR_GRM(plink_path_pgen)
                bim_bed_fam__GRM    = PGEN_TO_BED_CONVERT_FOR_GRM.out.bim_bed_fam
                plink_path_bed__GRM = PGEN_TO_BED_CONVERT_FOR_GRM.out.plink_path_bed
            }

        }else{
            log.info "--- BED already provided (WARNING - same BED will be Used for both associations and GRM) ---"
            plink_path_bed = Channel.from(params.genotypes.preprocessed_bed_file)
            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.bed", followLinks: true)
                .set { bed_files }

            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.bim", followLinks: true)
                .set { bim_files }

            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.fam", followLinks: true)
                .set { fam_files }
            bim_bed_fam = bim_files
                        .combine(bed_files)
                        .combine(fam_files) 
            bim_bed_fam__GRM    = bim_bed_fam
            plink_path_bed__GRM = plink_path_bed 
        }
    }


    if (params.genotypes.preprocessed_pgen_file==''){
        // USE VCF FILE
        log.info "--- preprocessed_pgen_file does not existm will create it ---"
        PGEN_CONVERT(plink_convert_input)
        plink_path_pgen = PGEN_CONVERT.out.plink_path
    }else{
        log.info "--- Using pgen file provided ${params.genotypes.preprocessed_pgen_file} ---"
        plink_path_pgen = Channel.from(params.genotypes.preprocessed_pgen_file)
    }

    // If use dosages we convert vcf to pgen
    // Otherwise we convert it to bed
    if(params.genotypes.use_gt_dosage){
        log.info "---We are using dosages, so pgen file will be used ---"
        plink_path = plink_path_pgen
    }else{
        log.info "---We are NOT using dosages, so plink1 file will be used ---"
        plink_path = plink_path_bed
    }

    // // GENOTYPE PCs
    if (params.covariates.genotype_pcs_file==''){
        log.info "--- Calculating genotype PCs ---"
        GENOTYPE_PC_CALCULATION(plink_path)
        genotype_pcs_file = GENOTYPE_PC_CALCULATION.out.gtpca_plink
    }else{
        log.info "--- Genotype PCs file already provided, will use this ---"
        genotype_pcs_file = Channel.fromPath(params.covariates.genotype_pcs_file)
    }

    // 4) Phenotype file preperation including PCs, normalisation
    genome_annotation = Channel.from(params.annotation_file)

    // Potentially add:
    // MBV method from QTLTools (PMID 28186259)  
    // RASCAL
    
    log.info "--- Normalising agregated penotype file  using: Filtering method: ${params.filter_method} Norm method: ${params.method} Inverse Transform: ${params.inverse_normal_transform} Norm method:  ${params.norm_method} Percent of Population expressed: ${params.percent_of_population_expressed} Sample PCa: ${params.use_sample_pca}---"
    NORMALISE_and_PCA_PHENOTYPE(phenotype_condition)
    log.info "--- testing these PCs: ${params.covariates.nr_phenotype_pcs} ---"
    Channel.of(params.covariates.nr_phenotype_pcs).splitCsv().flatten().set{pcs}
    NORMALISE_and_PCA_PHENOTYPE.out.for_bed.combine(pcs).set{test123}
         
    SUBSET_PCS(test123)

    if (params.covariates.adata_obs_covariate){
        log.info "--- Extracting aditional covariates from h5ad file: ${params.covariates.adata_obs_covariate} ---"
    
        test123_fixed = SUBSET_PCS.out.for_bed.map { row ->
            def pheno_full = row[0]
            def pheno_core = pheno_full.contains('__') ? pheno_full.split('__')[0..-2].join('__') : pheno_full
            tuple(pheno_core, row[0],row[1],row[2],row[3])  // [matching_key, full_row]
        }
        // test123_fixed.subscribe { println "test123_fixed: $it" }
       
        test123_fixed.combine(covariates_by_name,by:0).set{for_covs_merge}
        
        MERGE_COVARIATES(for_covs_merge)
        for_bed = MERGE_COVARIATES.out.for_bed_covs
    }else{
        for_bed = SUBSET_PCS.out.for_bed
    }


    if(params.genotypes.use_gt_dosage && (params.LIMIX.run || params.SAIGE.run )){
        if (params.genotypes.preprocessed_bgen_file==''){
            log.info "--- generating bgen file for Limix or SaigeQtl ---"
    
            BGEN_CONVERT(plink_convert_input)
            plink_path_limix = BGEN_CONVERT.out.plink_path
            genotypes_saige =  BGEN_CONVERT.out.plink_path
        }else{
            plink_path_limix = Channel.from(params.genotypes.preprocessed_bgen_file)
            genotypes_saige = Channel.from(params.genotypes.preprocessed_bgen_file)
        }
    }else if (params.LIMIX.run || params.SAIGE.run){
        plink_path_limix = plink_path_bed
        genotypes_saige = plink_path_bed
    }

    // LIMIX QTL mapping method
    if (params.LIMIX.run){
        log.info "--- Calculating Kinship matrix for Limix ---"
        KINSHIP_CALCULATION(plink_path)
        kinship_file = KINSHIP_CALCULATION.out.kinship_matrix

        filtered_pheno_channel = for_bed.map { tuple ->  
            [tuple[3], [tuple[0], tuple[1], tuple[2]]]
        }.flatten().collate(4)

        LIMIX_eqtls(
            filtered_pheno_channel,
            plink_path_limix,
            genome_annotation,
            kinship_file
        )
    }

    for_bed_channel = for_bed.map { tuple ->  [tuple[3],[[tuple[0],tuple[1],tuple[2]]]]}.flatten().collate(4)
    
    log.info "---Prepearing covariates and Bed files for QTL analysis for JaxQTL or TensorQtl---"
    PREPERE_COVARIATES(for_bed_channel,genotype_pcs_file)
    covs = PREPERE_COVARIATES.out.exp_bed

    PREPERE_EXP_BED(for_bed_channel,genome_annotation)
    beds = PREPERE_EXP_BED.out.exp_bed
    beds.combine(covs,by:0).set{tensorqtl_input}

    tensorqtl_input_ch_cleaned = tensorqtl_input.map { items ->
        def new_first = items[0].replaceFirst(/__[^_]+pcs\.tsv$/, '')
        [new_first, *items[1..-1]]
    }

    if (params.TensorQTL.run){
        TENSORQTL_eqtls(
            tensorqtl_input_ch_cleaned,
            plink_path,
        )
    }

    if (params.JAXQTL.run){
        JAXQTL_eqtls(
            tensorqtl_input,
            plink_path_bed,
        )
    }

    // SAIGE SCRNA QTL mapping method
    if (params.method=='single_cell'){
        if (params.SAIGE.run){
            log.info "---Running saigeQtls---"
    
            if (params.genotype_phenotype_mapping_file==''){
                saige_genotype_phenotype_mapping_file = Channel.of("$projectDir/assets/fake_file.fq")
            }else{
                saige_genotype_phenotype_mapping_file = Channel.fromPath(params.genotype_phenotype_mapping_file, followLinks: true, checkIfExists: true)
            }

            if (params.covariates.extra_covariates_file==''){
                extra_covariates_file = Channel.of("$projectDir/assets/fake_file.fq")
            }else{
                extra_covariates_file = Channel.fromPath(params.covariates.extra_covariates_file, followLinks: true, checkIfExists: true)
            }
            covs_Saige = PREP_SAIGE_COVS(genotype_pcs_file,extra_covariates_file)

            SAIGE_qtls(covs_Saige,adata,plink_path_bed__GRM,plink_path_bed,genotypes_saige,genome_annotation,saige_genotype_phenotype_mapping_file)
        }
    }

    // Generate plots of comparisons of eQTLs detected by all run methods.


}

