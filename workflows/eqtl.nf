

include {PREPROCESS_GENOTYPES} from '../modules/local/preprocess_genotypes/main' 
include {PLINK_CONVERT;PLINK_CONVERT__GRM;PGEN_CONVERT;BGEN_CONVERT; PGEN_TO_BED_CONVERT_FOR_QTLS; PGEN_TO_BED_CONVERT_FOR_GRM} from '../modules/local/plink_convert/main' 
include {SUBSET_GENOTYPE} from '../modules/local/subset_genotype/main' 
include {GENOTYPE_PC_CALCULATION} from '../modules/local/genotype_pc_calculation/main' 
include {SPLIT_PHENOTYPE_DATA} from '../modules/local/split_phenotype_data/main' 
include {NORMALISE_and_PCA_PHENOTYPE} from '../modules/local/normalise_and_pca/main' 
include {LIMIX_eqtls} from '../modules/local/limix/main'
include {PREPROCESS_SAMPLE_MAPPING} from '../modules/local/preprocess_sample_mapping/main'
include {NORMALISE_ANNDATA; REMAP_GENOTPE_ID} from '../modules/local/normalise_anndata/main'
include {AGGREGATE_UMI_COUNTS; COMBINE_AGGREGATES; AGGREGATE_UMI_COUNTS as AGGREGATE_UMI_COUNTS__CHUNKS; SPLIT_AGGREGATION_ADATA; ORGANISE_AGGREGATED_FILES; CHUNK_ADATA_FOR_AGGREGATION} from '../modules/local/aggregate_UMI_counts/main'
include {PREPERE_EXP_BED;PREPERE_COVARIATES; PREP_SAIGE_COVS} from '../modules/local/prepere_exp_bed/main'
include {TENSORQTL_eqtls} from '../modules/local/tensorqtl/main'
include {JAXQTL_eqtls} from '../modules/local/jaxqtl/main'
include {SAIGE_qtls} from '../modules/local/saige/main'
include {QUASAR} from '../modules/local/quasar/main'
include {SUBSET_PCS; MERGE_COVARIATES} from '../modules/local/covar_processing/main'
include {KINSHIP_CALCULATION} from "$projectDir/modules/local/kinship_calculation/main"
include{CREATE_SPARSE_GRM} from '../modules/local/sparse_grm_creation/main'

/*
===============================
    QTLight eQTL Workflow
===============================
Overview:
- Bulk mode: expects counts (TSV) + optional mapping file.
- Single-cell mode: expects AnnData, performs aggregation and normalisation.
- Config-driven: see `conf/base.config` for default parameter descriptions.
- Supported tools: SAIGE, TensorQTL, LIMIX, JAXQTL.

Workflow stages:
1) Genotype preprocessing (VCF/PGEN → PGEN/BED/BGEN).
2) Phenotype aggregation and normalisation (bulk or scRNA).
3) Covariate and PCA handling.
4) Run QTL tools.
*/

workflow EQTL {

    log.info 'Lets run eQTL mapping'
    genome_annotation = Channel.from(params.annotation_file)
    // if single cell data then have to prepere pseudo bulk dataset.
    if (params.method=='bulk'){
        log.info '------ Bulk analysis ------'
        // --- Bulk workflow ---
        // Prepare phenotype data directly from counts and mappings.
        // Skip single-cell aggregation steps, but still normalise and add PCs.
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
        // --- Check for pre-aggregated counts ---
        // If `params.pre_aggregated_counts_folder` is empty, run the full aggregation
        // pipeline (normalisation + aggregation + mapping).  
        // If a folder path is provided, the pipeline will instead reuse existing
        // pre-aggregated phenotype files and mapping TSVs, avoiding recomputation
        // across runs (important for large scRNA-seq datasets).
        if (params.pre_aggregated_counts_folder==''){    
            if (params.normalise_before_or_after_aggregation=='before'){
                log.info '------ We are normalising the entire andata together prior spliting by condition defined by  params.aggregation_columns------'
                // If user wants normalization before aggregation, normalise whole AnnData.
                // Otherwise normalisation happens per split group. default: after
                NORMALISE_ANNDATA(params.phenotype_file)
                pheno = NORMALISE_ANNDATA.out.adata
            }else{
                pheno = params.phenotype_file
            }

            // --- Handle split AnnData by aggregation columns ---
            // Controlled by `params.split_aggregation_adata` and `params.existing_split_adata_dir`.
            //
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
                // If user wants normalization before aggregation, normalise whole AnnData.
                // Otherwise normalisation happens per split group. default: after
                log.info '------ We are normalising the split andata defined by params.aggregation_columns------'
                NORMALISE_ANNDATA(adata)
                splits_h5ad = NORMALISE_ANNDATA.out.adata
            }else{
                splits_h5ad = adata
            }
        

            // Purpose:
            //   Prepare pseudo-bulk expression matrices from single-cell AnnData, grouped by
            //   `aggregation_columns` (e.g. cell type, stimulation condition). These aggregated
            //   counts are required for TensorQTL, LIMIX, and JAXQTL, but not for SAIGE when
            //   run alone (SAIGE can work directly on cell-level data).
            // Behaviour:
            //   - If at least one of TensorQTL, LIMIX, or JAXQTL is enabled, or SAIGE is disabled,
            //     run `AGGREGATE_UMI_COUNTS` to generate:
            //       • Pseudo-bulk phenotype matrices per donor/condition
            //       • Genotype–phenotype mapping files
            //       • Sample-level covariates (e.g. cell counts per donor)
            //   - If `params.analysis_subentry` is set, only aggregate specified subsets
            //     (e.g. "Mono,B,DC"); otherwise aggregate all cell types.
            //   - If only SAIGE is enabled, skip aggregation and pass empty channels.

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
                    splits_h5ad2 =  splits_h5ad
                }

                if (params.chunk_adata_for_aggregation){
                    CHUNK_ADATA_FOR_AGGREGATION(splits_h5ad2,params.sample_column,10)
                    AGGREGATE_UMI_COUNTS__CHUNKS(CHUNK_ADATA_FOR_AGGREGATION.out.split_phenotypes.flatten(), params.aggregation_columns, params.gt_id_column, params.sample_column, params.n_min_cells, 1, 'chunk')
                    grouped_chunks = AGGREGATE_UMI_COUNTS__CHUNKS.out.phenotype_genotype_file_with_sample_covs.groupTuple(by:0)
                    COMBINE_AGGREGATES(grouped_chunks,params.n_min_individ)
                    sample_covariates = COMBINE_AGGREGATES.out.sample_covariates
                    phenotype_genotype_file = COMBINE_AGGREGATES.out.phenotype_genotype_file
                    genotype_phenotype_mapping = COMBINE_AGGREGATES.out.genotype_phenotype_mapping.flatten()
                }else{
                    // This may be very slow and mem demanding for large adata objects. Might want to do a two pass aggregation where smaller subsets of donors h5ad_files are produced.
                    AGGREGATE_UMI_COUNTS(splits_h5ad2, params.aggregation_columns, params.gt_id_column, params.sample_column, params.n_min_cells, params.n_min_individ, 'full')
                    sample_covariates = AGGREGATE_UMI_COUNTS.out.sample_covariates
                    phenotype_genotype_file = AGGREGATE_UMI_COUNTS.out.phenotype_genotype_file
                    genotype_phenotype_mapping = AGGREGATE_UMI_COUNTS.out.genotype_phenotype_mapping.flatten()
                }


                covariates_by_name = sample_covariates.map { file ->
                    def fname = file.getBaseName().replaceAll(/___sample_covariates/, '')
                    return tuple(fname, file)
                }

            } else {
                log.info "Skipping AGGREGATE_UMI_COUNTS — only SAIGE is enabled"
                phenotype_genotype_file = Channel.empty()
                genotype_phenotype_mappings = Channel.empty()
                covariates_by_name = Channel.empty()
            }
        }else{
            // --- Reuse pre-aggregated files ---
            // If `params.pre_aggregated_counts_folder` is provided, the pipeline expects that
            // pseudo-bulk expression and mapping files already exist on disk.  
            // Specifically, it looks for:
            //   • `___phenotype_file.tsv` (aggregated expression per donor/condition)
            //   • `___genotype_phenotype_mapping.tsv` (bridging file linking samples to genotypes)
            //   • `___sample_covariates.tsv` (optional covariates per donor/sample)
            // These are collected into channels for downstream analysis instead of re-running
            // the aggregation/normalisation modules.
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


        // --- Organise aggregated phenotype/mapping files ---
        // If `split_aggregation_adata = false`, run `ORGANISE_AGGREGATED_FILES` to combine
        // and index phenotype + genotype mapping TSVs into a single structured table.  
        // If splitting was already applied (`split_aggregation_adata = true`), the
        // pre-split channel can be used directly without extra organisation.
        if (params.split_aggregation_adata==false){
            ORGANISE_AGGREGATED_FILES(phenotype_genotype_file)
            ORGANISE_AGGREGATED_FILES.out.phenotype_files_tsv.splitCsv(header: true, sep: params.input_tables_column_delimiter)
                .map{row->tuple(row.name, file(row.phen_file), file(row.gp_fi) )}
            .set{umi_counts_phenotype_genotype_file}
        }else{
            umi_counts_phenotype_genotype_file = phenotype_genotype_file
        }


        // --- Select cell types / conditions to analyse ---
        // Controlled by `params.analysis_subentry`.
        //
        // Purpose:
        //   Ensure that only the requested cell types (e.g. "Mono,B,DC") are analysed.
        //   This filtering step is repeated here because if the workflow is reusing
        //   pre-aggregated files (`pre_aggregated_counts_folder`), the earlier
        //   subsetting step is skipped. Without this block, all cell types in the
        //   pre-aggregated folder would be included by default.
        //
        // Behaviour:
        //   - If `analysis_subentry` is set, filter `umi_counts_phenotype_genotype_file`
        //     so that only matching entries remain.
        //   - If not set, analyse all available cell types.
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


        // --- Handle genotype–phenotype mapping ---
        // Controlled by `params.genotype_phenotype_mapping_file`.
        //
        // Purpose:
        //   Align donor IDs between genotype data (VCF/PLINK) and phenotype data
        //   (aggregated expression). In some datasets, the AnnData `gt_id_column`
        //   already matches the genotype IDs directly; in others, a separate
        //   genotype–phenotype mapping (bridging) file is required.
        //
        // Behaviour:
        //   - If a mapping file is provided, run `REMAP_GENOTPE_ID` to replace the
        //     phenotype IDs with the correct genotype IDs based on the mapping.
        //   - If no mapping file is provided, assume IDs are already harmonised and
        //     use them as-is.
        //
        // Notes:
        //   - The mapping file should include at least: `phenotype_id`, `genotype_id`,
        //     and `Sample_Category`.
        //   - This step ensures consistency between aggregated phenotype data and
        //     genotype matrices used downstream.
        if(params.genotype_phenotype_mapping_file!=''){
            // Here user has provided a genotype phenotype file where the provided gt_id_column is contaiming a mapping file instead of actual genotype
            log.info("------- Genotype-Phenotype mapping file provided, will remap using this file - ${params.genotype_phenotype_mapping_file} ------- ")
            mapping_ch = Channel.fromPath(params.genotype_phenotype_mapping_file)
            REMAP_GENOTPE_ID(out2.combine(mapping_ch))
 
            REMAP_GENOTPE_ID.out.remap_genotype_phenotype_mapping
                .map { sc, pheno, remap ->
                    def f = (remap instanceof List) ? remap[0] : remap
                    def n = f.newReader().lines().count() as long
                    tuple(sc, pheno, f, n)
                }
                .filter { sc, pheno, f, n -> n > params.n_min_individ }
                .map    { sc, pheno, f, n -> tuple(sc, pheno, f) }
                .set { phenotype_condition_remap }

            REMAP_GENOTPE_ID.out.genotype_phenotype_mapping
                .map { f -> tuple(f, f.newReader().lines().count()) }
                .filter { f, n -> n > params.n_min_individ }
                .map { f, n -> f }   // keep only the file, drop the count
                .set { genotype_phenotype_mapping_file_remap }

        }else{
            log.info("------- Genotype-Phenotype mapping file NOT provided, will asume that the ids are already correct ------- ")
            phenotype_condition_remap = out2
            genotype_phenotype_mapping_file_remap = genotype_phenotype_mappings.flatten()
        }

        genotype_phenotype_mapping_file_remap.splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.Genotype)}.distinct()
            .set{channel_input_data_table2}
        channel_input_data_table=channel_input_data_table2.collect()
    }


    // --- Genotype preparation ---
    // Purpose:
    //   Convert and harmonise genotype input into the formats required by
    //   downstream QTL tools. Depending on input type (VCF/PGEN/BED) and
    //   parameter settings, this block generates one or more of:
    //     • PGEN (dosages)  → required if use_gt_dosage = true
    //     • BED (hard calls) → required for SAIGE (GRM) and JAXQTL
    //     • BGEN            → required for LIMIX or SAIGE (if dosages are used)
    //     • Genotype PCs    → always calculated unless precomputed file provided
    //
    // Controlled by:
    //   params.input_vcf, params.genotypes.* (subset_genotypes_to_available,
    //   apply_bcftools_filters, use_gt_dosage, preprocessed_*_file), and
    //   which QTL methods are enabled (SAIGE, LIMIX, JAXQTL).
    //
    // Behaviour summary:
    //   - Start from VCF (or reuse provided PGEN/BED/BGEN).
    //   - Optionally subset to expression donors.
    //   - Optionally filter VCF with bcftools.
    //   - Convert into required formats for each downstream tool.
    //   - Calculate genotype PCs unless provided.

    // --- Genotype input handling (VCF vs PGEN) ---
    // Purpose:
    //   Determine the starting genotype source for the pipeline. Either:
    //     • Case 1: A raw VCF file is provided (most common).
    //     • Case 2: No VCF, but a preprocessed PGEN dataset is available.
    //   From here, the workflow ensures genotypes are subset, filtered, and
    //   converted into formats needed by downstream tools.
    //
    // Case 1: VCF provided
    //   - If `subset_genotypes_to_available = true`, restrict VCF to donors present
    //     in the phenotype data (saves memory/compute).
    //   - If `apply_bcftools_filters = true`, run `PREPROCESS_GENOTYPES` to apply
    //     bcftools filters (e.g. MAF, SNP-only, biallelic).
    //   - The resulting VCF (possibly subset and filtered) becomes the
    //     `plink_convert_input` for downstream conversion.
    //
    // Case 2: No VCF
    //   - If a preprocessed PGEN folder is provided, use that directly (prioritised).
    //   - Otherwise, initialise `plink_convert_input` as empty.
    //     (Downstream blocks will handle the missing input accordingly.)

    if (params.input_vcf){
        // VCF file processing
         // Case 1: VCF provided → optionally subset and filter, then use as input
        log.info "---- VCF file provided - ${params.input_vcf} Lets preprocess genotypes and aply any bcftools filters ---"
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
        // Case 2: No VCF → check for preprocessed PGEN, else empty input
        if (params.genotypes.preprocessed_pgen_file==''){
            plink_convert_input=Channel.of()
        }else{
            log.info "---- pgen file provided, will prioritise this over VCF---"
            plink_convert_input=Channel.from(params.genotypes.preprocessed_pgen_file)
        }
    }

    // This block is necesary for Saige as it needs BED for GRM construction and if we dont use dosages we also need this for the QTL matrix construction
    // --- Produce / load PLINK1 BED datasets (for QTL tests and GRM) ---
    // Triggered when:
    //   • SAIGE is enabled (needs BED for GRM construction), or
    //   • `use_gt_dosage = false` (hard calls required), or
    //   • JAXQTL is enabled (expects BED input).
    //
    // Purpose:
    //   Ensure we have PLINK1 BED datasets for two distinct uses:
    //     1) Association testing (general BED set) → `bim_bed_fam`, `plink_path_bed`
    //     2) GRM construction (often pruned/filtered differently) → `bim_bed_fam__GRM`, `plink_path_bed__GRM`
    //
    // Behaviour:
    //   - If no preprocessed BED folder is provided:
    //       • From VCF: run `PLINK_CONVERT` (assoc BED) and `PLINK_CONVERT__GRM` (GRM BED).
    //       • From PGEN: run `PGEN_TO_BED_CONVERT_FOR_QTLS` (assoc BED) and
    //                    `PGEN_TO_BED_CONVERT_FOR_GRM` (GRM BED).
    //     Rationale: GRM construction typically uses a differently filtered/pruned marker
    //     set than association testing (e.g., LD-pruned, MAF/HWE thresholds), so we
    //     create two BEDs.
    //   - If a preprocessed BED folder is provided:
    //       • Reuse it for both association and GRM paths (warning: same BED for both).
    //         This is acceptable but may not be optimal if GRM requires a different
    //         marker set than association testing.
    //
    // Outputs:
    //   - `bim_bed_fam`, `plink_path_bed` → used by association tools (SAIGE step2, JAXQTL, etc.).
    //   - `bim_bed_fam__GRM`, `plink_path_bed__GRM` → used for GRM/sparse GRM construction (SAIGE step1).

    if (params.SAIGE.run || (params.genotypes.use_gt_dosage == false) || params.JAXQTL.run || params.QUASAR.run) {
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
            // TODO: need GRM and association PLINK1 option
        }
    }

    // --- Ensure PGEN (dosage) availability & choose main genotype format ---
    //
    // Purpose:
    //   PGEN format is required when using dosages (`use_gt_dosage = true`).
    //   This block guarantees we always have a usable PGEN path, and then
    //   decides whether the main downstream genotype input (`plink_path`)
    //   should be PGEN (dosages) or BED (hard calls).
    //
    // Behaviour:
    //   - If no preprocessed PGEN is provided:
    //       • Convert from VCF → PGEN (`PGEN_CONVERT`).
    //   - If preprocessed PGEN exists:
    //       • Use it directly.
    //
    //   After that, select genotype format for downstream tools:
    //     • If `use_gt_dosage = true` → set `plink_path = plink_path_pgen`.
    //     • If `use_gt_dosage = false` → set `plink_path = plink_path_bed`.
    //
    // Notes:
    //   - BED is required for tools that don’t support dosages (e.g. JAXQTL, GRM for SAIGE).
    //   - PGEN is required for dosage-based analyses (TensorQTL, LIMIX, SAIGE with DS field).

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

    // --- Genotype principal components (PCs) ---
    //
    // Purpose:
    //   Control for population structure in QTL mapping by including genotype PCs
    //   as covariates. PCs are derived from the genotype matrix (PGEN or BED),
    //   unless a precomputed file is supplied.
    //
    // Behaviour:
    //   - If `params.covariates.genotype_pcs_file` is empty:
    //       • Run `GENOTYPE_PC_CALCULATION` on `plink_path` to compute PCs.
    //       • Output is a TSV with top genotype PCs per sample.
    //   - If a file path is provided:
    //       • Load and use the user-supplied PCs directly.
    //
    // Notes:
    //   - The number of PCs used downstream is controlled by
    //     `params.covariates.nr_genotype_pcs` (see config).
    //   - Precomputing PCs is useful for large cohorts or repeated runs, as
    //     calculating them from scratch can be time-consuming.

    if (params.covariates.genotype_pcs_file==''){
        log.info "--- Calculating genotype PCs ---"
        GENOTYPE_PC_CALCULATION(plink_path)
        genotype_pcs_file = GENOTYPE_PC_CALCULATION.out.gtpca_plink
    }else{
        log.info "--- Genotype PCs file already provided, will use this ---"
        genotype_pcs_file = Channel.fromPath(params.covariates.genotype_pcs_file)
    }

    // --- BGEN generation for LIMIX / SAIGE ---
    //
    // Purpose:
    //   Provide genotype input in the correct format for LIMIX and SAIGE analyses.
    //   Both tools can consume BGEN when using dosages, or BED when using hard calls.
    //
    // Behaviour:
    //   - If `use_gt_dosage = true` and either LIMIX or SAIGE is enabled:
    //       • If no preprocessed BGEN is provided, convert from VCF/PGEN → BGEN
    //         using `BGEN_CONVERT`.
    //       • If a BGEN path is provided, reuse it directly.
    //       • Resulting BGEN path is used for both LIMIX (`plink_path_limix`)
    //         and SAIGE (`genotypes_saige`).
    //   - If `use_gt_dosage = false` but LIMIX or SAIGE is enabled:
    //       • Fall back to BED input (`plink_path_bed`) for both tools.
    //
    // Notes:
    //   - BGEN is required for LIMIX regardless of dosage setting, but here we
    //     generate it only when using dosages.
    //   - SAIGE can run from either BGEN (dosages) or BED (hard calls), depending
    //     on `use_gt_dosage`.

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



    // --- Phenotype normalisation and phenotype PCs ---
    //
    // Purpose:
    //   Prepare expression phenotype files for QTL mapping by normalising counts
    //   and optionally extracting phenotype PCs as covariates. This step ensures
    //   comparability across samples and reduces technical variation.
    //
    // Behaviour:
    //   - Run `NORMALISE_and_PCA_PHENOTYPE` on the aggregated phenotype data
    //     (`phenotype_condition`). Normalisation method and filtering strategy are
    //     controlled by config parameters:
    //       • `filter_method` (None, HVG, filterByExpr)
    //       • `norm_method` (DESEQ, TMM, NONE)
    //       • `inverse_normal_transform`
    //       • `percent_of_population_expressed`
    //       • `use_sample_pca`
    //   - Emit phenotype matrices in BED-compatible format (`for_bed`).
    //   - Expand over the requested set of phenotype PC counts
    //     (`nr_phenotype_pcs`, e.g. "2,4") by combining each phenotype matrix with
    //     a different PC setting (`test123`).
    //   - Run `SUBSET_PCS` to prepare phenotype+PC combinations for downstream
    //     QTL tools.
    //
    // Notes:
    //   - This stage mirrors genotype PCs: adds phenotype-level PCs as covariates.
    //   - Supports multiple runs with different numbers of PCs to optimise model fit.
    //   - Future extensions could include MBV (QTLTools) or RASCAL methods.

    
    log.info "--- Normalising agregated penotype file  using: Filtering method: ${params.filter_method} Norm method: ${params.method} Inverse Transform: ${params.inverse_normal_transform} Norm method:  ${params.norm_method} Percent of Population expressed: ${params.percent_of_population_expressed} Sample PCa: ${params.use_sample_pca}---"
    NORMALISE_and_PCA_PHENOTYPE(phenotype_condition_remap)
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
        test123_fixed.combine(covariates_by_name,by:0).set{for_covs_merge}
        MERGE_COVARIATES(for_covs_merge)
        for_bed = MERGE_COVARIATES.out.for_bed_covs
    }else{
        for_bed = SUBSET_PCS.out.for_bed
    }

    // create GRM for SAIGE or QUASAR
    if (params.QUASAR.run || params.SAIGE.run){
                if (params.existing_sparse_grm){
            sparseGRM = Channel
                .fromPath(params.existing_sparse_grm + "/sparseGRM_*.mtx")
                .ifEmpty { error " No sparseGRM .mtx file found in ${params.existing_sparse_grm}" }

            sparseGRM_sample = Channel
                .fromPath(params.existing_sparse_grm + "/sparseGRM_*.sampleIDs.txt")
                .ifEmpty { error " No sparseGRM sample ID file found in ${params.existing_sparse_grm}" }

        }else{
            CREATE_SPARSE_GRM(plink_path_bed__GRM)
            sparseGRM = CREATE_SPARSE_GRM.out.sparseGRM
            sparseGRM_sample = CREATE_SPARSE_GRM.out.sparseGRM_sample
        }

    }

    // --- LIMIX QTL mapping ---
    // Controlled by: `params.LIMIX.run`.
    //
    // Purpose:
    //   Run eQTL mapping with the LIMIX linear mixed model framework. LIMIX
    //   accounts for relatedness and population structure by including a kinship
    //   matrix, making it suitable for datasets with complex sample structure.
    //
    // Behaviour:
    //   - Compute a kinship matrix from the genotype data (`KINSHIP_CALCULATION`),
    //     required as a random effect in LIMIX models.
    //   - Reformat the phenotype+covariate channel (`for_bed`) into LIMIX’s expected
    //     input structure (`filtered_pheno_channel`).
    //   - Run `LIMIX_eqtls` using:
    //       • Phenotype channel (expression + PCs + covariates)
    //       • Genotype path (`plink_path_limix`, usually BGEN or BED)
    //       • Genome annotation (GTF-based gene coordinates)
    //       • Kinship matrix
    //
    // Notes:
    //   - LIMIX uses permutations (configurable via `params.LIMIX.numberOfPermutations`)
    //     to estimate significance thresholds.
    //   - Kinship calculation can be computationally heavy on large cohorts;
    //     precomputed matrices may be supported in future extensions.

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

    // --- TensorQTL / JAXQTL ? QUASAR preparation and execution ---
    //
    // Purpose:
    //   Prepare covariate and expression BED files in the correct format for
    //   TensorQTL and JAXQTL, then run the selected tools. Both tools require
    //   expression matrices (per gene), sample covariates, and genotype input.
    //
    // Behaviour:
    //   - Reformat `for_bed` into `for_bed_channel`, grouping phenotypes and
    //     covariates in chunks.
    //   - Run `PREPERE_COVARIATES` to generate covariate files (includes genotype PCs).
    //   - Run `PREPERE_EXP_BED` to convert expression matrices into BED format
    //     using genome annotation.
    //   - Combine covariates and expression BEDs → `tensorqtl_input`.
    //   - Clean filenames in `tensorqtl_input` to remove PC suffixes
    //     (`tensorqtl_input_ch_cleaned`), required by TensorQTL.
    //
    // Execution:
    //   - If `params.TensorQTL.run = true`:
    //       • Run `TENSORQTL_eqtls` with cleaned input and chosen genotype path
    //         (`plink_path`, PGEN if dosages / BED otherwise).
    //   - If `params.JAXQTL.run = true`:
    //       • Run `JAXQTL_eqtls` with original input and BED genotypes.
    //
    // Notes:
    //   - TensorQTL supports both cis and trans-QTLs, conditional analysis,
    //     and interaction models (see config).
    //   - JAXQTL is optimised for speed and GPU usage, but only supports BED input.
    //   - This block harmonises covariates + expression formats so both tools
    //     can be run consistently from the same workflow.
    //   - QUASAR requires a specific header on the phenotype bed file, which is not zipped
    //   - QUASAR and SAIGE both require GRM  


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
    
    if (params.QUASAR.run){
        //the inputs we need for quasar are phenoype covariats, plink prefix, phenotype file,grm file (optional), mode and model 
        //we need the correct inputs depending on the model used (dSum or dMean)
        //create channels for Quasar inputs
        dSum_ch = tensorqtl_input.filter{ it[0].startsWith('dSum') }
        dMean_ch = tensorqtl_input.filter{ it[0].startsWith('dMean') }


        quasar_params = Channel.of(params.QUASAR.model, params.QUASAR.mode)
        dMean_ch_quasar = dMean_ch.combine(genome_annotation)
                                    .combine(bim_bed_fam)
                                    .combine(sparseGRM)
                                    .combine(sparseGRM_sample)
                                    .combine(quasar_params.collate(2).first())

        dSum_ch_quasar = dSum_ch.combine(genome_annotation)
                                    .combine(bim_bed_fam)
                                    .combine(sparseGRM)
                                    .combine(sparseGRM_sample)
                                    .combine(quasar_params.collate(2).first())

        if (params.QUASAR.model == 'lm' || params.QUASAR.model == 'lmm'){
            QUASAR(
                dMean_ch_quasar
            )
        } else {
            QUASAR(
                dSum_ch_quasar
            )
        }

    }

    // --- SAIGE single-cell QTL mapping ---
    // Controlled by: `params.SAIGE.run` and `params.method == 'single_cell'`.
    //
    // Purpose:
    //   Run SAIGE for single-cell eQTL mapping. SAIGE is a mixed-model approach
    //   designed to handle relatedness and case–control imbalance. In this workflow,
    //   it operates on pseudo-bulked single-cell expression matrices.
    //
    // Behaviour:
    //   - Only triggered in `single_cell` mode when `params.SAIGE.run = true`.
    //   - Genotype–phenotype mapping file:
    //       • If provided, load and use directly.
    //       • If not provided, create a dummy channel (placeholder required by module).
    //   - Extra covariates file:
    //       • If provided, include in covariates.
    //       • If not, substitute with a dummy file.
    //   - Run `PREP_SAIGE_COVS` to merge genotype PCs with extra covariates.
    //   - Execute `SAIGE_qtls` with:
    //       • Covariates (`covs_Saige`)
    //       • Expression data (`adata`)
    //       • GRM BED (`plink_path_bed__GRM`)
    //       • Association BED (`plink_path_bed`)
    //       • Dosage file (BGEN/PGEN, `genotypes_saige`)
    //       • Genome annotation
    //       • Genotype–phenotype mapping file
    //
    // Notes:
    //   - Dummy files (`fake_file.fq`) are used to satisfy SAIGE input channels
    //     when optional inputs are not provided.
    //   - Requires both association genotypes and GRM genotypes, as SAIGE uses
    //     the GRM during null model fitting (step 1).
    //   - Config parameters (e.g. minMAC, SPAcutoff, useSparseGRM*) control the
    //     detailed behaviour of SAIGE runs.

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

            SAIGE_qtls(covs_Saige,adata,sparseGRM,sparseGRM_sample,plink_path_bed,genotypes_saige,genome_annotation,saige_genotype_phenotype_mapping_file)
        }
    }

    // Generate plots of comparisons of eQTLs detected by all run methods.


}
