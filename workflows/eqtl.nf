/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowEqtl.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
// include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include {PREPROCESS_GENOTYPES} from '../modules/nf-core/modules/preprocess_genotypes/main' 
include {PLINK_CONVERT;PGEN_CONVERT} from '../modules/nf-core/modules/plink_convert/main' 
include {SUBSET_GENOTYPE} from '../modules/nf-core/modules/subset_genotype/main' 
include {KINSHIP_CALCULATION} from '../modules/nf-core/modules/kinship_calculation/main' 
include {GENOTYPE_PC_CALCULATION} from '../modules/nf-core/modules/genotype_pc_calculation/main' 
include {SPLIT_PHENOTYPE_DATA} from '../modules/nf-core/modules/split_phenotype_data/main' 
include {NORMALISE_and_PCA_PHENOTYPE} from '../modules/nf-core/modules/normalise_and_pca/main' 
include {LIMIX_eqtls} from '../modules/nf-core/modules/limix/main'
include {PREPROCESS_SAMPLE_MAPPING} from '../modules/nf-core/modules/preprocess_sample_mapping/main'
include {NORMALISE_ANNDATA; REMAP_GENOTPE_ID} from '../modules/nf-core/modules/normalise_anndata/main'
include {AGGREGATE_UMI_COUNTS; SPLIT_AGGREGATION_ADATA} from '../modules/nf-core/modules/aggregate_UMI_counts/main'
include {CHUNK_GENOME} from '../modules/nf-core/modules/chunk_genome/main'
include {PREPERE_EXP_BED} from '../modules/nf-core/modules/prepere_exp_bed/main'
include {TENSORQTL_eqtls} from '../modules/nf-core/modules/tensorqtl/main'
include {H5AD_TO_SAIGE_FORMAT} from '../modules/nf-core/modules/saige/main'
include {SAIGE_qtls} from '../modules/nf-core/modules/saige/main'
include {SUBSET_PCS} from '../modules/nf-core/modules/covar_processing/main'
// include {OPTIMISE_PCS} from '../modules/nf-core/modules/optimise_pcs/main'
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EQTL {

    log.info 'Lets run eQTL mapping'
    
    // if single cell data then have to prepere pseudo bulk dataset.
    if (params.method=='bulk'){
        log.info '------ Bulk analysis ------'
        genotype_phenotype_mapping_file=params.genotype_phenotype_mapping_file
        phenotype_file=params.phenotype_file

        input_channel = Channel.fromPath(genotype_phenotype_mapping_file, followLinks: true, checkIfExists: true)
        
        input_channel.splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.Genotype)}.unique()
            .set{channel_input_data_table}
        channel_input_data_table = channel_input_data_table.collect().flatten().distinct()
        input_channel.splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->row.Sample_Category}.unique().set{condition_channel}
        condition_channel.subscribe { println "condition_channel.out.adata: $it" }
        SPLIT_PHENOTYPE_DATA(genotype_phenotype_mapping_file,phenotype_file,condition_channel)

        phenotype_condition = SPLIT_PHENOTYPE_DATA.out.phenotye_file
        out2 = SPLIT_PHENOTYPE_DATA.out.phenotye_file

    }else if (params.method=='single_cell'){
        log.info '------ Scrna analysis ------'
        
        if (params.normalise_before_or_after_aggregation=='before'){
            // here we normalise the adata all together per splits
            NORMALISE_ANNDATA(params.phenotype_file)
            // NORMALISE_ANNDATA.out.adata.subscribe { println "NORMALISE_ANNDATA.out.adata: $it" }
            pheno = NORMALISE_ANNDATA.out.adata
        }else{
            pheno = params.phenotype_file
        }
        
        if (params.split_aggregation_adata){
            // pheno.subscribe { println "pheno: $it" }
            SPLIT_AGGREGATION_ADATA(pheno,params.aggregation_columns)
            adata = SPLIT_AGGREGATION_ADATA.out.split_phenotypes.flatten()
        }else{
            adata = pheno
        }   
        
        if (params.normalise_before_or_after_aggregation=='after'){
            // here we normalise the adata per splits
            // adata.subscribe { println "adata: $it" }
            NORMALISE_ANNDATA(adata)
            splits_h5ad = NORMALISE_ANNDATA.out.adata
        }else{
            splits_h5ad = adata
        }

        AGGREGATE_UMI_COUNTS(splits_h5ad,params.aggregation_columns,params.gt_id_column,params.sample_column,params.n_min_cells,params.n_min_individ)
        out2 = AGGREGATE_UMI_COUNTS.out.phenotype_genotype_file.map { data ->
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
        // out2.subscribe { println "out2: $it" }

        if(params.genotype_phenotype_mapping_file!=''){
            // Here user has provided a genotype phenotype file where the provided gt_id_column is contaiming a mapping file instead of actual genotype
            REMAP_GENOTPE_ID(out2,params.genotype_phenotype_mapping_file)
            phenotype_condition = REMAP_GENOTPE_ID.out.remap_genotype_phenotype_mapping
            genotype_phenotype_mapping_file = REMAP_GENOTPE_ID.out.genotype_phenotype_mapping
        }else{
            phenotype_condition = out2
            genotype_phenotype_mapping_file = AGGREGATE_UMI_COUNTS.out.genotype_phenotype_mapping.flatten()
        }

        // genotype_phenotype_mapping_file.subscribe { println "genotype_phenotype_mapping_file: $it" }
        genotype_phenotype_mapping_file.splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.Genotype)}.distinct()
            .set{channel_input_data_table2}
        channel_input_data_table=channel_input_data_table2.collect()
    }

    
    if (params.genotypes.preprocessed_bed_file=='' && params.genotypes.preprocessed_pgen_file==''){
        // USE VCF FILE
        log.info 'Lets preprocess genotypes'
        donorsvcf = Channel.from(params.input_vcf)
        if (params.genotypes.subset_genotypes_to_available){
            // Subset genotypes to available in expression data
            SUBSET_GENOTYPE(donorsvcf,channel_input_data_table.collect())
            subset_genotypes = SUBSET_GENOTYPE.out.samplename_subsetvcf
        }else{
            subset_genotypes = donorsvcf
        }

        if (params.genotypes.apply_bcftools_filters){
            // preprocess vcf files to be in the right format for plink
            PREPROCESS_GENOTYPES(subset_genotypes)
            plink_convert_input=PREPROCESS_GENOTYPES.out.filtered_vcf
        }else{
            plink_convert_input=subset_genotypes  
        }
    }else if (params.genotypes.preprocessed_pgen_file!=''){
        plink_path = Channel.from(params.genotypes.preprocessed_pgen_file)
    }else{
        plink_path = Channel.from(params.genotypes.preprocessed_pgen_file)
    }


    // If use dosages we convert vcf to pgen
    // Otherwise we convert it to bed
    if(params.genotypes.use_gt_dosage){
        if(params.genotypes.preprocessed_pgen_file==''){
            PGEN_CONVERT(plink_convert_input)
            plink_path = PGEN_CONVERT.out.plink_path
        }else{
            plink_path = Channel.from(params.genotypes.preprocessed_pgen_file)
        }

        // Saige S1 needs bed file
        if (params.SAIGE.run){
            PLINK_CONVERT(plink_convert_input)
            bim_bed_fam = PLINK_CONVERT.out.bim_bed_fam
        }

    }else{
        if(params.genotypes.preprocessed_bed_file!=''){
            // USE BED FILE
            plink_path = Channel.from(params.genotypes.preprocessed_bed_file)
            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.bed", followLinks: true)
                .set { bed_files }

            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.bim", followLinks: true)
                .set { bim_files }

            Channel.fromPath("${params.genotypes.preprocessed_bed_file}/*.fam", followLinks: true)
                .set { fam_files }
            bim_bed_fam = bim_files
                        .combine(bed_files)
                        .combine(fam_files)   
        }else{
            // here we dont use dosage for analysis
            PLINK_CONVERT(plink_convert_input)
            plink_path = plink_path
            bim_bed_fam = PLINK_CONVERT.out.bim_bed_fam
        }         

    }
        
    





        // we prioritise vcf to convert to pgen
    //     if(params.genotypes.preprocessed_pgen_file==''){
    //         if(params.genotypes.use_gt_dosage){
    //             if (params.input_vcf!=''){
    //                 PGEN_CONVERT(plink_convert_input)
    //                 plink_path = PGEN_CONVERT.out.plink_path
    //             }
    //         }else{
    //             // here we dont use dosage for analysis
    //             PLINK_CONVERT(plink_convert_input)
    //             plink_path = plink_path
    //         }
    //     }else{
    //         if(params.genotypes.use_gt_dosage){
    //             plink_path = Channel.from(params.genotypes.preprocessed_pgen_file)
    //         }else{
    //             // here we dont use dosage for analysis
    //             plink_convert_input = Channel.from(params.genotypes.preprocessed_pgen_file)
    //             PLINK_CONVERT(plink_convert_input)
    //             plink_path = plink_path
    //         }
    //     }

    // }
    
    // if(params.genotypes.preprocessed_pgen_file!=''){
    //     // USE PGEN FILE
    //     plink_convert_input=Channel.from(genotypes.preprocessed_pgen_file)
    //     PLINK_CONVERT(plink_convert_input)
    //     if(params.genotypes.use_gt_dosage){
    //         plink_path = plink_convert_input
    //     }else{
    //         // here we dont use dosage for analysis
    //         plink_path = PLINK_CONVERT.out.plink_path
    //     }
    //     bim_bed_fam = PLINK_CONVERT.out.bim_bed_fam
    // }
    
    
    // if (params.genotypes.preprocessed_bed_file=='' || params.genotypes.preprocessed_pgen_file==''){

    //     // 2) Generate the PLINK file
    //     PLINK_CONVERT(plink_convert_input)
    //     plink_path = PLINK_CONVERT.out.plink_path
    //     bim_bed_fam = PLINK_CONVERT.out.bim_bed_fam

    //     // we prioritise vcf to convert to pgen
    //     if(params.genotypes.preprocessed_pgen_file==''){
    //         if(params.genotypes.use_gt_dosage){
    //             if (params.input_vcf!=''){
    //                 PGEN_CONVERT(plink_convert_input)
    //                 plink_path = PGEN_CONVERT.out.plink_path
    //             }
    //         }else{
    //             // here we dont use dosage for analysis
    //             plink_path = plink_path
    //         }
    //     }
    // }

    // // // 
    // if (params.LIMIX.run){
        
    // }
    

    // if(params.genotypes.preprocessed_pgen_file==''){
    //     if(params.genotypes.use_gt_dosage){
    //         // Here we are using dosage for anaysis.
    //         if (params.input_vcf==''){
    //             PGEN_CONVERT(plink_path)
    //             plink_path = PGEN_CONVERT.out.plink_path
    //         }

    //     }else{
    //         // here we dont use dosage for analysis
    //         plink_path = plink_path
    //     }
    // }
    if (params.covariates.genotype_pcs_file==''){
        GENOTYPE_PC_CALCULATION(plink_path)
        genotype_pcs_file = GENOTYPE_PC_CALCULATION.out.gtpca_plink
    }else{
        genotype_pcs_file = Channel.fromPath(params.covariates.genotype_pcs_file)
    }
    

    // condition_channel = condition_channel.unique() 

    // 4) Phenotype file preperation including PCs, normalisation
    genome_annotation = Channel.from(params.annotation_file)
    // Prepeare chunking file
    
    // MBV method from QTLTools (PMID 28186259)  
    // RASCAL
    

    NORMALISE_and_PCA_PHENOTYPE(phenotype_condition)

    Channel.of(params.covariates.nr_phenotype_pcs).splitCsv().flatten().set{pcs}
    NORMALISE_and_PCA_PHENOTYPE.out.for_bed.combine(pcs).set{test123}

    SUBSET_PCS(test123)

    // LIMIX QTL mapping method
    if (params.LIMIX.run){
        // 3) Generate the kinship matrix and genotype PCs
        KINSHIP_CALCULATION(plink_path)
        filtered_pheno_channel =SUBSET_PCS.out.for_bed.map { tuple ->  [tuple[3],[[tuple[0],tuple[1],tuple[2]]]].combinations()}.flatten().collate(4)
        CHUNK_GENOME(genome_annotation,filtered_pheno_channel)
        limix_cunks_for_each_cond = CHUNK_GENOME.out.limix_condition_chunking.take(2)
        // limix pipeline is curently not correctly chunked. 
        // Genes should be batched and the regions that they need to be tested on also chunked. 
        // Curently we are testing all the genes for all the possible gene cis windoes.
        LIMIX_eqtls(
            limix_cunks_for_each_cond,
            plink_path,
            KINSHIP_CALCULATION.out.kinship_matrix,
            CHUNK_GENOME.out.filtered_chunking_file
        )
    }


    // TensorQTL QTL mapping method
    if (params.TensorQTL.run){
        for_bed_channel = SUBSET_PCS.out.for_bed.map { tuple ->  [tuple[3],[[tuple[0],tuple[1],tuple[2]]]].combinations()}.flatten().collate(4)
        PREPERE_EXP_BED(for_bed_channel,genome_annotation,genotype_pcs_file)

        TENSORQTL_eqtls(
            PREPERE_EXP_BED.out.exp_bed,
            plink_path,
        )
    }

    // SAIGE SCRNA QTL mapping method
    if (params.method=='single_cell'){
        if (params.SAIGE.run){
            SAIGE_qtls(genotype_pcs_file,adata,bim_bed_fam,genome_annotation)
        }
    }

    // Then run a LIMIX and/or TensorQTL - here have to combine the inputs.
    
    // Generate plots of comparisons of eQTLs detected by both methods.


}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {

    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
