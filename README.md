## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**QTLight** is a bioinformatics best-practice analysis pipeline for eqtl analysis with TensorQTL, SaigeQTL, LIMIX, quasar, jaxQTL. 
It takes your vcf files (or pgen/bed) alongside flat quantification data (such as bulk RNAseq expression files, ATACseq qantification data, Splicing Quantification data) or a scRNA h5ad file and performs relevant QTL analysis.


<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Logo.png" width="60%"/>
</p>

This pipeline is running TensorQTL and/or LIMIX and/or jaxQTL on bulk and/or SAIGE-qtl on single cell RNA seq datasets and assessed the overlap of the eGenes identified by both methodologies. While TensorQTL is very fast, this methodology uses linear regression which may not be capable in adequately represent the underlying population structure and other covariates, whereas Limix, while very computationally intensive is based on the linear mixed models (LMM) where the kinship matrices can be provided and hence accounting for random effects in a better manner. 


The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.

### QTLight workflow overview
<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/main/assets/images/eqtl_workflow.png" width="100%"/>
</p>

> **Figure 1. Overview of the QTLight workflow.**  
> Input genotypes (VCF, PLINK binary, or BGEN) and phenotype matrices (e.g. single-cell or bulk RNA-seq counts, ATAC-seq peak counts, or proteomics intensities) are processed through modular steps for filtering, normalisation, covariate integration, and format conversion.  
> Outputs are directed into five mapping backends:  
> 🟢 **SAIGE-QTL** – Poisson mixed models robust to case–control imbalance and rare variants  
> 🟠 **TensorQTL** – fast regression framework for large-scale cis/trans scans  
> 🧊 **Limix** – flexible mixed-model inference  
> 🔵 **quasar** – fast C++ QTL mapper supporting quantitative and count-based traits  
> 🔴 **JaxQTL** – GPU-accelerated mapping for high-throughput contexts  
> Coloured lines in the diagram correspond to these engines, indicating the data paths to each backend.




## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Genotype preperation, filtering and subsetting ([`bcftools`]( https://github.com/single-cell-genetics/limix_qtl ))
2. Genotype conversion to PLINK format and filtering ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
3. Genotype kinship matrix calculation ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
4. Genotype and Phenotype PC calculation and QTL mapping with various number of PCs ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
5. LIMIX eqtl mapping ([`LIMIX`]( https://github.com/single-cell-genetics/limix_qtl ))
6. TensorQTL qtl mapping ([`TensorQTL`](https://github.com/broadinstitute/tensorqtl))
7. SAIGE-QTL mapping ([`SAIGE-QTL`](https://github.com/weizhou0/qtl))
7. jaxQTL mapping ([`jaxQTL`](https://github.com/mancusolab/jaxqtl))
8. quasar mapping ([`quasar`](https://jeffreypullin.github.io/quasar/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run /path/to/cloned/QTLight -profile test_bulk,<docker/singularity/institute>
    ```

4. Prepeare the input.nf parameters file:
     ```console
        params {
            method = 'single_cell' 
            // Options: 'single_cell' or 'bulk'
            // - If 'single_cell': phenotype_file must be a .h5ad file (AnnData object)
            // - If 'bulk': phenotype_file should point to raw count matrices (e.g., STAR/featureCounts outputs)
        
            input_vcf = false
            // Optional if using preprocessed genotypes.
            // Leave as false or empty if providing one of:
            //   - params.genotypes.preprocessed_pgen_file
            //   - params.genotypes.preprocessed_bed_file
            //   - params.genotypes.preprocessed_bgen_file
        
            genotype_phenotype_mapping_file = '/path/to/geno_pheno_mapping.tsv'
            // Required. TSV file with:
            //   [Genotype_ID    Phenotype_ID    Sample_Category]
            // - Genotype_ID: must match PLINK IID (in .psam/.fam/.pvar)
            // - Phenotype_ID: must match sample ID in h5ad `.obs`
            // - Sample_Category: optional grouping label (e.g., 'default', 'stimA')
        
            annotation_file = '/path/to/annotation.gtf'
            // Required. Gene annotation in GTF format OR custom 4-column TSV:
            //   [feature_id  start  end  chromosome]
            // The coordinate used (TSS vs midpoint) is controlled by `position`
        
            phenotype_file = '/path/to/input_expression.h5ad'
            // For 'single_cell': must be an .h5ad file with raw or normalized counts
            // For 'bulk': a gene expression matrix (TSV)
        
            aggregation_columns = 'cell_type'
            // Comma-separated column(s) in `.obs` used for pseudobulk aggregation
            // E.g., 'cell_type', 'Azimuth:predicted.celltype.l2'
        
            aggregation_subentry = ''
            // Optional. If provided, restricts analysis to these sublevels within aggregation_columns
            // E.g., 'Mono,B,Platelet'
        
            aggregation_method = 'dMean,dSum'
            // Aggregation methods to apply: dMean = average expression, dSum = summed counts
            // Can provide both, comma-separated
        
            split_aggregation_adata = true
            // Whether to split .h5ad by Sample_Category before aggregating
        
            gt_id_column = 'Vacutainer ID'
            // Column in `.obs` with the **donor/genotype ID**.
            // Must match the VCF/PLINK ID or the `RNA` column in the genotype–phenotype mapping file.

            sample_column = 'pheno_id'
            // Column in `.obs` with the **sample/library ID**.
            // Distinguishes multiple measurements from the same donor.
            // Can be the same as `gt_id_column` if each sample maps to one donor.
        
            norm_method = 'NONE'
            // Normalisation strategy for bulk datasets: DESEQ | TMM | NONE

            dMean_norm_method = 'cp10k'
            // Normalization method to apply before dMean aggregation.
            // Options:
            //   - 'cp10k'         : Total-count normalize to 10,000 UMIs/cell, then log1p
            //   - 'pf_log1p_pf'   : Pseudofactor normalization → log1p → pseudofactor again
            //   - 'NONE'          : No normalization; original file passed through unchanged

            //
            // Notes:
            // - Raw count matrix is expected to be in `adata.X` or `adata.layers['counts']`
            // - If not present, the pipeline assumes `adata.X` is raw and warns the user

     
            filter_method = 'None'
            // Gene filtering strategy before PCA/QTL: HVG | filterByExpr | None
        
            inverse_normal_transform = 'FALSE'
            // Whether to apply inverse normal transform post-normalisation
        
            windowSize = 500000
            // Window size (+/- bp) around gene TSS or midpoint for cis-QTL
        
            percent_of_population_expressed = 0.05
            // Minimum fraction of individuals in which gene must be expressed

           inverse_normal_transform = 'FALSE'
            // Apply inverse normal transformation to data after normalization (if TRUE)
     
            n_min_cells = '5'
            // Minimum cells per individual per celltype to include in QTL
        
            n_min_individ = '25'
            // Minimum individuals with valid expression to include gene
        
            maf = 0.01
            hwe = 0.000001
            numberOfPermutations = 1000
            
            covariates {
                nr_phenotype_pcs = '2,4' 
                // Comma-separated values. Each entry defines how many phenotype PCs to use per model.
            
                nr_genotype_pcs = 4 
                // Number of genotype PCs to include in the model for population structure correction.
            
                genotype_pc_filters = '--indep-pairwise 50 5 0.2'
                // PLINK2 parameters used to calculate genotype PCs if not provided.
            
                genotype_pcs_file = ''
                // Optional. Path to precomputed genotype PCs (TSV)
                // Format: rows = PC names, columns = sample IDs (must match .psam IIDs)
                // Ensure it includes at least `nr_genotype_pcs` components.
            
                extra_covariates_file = ''
                // Optional. Path to a TSV file with additional covariates (numeric only!)
                // These will be added to the model along with PCs.
                //
                // Format:
                //     covariate   S1   S2   S3 ...
                //     Age         35   40   29
                //     BMI         22   27   24
                //
                // - First column: covariate names
                // - First row: header with sample IDs (must match genotype IIDs)
                // - All values must be strictly numeric (no categories, booleans, or NA)
                // - Missing values are not allowed — impute or remove samples upstream.
            }
                    
            genotypes {
                subset_genotypes_to_available = false
                // If true: subset genotype data to only individuals found in expression data
                // (useful for large genotype datasets)
            
                use_gt_dosage = true
                // If true: use genotype dosages (DS field in VCF or PGEN format)
                // If false: use hard-called genotypes (GT field from VCF or PLINK BED)
            
                preprocessed_pgen_file = '/path/to/pgen_dir/'
                // Path to directory containing a PLINK2 dataset: .pgen, .psam, .pvar
                // This should be a clean folder with only one PLINK2 trio.
            
                preprocessed_bed_file = ''
                // Optional: path to PLINK1 dataset (BED format)
                // Folder should contain matching .bed, .bim, .fam
            
                preprocessed_bgen_file = ''
                // Optional: path to BGEN file (for LIMIX only)
                // Must include .bgen, .sample, and .bgi index
            }
        }
    ```
    example genotype_phenotype_mapping_file
    |Genotype	|RNA	|Sample_Category|
    |-----------------|----------|-------------------|
    |HPSI0713i-aehn_22|	MM_oxLDL7159503|	M0_Ctrl|
    |HPSI0713i-aehn_22|	MM_oxLDL7159504|	M0_oxLDL|
    |HPSI0713i-aehn_22	|MM_oxLDL7159505	|M1_oxLDL|



4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```console
    nextflow run /path/to/cloned/QTLight -profile sanger -resume -c input.nf
    ```

## Documentation

The nf-core/eqtl pipeline comes with documentation about the pipeline [usage](./docs/usage.md) and [output](./docs/output.md).

## Credits

QTLight was developed by Matiss Ozols, Tobi Alegbe, Marc Jan Bonder, Hannes Ponstingl, Bradley Harris, Haerin Jang, Vivek Iyer, Nicole Soranzo.
<!-- 
## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#eqtl` channel](https://nfcore.slack.com/channels/eqtl) (you can join with [this invite](https://nf-co.re/join/slack)). -->

## Citations


If you use  nf-core/eqtl for your analysis, please cite it using the following doi: [10.5281/zenodo.15601494](https://doi.org/10.5281/zenodo.15601494)
> Ozols, M. et al. QTLight (Quantitative Trait Loci mapping pipeline): GitHub. https://github.com/wtsi-hgi/QTLight. [![DOI](https://zenodo.org/badge/454434698.svg)](https://doi.org/10.5281/zenodo.15601493)

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->


<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

