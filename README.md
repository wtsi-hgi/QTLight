## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**QTLight** is a bioinformatics best-practice analysis pipeline for eqtl analysis with TensorQTL, SaigeQTL, LIMIX. 
It takes your vcf files (or pgen/bed) alongside flat quantification data (such as bulk RNAseq expression files, ATACseq qantification data, Splicing Quantification data) or a scRNA h5ad file and performs relevant QTL analysis.


<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Logo.png" width="60%"/>
</p>

This pipeline is running TensorQTL and/or LIMIX and/or jaxQTL on bulk and/or SAIGE-qtl on single cell RNA seq datasets and assessed the overlap of the eGenes identified by both methodologies. While TensorQTL is very fast, this methodology uses linear regression which may not be capable in adequately represent the underlying population structure and other covariates, whereas Limix, while very computationally intensive is based on the linear mixed models (LMM) where the kinship matrices can be provided and hence accounting for random effects in a better manner. 


The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/eqtl/results).

### TensorQTL
<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Yascp_workflow-03eqtl_workflow.png" width="100%"/>
</p>

### SaigeQTL
<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Yascp_workflow-03eqtl_workflow_saige.png" width="100%"/>
</p>


### LIMIX
<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Yascp_workflow-03eqtl_workflow_limix.png" width="100%"/>
</p>

### jaxQTL
<p align="center">
  <img src="https://github.com/wtsi-hgi/QTLight/blob/v1.80/assets/images//Yascp_workflow-03eqtl_workflow_jax.png" width="100%"/>
</p>

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
      // - If 'bulk': phenotype_file should point to gene expression count tables (e.g., STAR featureCounts outputs)
  
      input_vcf = '/path/to/genotype.vcf'
      // Optional if using preprocessed genotypes (e.g., PGEN, BED, or BGEN)
      // Leave blank if providing `preprocessed_pgen_file`, `preprocessed_bed_file`, or `preprocessed_bgen_file`
  
      genotype_phenotype_mapping_file = '' 
      // A TSV with three columns: [Genotype_ID    Phenotype_ID    Sample_Category]
      // - Genotype_ID: must match the IID in PLINK .psam/.fam/.pvar
      // - Phenotype_ID: sample ID from expression data
      // - Sample_Category: optional grouping column (e.g., stimulation/timepoint); if not needed, set all values to 'default'
  
      annotation_file = './path/to/annotation.gtf'
      // Required. Defines genomic coordinates of features (e.g., genes, peaks).
      // Accepts either:
      // - A standard GTF file (recommended for gene-level QTLs)
      // - A custom 4-column TSV with no header, containing:
      //   [feature_id   start   end   chromosome]
      //   Example:
      //     ENSG00000160072   1471765   1497848   1
      // The pipeline will extract TSS/midpoint depending on the `position` setting.
  
      phenotype_file = 'path/to/adata.h5ad'
      // - For 'single_cell': must be a .h5ad file containing raw or normalized counts
      // - For 'bulk': can point to a folder with STAR/featureCounts matrices (one per sample)
  
      aggregation_collumn = 'Azimuth:predicted.celltype.l2'
      // Used when method = 'single_cell'
      // This should match a column in the `.obs` of the h5ad file
      // Defines how cells are grouped for pseudobulk aggregation (e.g., by cell type or cluster)
  
      extra_covariates_file = ''
      // Optional: path to a TSV file with additional covariates.
      // - Format: rows = covariate names, columns = sample IDs (matching genotype IDs, i.e., IID)
      // - These covariates will be included alongside principal components in SAIGE/TensorQTL models.
      // Example:
      //
      //         sample   682_683  683_684  684_685  685_686  686_687  687_688  688_689  689_690  690_691  691_692  692_693  693_694
      //         cov1          1        2        0        2        2        3        2        1        0        0        0        1
      //
      // - Sample IDs must match the IID column in the PLINK .psam or genotype_phenotype_mapping_file.
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

QTLight was developed by Matiss Ozols, Anna Cuomo, Marc Jan Bonder, Hannes Ponstingl, Tobi Alegbe, Bradley Harris, Haerin Jang, Vivek Iyer, Nicole Soranzo.
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

