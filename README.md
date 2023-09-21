## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/eqtl** is a bioinformatics best-practice analysis pipeline for eqtl analysis.

<p align="center">
  <img src="https://github.com/wtsi-hgi/eqtl/blob/merged_mo11/assets/images/Logo.png" width="60%"/>
</p>

This pipeline is running TensorQTL and/or LIMIX on bulk and single cell RNA seq datasets and assessed the overlap of the eGenes identified by both methodologies. While TensorQTL is very fast, this methodology uses linear regression which may not be capable in adequately represent the underlying population structure and other covariates, whereas Limix, while very computationally intensive is based on the linear mixed models (LMM) where the kinship matrices can be provided and hence accounting for random effects in a better manner. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/eqtl/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Genotype preperation, filtering and subsetting ([`bcftools`]( https://github.com/single-cell-genetics/limix_qtl ))
2. Genotype conversion to PLINK format and filtering ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
3. Genotype kinship matrix calculation ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
4. Genotype PC calculation ([`PLINK2`]( https://github.com/single-cell-genetics/limix_qtl ))
5. LIMIX eqtl mapping ([`LIMIX`]( https://github.com/single-cell-genetics/limix_qtl ))
6. TensorQTL eqtl mapping ([`TensorQTL`](https://github.com/broadinstitute/tensorqtl))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/eqtl -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.
4. Prepeare the input.nf parameters file:
    ```console
    params{
        method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
        input_vcf ='/path/to/genotype/vcf/file.vcf'
        genotype_phenotype_mapping_file = '' //if bulk RNA seq data is fed in then need a tsv file with 3 columns - [Genotype	RNA	Sample_Category]
        annotation_file = './assets/annotation_file.txt'
        phenotype_file = 'path/to/adata.h5ad' //this should point to h5ad file in a single cell experiments or a star counts matrices for the bulk rna seq data
        aggregation_collumn='Azimuth:predicted.celltype.l2' // for scRNA seq data since we feed in the h5ad file we specify here which obs entry to account for for aggregating data.
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
    nextflow run /path/to/cloned/eqtl -profile sanger -resume -c input.nf
    ```

## Documentation

The nf-core/eqtl pipeline comes with documentation about the pipeline [usage](https://nf-co.re/eqtl/usage), [parameters](https://nf-co.re/eqtl/parameters) and [output](https://nf-co.re/eqtl/output).

## Credits

nf-core/eqtl was originally written by Matiss Ozols.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#eqtl` channel](https://nfcore.slack.com/channels/eqtl) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/eqtl for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

