

# Running bulk qtls
The pipeline supports **bulk QTL mapping** (e.g. bulk RNA-seq, splicing QTLs, ATAC-QTLs, Protein-QTLs).  
Unlike single-cell mode, no pseudobulk aggregation is performed. Users must provide phenotype matrices directly.
## Input Requirements

1. **Phenotype file**  
   - A tab-separated (TSV) expression matrix.  
   - Rows = features (e.g. genes, peaks, splice junctions), columns = samples.  

2. **Genotype–phenotype mapping file**  
   - TSV with at least three columns:  
     ```
     Genotype_ID    Phenotype_ID    Sample_Category
     ```
     - `Genotype_ID`: must match genotype IDs (IID in PLINK).  
     - `Phenotype_ID`: sample ID from phenotype file.  
     - `Sample_Category`: optional grouping (e.g. tissue, condition).  

3. **Annotation file**  
   - For **gene-level QTLs**: GTF file.  
   - For **splicing QTLs**: junction/feature annotation (TSV with [feature_id, start, end, strand]).  
   - For **ATAC-QTLs**: peak intervals with start/end coordinates.  

4. **Genotype file**  
   - Input can be **VCF/BCF**, or preprocessed **PGEN / BED / BGEN**.  

5. **Optional covariates file**  
   - TSV with sample IDs as rows, covariates as columns (e.g. sex, batch).  
   - The pipeline will automatically compute **genotype PCs** and **phenotype PCs**; user covariates can be merged.  

6. **Optional interaction file**  
   - If running interaction QTLs, provide a TSV specifying interaction terms per sample.  

## Example Configs

### TensorQTL (bulk mode)

```groovy
params {
  method = 'bulk'
  phenotype_file = '/path/to/expression.tsv'
  genotype_phenotype_mapping_file = '/path/to/sample_mappings.tsv'
  annotation_file = '/path/to/Homo_sapiens.GRCh38.99.gtf'
  input_vcf = '/path/to/genotypes.vcf'
  norm_method = 'DESEQ'
  inverse_normal_transform = 'FALSE'
  windowSize = 1000000
  outdir = 'results_bulk'
  position = 'TSS' // or Mid
  covariates.nr_phenotype_pcs = '0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20'
  covariates.nr_genotype_pcs = 2
}
```

# Running scRNA qtls
This pipeline supports QTL mapping from single-cell RNA-seq data using multiple backends (TensorQTL, SAIGE-QTL, JAXQTL, Limix).
It handles pseudobulk aggregation automatically if not already done by the user.
## Input Requirements

When running in `single_cell` mode, the following inputs are required:

1. **Phenotype file**  
   - `.h5ad` file (AnnData object) containing raw or normalized counts.  

2. **Genotype–phenotype mapping file**  
   - TSV file with three columns:  
     ```
     Genotype_ID    Phenotype_ID    Sample_Category
     ```
     - `Genotype_ID`: matches IID in PLINK `.psam` / `.fam` / `.pvar`  
     - `Phenotype_ID`: sample ID in expression data  
     - `Sample_Category`: optional grouping (e.g. timepoint). If unused, set to `"default"`.  

3. **Annotation file**  
   - GTF file (recommended), or  
   - Custom 4-column TSV: `[feature_id   start   end   chromosome]`.  

4. **Genotype data**  
   - VCF/BCF or preprocessed formats (PGEN, BED, BGEN).  
   - If preprocessed is used, leave `input_vcf` empty.  

5. **Aggregation column**  
   - `.obs` column from the h5ad used for pseudobulk (e.g. `cell_type`).  

6. **ID columns**  
   - `gt_id_column`: column with individual ID (matching genotype IDs).  
   - `sample_column`: column with sample name (can be identical to `gt_id_column`). 

## Engines

### TensorQTL


Fast cis/trans QTL mapping (GPU/CPU).

```
params {
  method = 'single_cell'
  phenotype_file = '/path/to/data.h5ad'
  annotation_file = '/path/to/genes.gtf'
  genotype_phenotype_mapping_file = '/path/to/geno_pheno_mapping.tsv'
  aggregation_columns = 'cell_type'
  gt_id_column = 'Vacutainer ID'
  sample_column = 'pheno_id'

  TensorQTL.run = true
  aggregation_method = 'dMean,dSum'
  inverse_normal_transform = true
  windowSize = 500000
  numberOfPermutations = 1000
}
```

## SaigeQtl
```
params {
  method = 'single_cell'
  phenotype_file = '/path/to/data.h5ad'
  annotation_file = '/path/to/genes.gtf'
  genotype_phenotype_mapping_file = '/path/to/geno_pheno_mapping.tsv'
  aggregation_columns = 'cell_type'

  SAIGE {
    run = true
    nr_expression_pcs = 5
    minMAF = 0.05
    minMAC = 20
    SPAcutoff = 10000
    cis_trans_mode = 'cis'
  }
}
```

## JAXqtl
```
params {
  method = 'single_cell'
  phenotype_file = '/path/to/data.h5ad'
  annotation_file = '/path/to/genes.gtf'
  genotype_phenotype_mapping_file = '/path/to/geno_pheno_mapping.tsv'
  aggregation_columns = 'cell_type'
  dMean_norm_method = 'NONE'
  aggregation_method = 'dSum' // Saige needs summed counts for the underlying model.
  inverse_normal_transform = 'FALSE'
  n_min_cells = '5'
  n_min_individ = '25'
  covariates.nr_phenotype_pcs = '4'
  covariates.nr_genotype_pcs = 10
  covariates.adata_obs_covariate = 'Sequencing time'

  JAXQTL {
    run = true
    use_gpu = true
    number_of_genes_per_chunk = 2000
    analysis_subentry = 'CD14_mono,CD4_T_CM'
  }
}
```

## Limix
