params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='path/to/vcf/input.vcf.gz'
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = './assets/annotation_file.txt' //assets file that has start and end positions for the ghenes, this one is using hg38
    phenotype_file = '/path/to/star/counts/or/for/scrna/h5ad/file/.tsv_.h5ad' //this should point to h5ad file in a single cell experiments, ensuring adata.X is raw counts
    aggregation_columns='Azimuth:predicted.celltype.l2,Celltypist:Immune_All_High,Celltypist:Immune_All_Low' //for the scrna h5ad file define which annotations to use when aggregating the data. Can be one value or multiple comma separated values
    gt_id_column='donor_id' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='convoluted_samplename' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}