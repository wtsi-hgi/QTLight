# Running bulk qtls

For bulk qtls useers must provide:
1) phenotype file - a tsv seperated file
2) genotype - phenotype mapping file which also contains a group
3) annotation files - which can be a gtf file or for splicing or for example for ATAC QTLs these would be interval start end and strand informations.
4) Genotype file in vcf/bcf format.
5) Optionally you could also provide an extra covariates file, pipeline will calculate genotype and phenotype PCs that will be included, but you may also add aditional covariates succh as sex.
6) Provide an aditional file with interactions - If you are running interaction QTLs 

# Running scRNA qtls

For scRNA analysis if you havent aready performed pseudobulk then pipeline can perform this for the users. In scRNA mode of pipeline users must provide:
1) phenotype file which in this case is a h5ad file.
2) genotype - phenotype mapping file
3) anotation file - a gtf file that lists the gene start, end and strand information.
4) Genotype file in vcf/bcf format.
5) Aggregation column information - which adata.obs column to use for agregation - here each of the unique labels in the particular column will be used as a category to test the qtls for.
6) gt_id_column - defines which column has individual level id (corresponding to the VCF file)
7) defines which column has the sample name (can be identical to gt_id_column if sample=individual)