params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='path/to/vcf/input.vcf.gz'
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = './assets/annotation_file.txt' //assets file that has start and end positions for the ghenes, this one is using hg38
    phenotype_file = '/path/to/star/counts/or/for/scrna/h5ad/file/.tsv_.h5ad' //this should point to h5ad file in a single cell experiments.
    aggregation_collumn='Azimuth:predicted.celltype.l2' //for the scrna h5ad file define which collum to use to aggregate the data based on.
}