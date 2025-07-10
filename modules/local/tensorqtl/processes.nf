tensor_label = params.utilise_gpu ? 'gpu' : "process_medium"   

process TENSORQTL {  
    label "${tensor_label}"
    tag "$condition, $interaction, $nr_phenotype_pcs"
    // cache false
    
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.eqtl_container}"
  } else {
    container "${params.eqtl_docker}"
  }
  

  input:
    tuple(
        val(condition),
        path(aggrnorm_counts_bed),
        path(covariates_tsv),
        val(nr_phenotype_pcs),
        path(interaction_file)
    )
    each path(plink_files_prefix)
    val(skip_nominal)
    val(preprocess_bed)
    val(skip_gwas)
    val(map_independent_qtls)

  output:
    tuple val(condition), val(interaction), path("${outpath}"), emit: pc_qtls_path

  script:
    // If a file with interaction terms is provided, use the interaction script otherwise use the standard script   
    if ("${interaction_file}" != 'fake_file.fq') {
      tensor_qtl_script = "tensorqtl_analyse_interaction.py -inter ${interaction_file} --interaction_maf ${params.TensorQTL.interaction_maf}"
      inter_name = file(interaction_file).baseName
      interaction = "${inter_name}"
      outpath = "${nr_phenotype_pcs}/interaction_output/${inter_name}"
      chrom_to_map_trans = ""
    } else {
      tensor_qtl_script = "tensorqtl_analyse.py -nperm ${params.numberOfPermutations}"
      outpath = "${nr_phenotype_pcs}/base_output/base"
      interaction = 'base'
      
      if (params.TensorQTL.chrom_to_map_trans == ''){
        chrom_to_map_trans = ""
      }else{
        chrom_to_map_trans = "--chrom_to_map_trans ${params.TensorQTL.chrom_to_map_trans} "
      }
    }

    if (skip_gwas){
      chrom_to_map_trans = ""
    }

    if (skip_nominal) {
      map_nominal_flag = ""
    } else {
      map_nominal_flag = "--map_nominal"
    }
    if (map_independent_qtls==true) {
      map_independent_qtls = "--map_independent_qtls"
    } else {
      map_independent_qtls = ""
    }

    if (preprocess_bed) {
      preprocess_bed = "bedtools sort -i ${aggrnorm_counts_bed} -header > Expression_Data.sorted.bed; sed -i 's/^chr//' Expression_Data.sorted.bed"
    } else {
      preprocess_bed = ""
    }

    if (params.genotypes.use_gt_dosage) {
      dosage = "--dosage"
    }else{
      dosage = ""
    }
    """
      plink_dir="${plink_files_prefix}"
      base_name=""
      if ls "\$plink_dir"/*.psam 1> /dev/null 2>&1; then
          base_name=\$(basename \$(ls "\$plink_dir"/*.psam | head -n 1) .psam)
          pgen_or_bed="--pfile"
      elif ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
          base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
          pgen_or_bed="--bfile"
      else
          echo "No .psam or .bed file found in \$plink_dir"
          exit 1
      fi
      echo " Detected base name: \$base_name"
      echo "Using mode: \$pgen_or_bed"

      echo ${map_independent_qtls}
      ${preprocess_bed}
      ${tensor_qtl_script} --plink_prefix_path "\$plink_dir/\$base_name" --expression_bed Expression_Data.sorted.bed --covariates_file ${covariates_tsv} -window ${params.windowSize} ${dosage} --maf ${params.maf} --outdir ${outpath} ${map_nominal_flag} ${chrom_to_map_trans} ${map_independent_qtls}
      cd ${outpath} && ln ../../../${covariates_tsv} ./ && ln ../../../Expression_Data.sorted.bed ./ && ln ../../../${interaction_file} ./
      if [[ "$map_nominal_flag" == "--map_nominal" ]]; then
          cat ${outpath}/cis_nominal1.cis_qtl_pairs.1.tsv >> ${outpath}/${condition}_full_sumstats.tsv
          cat ${outpath}/cis_nominal1.cis_qtl_pairs.{2..22}.tsv | grep -v variant_id >> ${outpath}/${condition}_full_sumstats.tsv
      fi
    """
}