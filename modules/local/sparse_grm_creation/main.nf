
process CREATE_SPARSE_GRM {

    tag "$condition, $pcs"

label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.saige_grm_container}"
    } else {
        container "${params.saige_grm_docker}"
    }    


    publishDir  path: "${params.outdir}/SPARSE_GRM",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path(plink)
        
    output:
        path("sparseGRM_*.mtx"), emit: sparseGRM
        path("sparseGRM_*.sampleIDs.txt"), emit: sparseGRM_sample

    script:
    """

        plink_dir="${plink}"
        base_name=""

        if ls "\$plink_dir"/*.bed 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bed | head -n 1) .bed)
        elif ls "\$plink_dir"/*.bim 1> /dev/null 2>&1; then
            base_name=\$(basename \$(ls "\$plink_dir"/*.bim | head -n 1) .bim)
        else
            echo "No .bed or .bim file found in \$plink_dir"
            exit 1
        fi

        echo "Detected base name: \$base_name"


        createSparseGRM.R \
            --nThreads=${task.cpus} \
            --outputPrefix=sparseGRM_output \
            --numRandomMarkerforSparseKin=${params.SAIGE.numRandomMarkerforSparseKin} \
            --relatednessCutoff ${params.SAIGE.relatednessCutoff} \
            --famFile  "\$plink_dir/\$base_name.fam" \
            --bimFile  "\$plink_dir/\$base_name.bim" \
            --bedFile  "\$plink_dir/\$base_name.bed" 
    """

} 