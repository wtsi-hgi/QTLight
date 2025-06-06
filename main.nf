#!/usr/bin/env nextflow
/*
========================================================================================
    QTLight
========================================================================================
    Github : https://github.com/wtsi-hgi/QTLight

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2
include { EQTL } from './workflows/eqtl'

workflow QTLight {
    EQTL ()
}

workflow {
    QTLight ()
}

workflow.onComplete {

    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}