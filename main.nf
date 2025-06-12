#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/**
===============================
BarBeQuE Pipeline - Benchmarking metabarcoding primers
===============================

This Pipeline performs benchmarking of metabarcoding primers

### Homepage / git
git@github.com:bio-raum/barbeque.git

**/

// Pipeline version
params.version = workflow.manifest.version

include { BARBEQUE }            from './workflows/barbeque'
include { BUILD_REFERENCES }    from './workflows/build_references'
include { paramsSummaryLog }    from 'plugin/nf-schema'

workflow {

    def summary = [:]

    multiqc_report = Channel.from([])
    if (!workflow.containerEngine) {
        log.warn "NEVER USE CONDA FOR PRODUCTION PURPOSES!"
    }

    WorkflowMain.initialise(workflow, params, log)
    WorkflowPipeline.initialise(params, log)

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    if (params.build_references) {
        BUILD_REFERENCES()
    } else {
        BARBEQUE()
        multiqc_report = multiqc_report.mix(BARBEQUE.out.qc).toList()
    }
    
    def emailFields = [:]
    emailFields['version'] = workflow.manifest.version
    emailFields['session'] = workflow.sessionId
    emailFields['success'] = workflow.success
    emailFields['dateStarted'] = workflow.start
    emailFields['dateComplete'] = workflow.complete
    emailFields['duration'] = workflow.duration
    emailFields['exitStatus'] = workflow.exitStatus
    emailFields['errorMessage'] = (workflow.errorMessage ?: 'None')
    emailFields['errorReport'] = (workflow.errorReport ?: 'None')
    emailFields['commandLine'] = workflow.commandLine
    emailFields['projectDir'] = workflow.projectDir
    emailFields['script_file'] = workflow.scriptFile
    emailFields['launchDir'] = workflow.launchDir
    emailFields['user'] = workflow.userName
    emailFields['Pipeline script hash ID'] = workflow.scriptId
    emailFields['manifest'] = workflow.manifest
    emailFields['summary'] = summary

    email_info = ''
    emailFields.each { s ->
        email_info += "\n${s.key}: ${s.value}"
    }

    outputDir = new File("${params.outdir}/pipeline_info/")
    if (!outputDir.exists()) {
        outputDir.mkdirs()
    }

    outputTf = new File(outputDir, 'pipeline_report.txt')
    outputTf.withWriter { w -> w << email_info }
}