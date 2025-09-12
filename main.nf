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
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_completion'
include { paramsSummaryLog }    from 'plugin/nf-schema'

workflow {

    multiqc_report = Channel.from([])
    if (!workflow.containerEngine) {
        log.info "\033[1;31mRunning with Conda is not recommended in production!\033[0m\n\033[0;31mConda environments are not guaranteed to be reproducible - for a discussion, see https://pubmed.ncbi.nlm.nih.gov/29953862/.\033[0m"
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
    
    PIPELINE_COMPLETION()

}