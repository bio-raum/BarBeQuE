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

include { BARBEQUE } from './workflows/barbeque'
include { BUILD_REFERENCES }    from './workflows/build_references'
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_completion'
include { DATABASE } from './subworkflows/database'
include { INTERACTIVE_RESULTS } from './workflows/interactive_results'
include { paramsHelp; paramsSummaryLog; validateParameters } from 'plugin/nf-schema'

workflow {


  if (!workflow.containerEngine) {
    log.info("\033[1;31mRunning with Conda is not recommended in production!\033[0m\n\033[0;31mConda environments are not guaranteed to be reproducible - for a discussion, see https://pubmed.ncbi.nlm.nih.gov/29953862/.\033[0m")
  }

  if (WorkflowPipeline.enabled(params.help)) {
      log.info paramsHelp(command: "nextflow run main.nf")
      System.exit(0)
  }
  if (WorkflowPipeline.enabled(params.helpFull)) {
      System.exit(0)
  }

  validateParameters()
  WorkflowMain.initialise(workflow, params, log)
  if (WorkflowPipeline.enabled(params.list_primers)) {
      PrimerCatalog.printCatalog(log)
  }
  WorkflowPipeline.initialise(params, log)

  // Print summary of supplied parameters
  log.info(paramsSummaryLog(workflow))

  if (WorkflowPipeline.enabled(params.build_references)) {
      BUILD_REFERENCES()
  } else {
      DATABASE()
      BARBEQUE(
          DATABASE.out.db,
          DATABASE.out.versions,
          DATABASE.out.taxdump,
          DATABASE.out.accession_taxonomy,
      )
      if (WorkflowPipeline.enabled(params.interactive)) {
          BARBEQUE.out.consensus.collect() | map { "${params.outdir}" } | INTERACTIVE_RESULTS
      }
  }

  PIPELINE_COMPLETION()
}
