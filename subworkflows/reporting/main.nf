/*
Include Modules
*/
include { BARBEQUE_SETTINGS } from './../../modules/helper/barbeque_settings/main'
include { MULTIQC } from './../../modules/multiqc/main'
include { MULTIQC as MULTIQC_COMBINED } from './../../modules/multiqc/main'

workflow REPORTING {

    take:
    ch_multiqc_files     // [meta, file] - every table to show in the reports
    ch_primers_with_db   // [meta, db_fasta] - one per primer x database pair
    ch_software_versions // software_versions_mqc.yml, plus other files the same for every report

    main:
    // MULTIQC declares config and logo as path inputs, which cannot be null, so an
    // unset one is passed as an empty list and stages no file.
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])

    ch_settings_meta = ch_primers_with_db.map { meta, db ->
        meta
    }
    // One settings table and one sample sheet for the whole run
    BARBEQUE_SETTINGS(ch_settings_meta.collect())

    // Files that are the same for every report: software versions, settings, sample sheet.
    // Passed next to the per-primer files, so each report stages them exactly once.
    ch_run_files = ch_software_versions
        .mix(BARBEQUE_SETTINGS.out.mqc, BARBEQUE_SETTINGS.out.samplesheet)
        .collect()

    // Group every collected table by its meta map, so MULTIQC emits one report
    // per primer-database combination rather than one for the whole run.
    MULTIQC(
        ch_multiqc_files.groupTuple(by: 0),
        ch_run_files,
        ch_multiqc_config,
        ch_multiqc_logo,
    )

    // A second report with every primer and database together. All file names already
    // start with <primer>_<db>, so they cannot collide when staged into one task.
    MULTIQC_COMBINED(
        ch_multiqc_files
            .filter { _meta, file ->
                !file.name.endsWith('.taxon_search_mqc.html') &&
                !file.name.endsWith('.amplified_taxa_mqc.tsv')
            }
            .map { _meta, files -> files }
            .flatten()
            .collect()
            .map { files -> [[id: 'all_primers'], files] },
        ch_run_files,
        ch_multiqc_config,
        ch_multiqc_logo,
    )

    emit:
    html = MULTIQC.out.html.mix(MULTIQC_COMBINED.out.html)
}
