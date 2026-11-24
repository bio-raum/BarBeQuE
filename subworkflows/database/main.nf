/*
Include Modules
*/
include { DB_FILTER } from './../../modules/seqkit/db_filter'
//experimental module
include { TAXID_DB_FILTER } from './../../modules/helper/taxid_db_filter'

workflow DATABASE {

    main:
    ch_versions = channel.from([])

    // Resolve the taxonomy resources once for both database filtering and the
    // downstream BARBEQUE workflow. BUILD_REFERENCES installs these two paths
    // together; explicit parameters remain available for custom installations.
    def taxdump_path = params.taxdump ?: params.references.taxdump
    def accession_path = params.accession_taxonomy ?: (params.reference_base
        ? "${params.reference_base}/barbeque/${params.reference_version}/taxonomy/nucl_gb.accession2taxid"
        : null)

    if (!accession_path) {
        log.error("No accession-to-taxonomy mapping found - provide --accession_taxonomy <file>, or rebuild the taxonomy references")
        System.exit(1)
    }

    ch_taxdump = channel.value(
        file(taxdump_path, checkIfExists: true)
    )
    ch_accession_taxonomy = channel.value(
        file(accession_path, checkIfExists: true)
    )

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    these_dbs = []
    if (params.custom_db) {
        // One FASTA, a comma-separated list, or a glob. `files()` always returns a list, so
        // each custom database becomes its own entry and is benchmarked separately - the same
        // shape --dbs produces below.
        params.custom_db
            .split(",")
            .each { pattern ->
                files(pattern.trim(), checkIfExists: true).each { custom_file ->
                    these_dbs << [["id": custom_file.baseName], custom_file]
                }
            }
    }
    else if (params.dbs) {
        valid_databases = params.references.databases.keySet()
        params.dbs
            .split(",")
            .collect { db_str -> db_str.toLowerCase() }
            .each { db ->
                if (!valid_databases.contains(db)) {
                    log.info("Not a valid database: ${db}\nValid options are: ${valid_databases}\n")
                    System.exit(1)
                }
                def info = params.references.databases[db]
                if (info.prebuilt) {
                    log.info("'${db}' is a pre-built BLAST index, not an in-silico-PCR reference - it cannot be used with --dbs.")
                    System.exit(1)
                }
                these_dbs << [
                    ["id": db],
                    file(info.db, checkIfExists: true),
                ]
            }
    }
    // Database names go into every output path, so a collision would overwrite results.
    def duplicate_ids = these_dbs.countBy { entry -> entry[0].id }.findAll { _id, count -> count > 1 }.keySet()
    if (duplicate_ids) {
        log.error("Duplicate database name(s): ${duplicate_ids.join(', ')} - every database needs a unique file basename")
        System.exit(1)
    }

    ch_dbs = channel.fromList(these_dbs)

    // Restrict every selected/custom db to a single taxon (and its descendants) before anything
    // else runs against it - lets --taxid be combined with --dbs or --custom_db to analyse/check
    // just that taxon instead of the whole database.
    if (params.taxid) {
        TAXID_DB_FILTER(ch_dbs, ch_taxdump, ch_accession_taxonomy, params.taxid)
        ch_versions = ch_versions.mix(TAXID_DB_FILTER.out.versions)
        ch_dbs = TAXID_DB_FILTER.out.fasta
    }

    // Apply the configured sequence filters to the selected database.
    if (WorkflowPipeline.enabled(params.db_filter)) {
        DB_FILTER(ch_dbs)
        ch_versions = ch_versions.mix(DB_FILTER.out.versions)
        ch_dbs = DB_FILTER.out.fasta
    }

    emit:
    db = ch_dbs
    versions = ch_versions
    taxdump = ch_taxdump
    accession_taxonomy = ch_accession_taxonomy
}
