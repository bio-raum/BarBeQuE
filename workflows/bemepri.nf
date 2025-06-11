// Modules
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'
include { CRABS_INSILICOPCR }           from './../modules/crabs/insilico_pcr'
include { CRABS_DEREPLICATE }           from './../modules/crabs/dereplicate'
include { CRABS_FILTER }                from './../modules/crabs/filter'

workflow BEMEPRI {

    main:

    ch_multiqc_config = params.multiqc_config   ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
    ch_multiqc_logo   = params.multiqc_logo     ? Channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : Channel.value([])

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])

    samplesheet = params.input ? Channel.fromPath(file(params.input, checkIfExists:true)) : Channel.value([])

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    ch_dbs = Channel.from([])
    if (params.custom_db) {
        ch_dbs = Channel.from(
            [[ "id": "custom", "db": file(params.custom_db, checkIfExists: true) ]]
        )
    } else if (params.dbs) {
        these_dbs = []
        valid_databases = params.references.keySet()
        params.dbs.split(",").collect{ it.toLowerCase()}.each { db ->
            if (!valid_databases.contains(db)) {
                log.info "Not a valid database: ${db}\nValid options are: ${valid_databases}\n"
                System.exit(1)
            }
            these_dbs << [ ["id": db, ], file(params.references[db].db, checkIfExists: true)  ]
        }
        ch_dbs = Channel.fromList(these_dbs)
    }

    // Check if the samplesheet is valid
    INPUT_CHECK(samplesheet)

    /*
     Combine each primer set with all requested database
     [ meta, database_meta, database_path ]
    */

    INPUT_CHECK.out.primers.combine(ch_dbs).map { m,n,d ->
        [
            [ 
                primer: m.primer,
                fwd: m.fwd,
                rev: m.rev,
                min: m.min,
                max: m.max,
                db: n.id
            ], d
        ]
    }.set { ch_primers_with_db }

    ch_primers_with_db.view()

    // perform insilico pcr
    CRABS_INSILICOPCR(
        ch_primers_with_db
    )

    // dereplicate in-silico amplicons
    CRABS_DEREPLICATE(
        CRABS_INSILICOPCR.out.txt
    )

    // Filter hits
    CRABS_FILTER(
        CRABS_DEREPLICATE.out.txt
    )

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    MULTIQC(
        multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )

    emit:
    qc = MULTIQC.out.html
}
