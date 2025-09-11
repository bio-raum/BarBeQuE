// Modules
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'
include { CRABS_INSILICOPCR }           from './../modules/crabs/insilico_pcr'
include { CRABS_DEREPLICATE }           from './../modules/crabs/dereplicate'
include { CRABS_FILTER }                from './../modules/crabs/filter'
include { CRABS_SUBSET }                from './../modules/crabs/subset'
include { CRABS_DIVERSITY_FIGURE }      from './../modules/crabs/diversity_figure'
include { VSEARCH_CLUSTER_FAST }        from './../modules/vsearch/cluster_fast'
include { CRABS_AMPLIFICATION_EFFICENCY_FIGURE } from './../modules/crabs/amplification_efficency_figure'
include { CRABS_AMPLICON_LENGTH_FIGURE }from './../modules/crabs/amplicon_length_figure'
include { HELPER_CLUSTER_CONSENSUS }    from './../modules/helper/cluster_consensus'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file'
include { HELPER_TAXONOMIC_COVERAGE }   from './../modules/helper/taxonomic_coverage'

workflow BARBEQUE {

    main:

    ch_multiqc_config = params.multiqc_config   ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
    ch_multiqc_logo   = params.multiqc_logo     ? Channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : Channel.value([])

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])

    samplesheet = params.input ? Channel.fromPath(file(params.input, checkIfExists:true)) : Channel.value([])

    // The pre-installed taxdump folder
    ch_taxdump = file(params.references.taxdump)

    pipeline_settings = Channel.fromPath(dumpParametersToJSON(params.outdir)).collect()

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    ch_dbs = Channel.from([])
    if (params.custom_db) {
        ch_dbs = Channel.from(
            [[ "id": "custom", "db": file(params.custom_db, checkIfExists: true) ]]
        )
    } else if (params.dbs) {
        these_dbs = []
        valid_databases = params.references.databases.keySet()
        params.dbs.split(",").collect{ it.toLowerCase()}.each { db ->
            if (!valid_databases.contains(db)) {
                log.info "Not a valid database: ${db}\nValid options are: ${valid_databases}\n"
                System.exit(1)
            }
            these_dbs << [ ["id": db, ], file(params.references.databases[db].db, checkIfExists: true)  ]
        }
        ch_dbs = Channel.fromList(these_dbs)
    }

    // Check if the samplesheet is valid
    INPUT_CHECK(samplesheet)

    // Copy the samplesheet to the results folder
    STAGE_SAMPLESHEET(samplesheet)

    /*
     Combine each primer set with all requested databases
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

    // perform insilico pcr, takes: [meta, database]
    CRABS_INSILICOPCR(
        ch_primers_with_db
    )
    ch_versions = ch_versions.mix(CRABS_INSILICOPCR.out.versions)

    // dereplicate in-silico amplicons, takes [meta, txt]
    CRABS_DEREPLICATE(
        CRABS_INSILICOPCR.out.txt
    )
    ch_versions = ch_versions.mix(CRABS_DEREPLICATE.out.versions)

    // Filter hits, takes [meta, txt]
    CRABS_FILTER(
        CRABS_DEREPLICATE.out.txt
    )
    ch_versions = ch_versions.mix(CRABS_FILTER.out.versions)

    // fast clustering of crabs OTUs
    VSEARCH_CLUSTER_FAST(
        CRABS_FILTER.out.fasta
    )
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)

    // Cluster consensus
    HELPER_CLUSTER_CONSENSUS(
        VSEARCH_CLUSTER_FAST.out.uc.join(CRABS_FILTER.out.txt),
        ch_taxdump
    )
    ch_versions = ch_versions.mix(HELPER_CLUSTER_CONSENSUS.out.versions)

    // If a taxon is provided, perform additional visualisation/filtering
    if (params.taxon) {

        // Analyse the coverage of the desired taxonomic level
        HELPER_TAXONOMIC_COVERAGE(
            HELPER_CLUSTER_CONSENSUS.out.txt.map { m,t ->
                tuple(m.db, m, t)
            }.combine(
                ch_dbs.map { m,d ->
                    tuple(m.id,d)
                }, by: 0
            ).map { k, m, s, d ->
                tuple(m,s,d)
            },
            params.taxon
        )

        // Generate a subset based on the --taxon argument
        CRABS_SUBSET(
            CRABS_FILTER.out.txt,
            params.taxon
        )
        
        // Visualize the length distribution of putative amplicons
        CRABS_AMPLICON_LENGTH_FIGURE(
            CRABS_SUBSET.out.txt
        )

        // Visualize diversity of amplicons
        CRABS_DIVERSITY_FIGURE(
            CRABS_SUBSET.out.txt
        )

        // Combine each subset with the correct database
        CRABS_SUBSET.out.txt.map {m, s ->
            tuple(m.db,m,s)
        }.combine(
            ch_dbs.map { m,d ->
                tuple(m.id,d)
            }, by: 0
        ).map { k, m, s, d ->
            tuple(m,s,d)
        }.set { ch_amplicons_with_db }

        // visualize amplification efficency
        CRABS_AMPLIFICATION_EFFICENCY_FIGURE(
            ch_amplicons_with_db,
            params.taxon
        )
    }

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

// turn the params map to a JSON file
def dumpParametersToJSON(outdir) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def filename  = "params_${timestamp}.json"
    def temp_pf   = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr   = groovy.json.JsonOutput.toJson(params)
    temp_pf.text  = groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
    return file("${outdir}/pipeline_info/params_${timestamp}.json")
}