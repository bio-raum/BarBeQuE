/*
Include Modules
*/
include { CRABS_DOWNLOADTAXONOMY }                      from './../modules/crabs/download_taxonomy'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_LRRNA }   from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CYTB }    from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CO1 }     from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MITOFISH }       from './../modules/crabs/download_db'
include { CRABS_IMPORT as CRABS_IMPORT_MIDORI   }       from './../modules/crabs/import'
include { CRABS_IMPORT as CRABS_IMPORT_MITOFISH   }     from './../modules/crabs/import'
include { UNTAR as UNTAR_TAXDUMP }                      from './../modules/untar'

workflow BUILD_REFERENCES {

    main:

    ch_midori_dbs = Channel.from([])

    taxdump = file(params.references.taxdump_url, checkIfExists: true)

    // Untar the downloaded taxdump archive
    UNTAR_TAXDUMP(
        taxdump
    )

    // Download the NCBI taxonomy
    CRABS_DOWNLOADTAXONOMY(
        Channel.from([ build: params.reference_version])
    )

    // Download Midori lrRNA database
    DOWNLOAD_MIDORI_LRRNA(
        Channel.from([ db: "midori_lrrna" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_LRRNA.out.db)

    // Download Midori Cytb database
    DOWNLOAD_MIDORI_CYTB(
        Channel.from([ db: "midori_cytb" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_CYTB.out.db)

    // Download Midori CO1 database
    DOWNLOAD_MIDORI_CO1(
        Channel.from([ db: "midori_co1" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_CO1.out.db)

    DOWNLOAD_MITOFISH(
        Channel.from([ db: "mitofish" ])
    )

    // Import midori databases into crabs format
    CRABS_IMPORT_MIDORI(
        ch_midori_dbs,
        CRABS_DOWNLOADTAXONOMY.out.taxonomy.collect()
    )

    // Import mitofish database into crabs format
    CRABS_IMPORT_MITOFISH(
        DOWNLOAD_MITOFISH.out.db,
        CRABS_DOWNLOADTAXONOMY.out.taxonomy.collect()
    )

    if (params.build_references) {
        workflow.onComplete = {
            log.info 'Installation complete - deleting staged files. '
            workDir.resolve("stage-${workflow.sessionId}").deleteDir()
        }
    }
}

