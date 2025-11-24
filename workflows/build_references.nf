/*
Include Modules
*/
include { CRABS_DOWNLOADTAXONOMY }                      from './../modules/crabs/download_taxonomy'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_LRRNA }   from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_SRRNA }   from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CYTB }    from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CO1 }     from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CO2 }     from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MIDORI_CO3 }     from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_MITOFISH }       from './../modules/crabs/download_db'
include { CRABS_DOWNLOADDB as DOWNLOAD_METAFISH }       from './../modules/crabs/download_db'
include { CRABS_IMPORT as CRABS_IMPORT_MIDORI   }       from './../modules/crabs/import'
include { CRABS_IMPORT as CRABS_IMPORT_MITOFISH   }     from './../modules/crabs/import'
include { CRABS_IMPORT as CRABS_IMPORT_METAFISH   }     from './../modules/crabs/import'
include { CRABS_IMPORT as CRABS_IMPORT_REFSEQ   }       from './../modules/crabs/import'
include { UNTAR as UNTAR_TAXDUMP }                      from './../modules/untar'
include { GUNZIP as GUNZIP_REFSEQ }                     from './../modules/gunzip'

workflow BUILD_REFERENCES {

    main:

    ch_midori_dbs = Channel.from([])

    taxdump = Channel.fromPath(file(params.references.taxdump_url, checkIfExists: true))
    refseq_mito = Channel.fromPath(file("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz", checkIfExists: true))

    GUNZIP_REFSEQ(
        refseq_mito.map {f ->
            [
                [ db: "mitochondrion.1.1.genomic"],
                f
            ]
        }
    )

    // Untar the downloaded taxdump archive
    UNTAR_TAXDUMP(
        taxdump.map { f ->
            def meta = [:]
            meta.id = f.getSimpleName()
            tuple(meta, f)
        }
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

    // Download Midori lrRNA database
    DOWNLOAD_MIDORI_SRRNA(
        Channel.from([ db: "midori_srrna" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_SRRNA.out.db)

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

    // Download Midori CO2 database
    DOWNLOAD_MIDORI_CO2(
        Channel.from([ db: "midori_co2" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_CO2.out.db)

    // Download Midori CO3 database
    DOWNLOAD_MIDORI_CO3(
        Channel.from([ db: "midori_co3" ])
    )
    ch_midori_dbs = ch_midori_dbs.mix(DOWNLOAD_MIDORI_CO3.out.db)

    DOWNLOAD_MITOFISH(
        Channel.from([ db: "mitofish" ])
    )

    // Download metafish-lib
    DOWNLOAD_METAFISH(
        Channel.from([ db: "metafish" ])
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
    // Import metafish database into crabs format
    CRABS_IMPORT_METAFISH(
        DOWNLOAD_METAFISH.out.db,
        CRABS_DOWNLOADTAXONOMY.out.taxonomy.collect()
    )

    // Import RefSeq Mito
    CRABS_IMPORT_REFSEQ(
        GUNZIP_REFSEQ.out.gunzip,
        CRABS_DOWNLOADTAXONOMY.out.taxonomy.collect()
    )

    if (params.build_references) {
        workflow.onComplete = {
            log.info 'Installation complete - deleting staged files. '
            workDir.resolve("stage-${workflow.sessionId}").deleteDir()
        }
    }
}

