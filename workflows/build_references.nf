include { DOWNLOAD_AND_COMBINE_DB } from './../modules/download_and_combine/main'
include { DOWNLOAD_MITOFISH; CONVERT_MITOFISH } from './../modules/download_mitofish/main'
include { DOWNLOAD_MIDORI; EXTRACT_MIDORI } from './../modules/download/main'
include { UNTAR as UNTAR_TAXDUMP } from './../modules/untar'
include { GUNZIP as GUNZIP_GENBANK } from './../modules/gunzip'
include { INSTALL_CORE_NT } from './../modules/install_core_nt/main'

workflow BUILD_REFERENCES {

    // Download every database that has URLs in conf/resources.config.
    def downloads = []
    def mitofishDownloads = []
    def midoriMembers = []
    params.references.databases.each { id, db ->
        if (db.archive_member) {
            midoriMembers << [db: id, member: db.archive_member]
        } else if (db.urls) {
            def item = [
                [
                    db: id,
                    dir: "databases/${id}",
                    filename: "${id}.fasta",
                    format: db.format,
                    release: db.release,
                ],
                db.urls,
            ]
            if (id == 'mitofish') {
                if (db.urls.size() != 1) {
                    error 'MitoFish requires exactly one FASTA URL'
                }
                mitofishDownloads << [item[0], db.urls[0]]
            } else {
                downloads << item
            }
        }
    }

    // Install the primer FASTAs used by FooDMe2.
    PrimerCatalog.fetchFastas().each { name, url ->
        downloads << [
            [
                db: "primer_${name}",
                dir: 'primers',
                filename: name,
                format: 'fasta',
                release: "FooDMe2 ${PrimerCatalog.REVISION}",
            ],
            [url],
        ]
    }

    DOWNLOAD_AND_COMBINE_DB(channel.fromList(downloads))
    DOWNLOAD_AND_COMBINE_DB.out.warnings
        .splitText()
        .filter { it.trim() }
        .subscribe { message -> log.warn(message.trim()) }

    DOWNLOAD_MITOFISH(channel.fromList(mitofishDownloads))
    DOWNLOAD_MITOFISH.out.warnings
        .splitText()
        .filter { it.startsWith('WARNING:') }
        .subscribe { message -> log.warn(message.trim()) }
    CONVERT_MITOFISH(DOWNLOAD_MITOFISH.out.gz)

    DOWNLOAD_MIDORI(channel.value(params.references.midori_archive_url))
    DOWNLOAD_MIDORI.out.warnings
        .splitText()
        .filter { it.startsWith('WARNING:') }
        .subscribe { message -> log.warn(message.trim()) }

    def midoriInputs = channel.fromList(midoriMembers)
        .combine(DOWNLOAD_MIDORI.out.archive)
    EXTRACT_MIDORI(midoriInputs)
    EXTRACT_MIDORI.out.warnings
        .splitText()
        .filter { it.startsWith('WARNING:') }
        .subscribe { message -> log.warn(message.trim()) }

    // Install both NCBI taxonomy files together.
    if (WorkflowPipeline.enabled(params.install_taxdump)) {
        UNTAR_TAXDUMP(channel.of([
            [id: 'new_taxdump'],
            file(params.reference_sources.ncbi_taxdump.url),
        ]))
        GUNZIP_GENBANK(channel.of([
            [id: 'nucl_gb.accession2taxid'],
            file(params.reference_sources.ncbi_genbank_accession2taxid.url),
        ]))
    }

    if (WorkflowPipeline.enabled(params.install_core_nt)) {
        INSTALL_CORE_NT(channel.value(
            "${params.reference_base}/barbeque/${params.reference_version}/core_nt"
        ))
    }

    workflow.onComplete = {
        if (workflow.success) {
            log.info("References installed in ${params.reference_base}")
        }
    }
}
