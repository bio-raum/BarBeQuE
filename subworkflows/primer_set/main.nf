/*
Include Modules
*/
include { PARSE_PRIMERS } from './../../modules/parse_primers/main'

workflow PRIMER_SET {

    main:
    ch_versions = channel.from([])

    // Resolve each requested --primer_set name against the live FooDMe2 catalog (PrimerCatalog,
    // lib/PrimerCatalog.groovy) - same catalog --list_primers prints, always current, nothing
    //update the commits as Foodme2 changes or can be out of sync
    def catalog = PrimerCatalog.fetchCatalog()
    def valid_names = catalog.keySet()

    def these_primers = []
    params.primer_set
        .split(',')
        .collect { name -> name.trim() }
        .each { name ->
            if (!valid_names.contains(name)) {
                log.info("Not a valid primer set: ${name}\nValid options are: ${valid_names}\n")
                System.exit(1)
            }
            def sub = PrimerCatalog.fetchSubConfig(catalog[name].config)
            if (!sub.fasta) {
                log.info("Primer set '${name}' has no fasta reference in its upstream config - cannot use it.\n")
                System.exit(1)
            }
            these_primers << [
                [id: name, min: sub.min, max: sub.max],
                file(PrimerCatalog.fastaUrl(sub.fasta), checkIfExists: true),
            ]
        }
    ch_primer_fasta = channel.fromList(these_primers)

    // Each downloaded set becomes a samplesheet, exactly as a --input FASTA would, so
    // workflows/barbeque.nf validates and stages every input route the same way.
    PARSE_PRIMERS(ch_primer_fasta)
    ch_versions = ch_versions.mix(PARSE_PRIMERS.out.versions)

    emit:
    samplesheet = PARSE_PRIMERS.out.samplesheet.map { _meta, tsv -> tsv }
    versions    = ch_versions
}
