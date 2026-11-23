//
// This file holds functions to validate user-supplied arguments
//

class WorkflowPipeline {

    //
    // Is a true/false parameter switched on? Nextflow passes `--mask false` from the
    // command line as the text "false", and any non-empty text counts as true in an
    // if, so `if (params.mask)` would switch masking ON. Converting through the text
    // works for real true/false values from nextflow.config and for command-line text.
    //
    public static boolean enabled(value) {
        return value.toString().toBoolean()
    }

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (enabled(params.build_references) && !params.reference_base) {
            log.error '--build_references requires --reference_base'
            System.exit(1)
        }
        if (enabled(params.list_dbs)) {
            println('Available databases:')
            println('===========================')
            params.references.databases.keySet().each { db ->
                def info = params.references.databases[db].description
                println("Name: ${db}\tSource: ${info}")
                println('---------------------------')
            }
            System.exit(1)
        }
        if (!params.dbs && !params.custom_db && !enabled(params.build_references)) {
            def options = params.references.databases.findAll { _, db -> !db.prebuilt }.keySet()
            log.error("No database selected. Choose one with --dbs <name>: ${options.join(', ')}\nOr provide your own FASTA with --custom_db <file.fasta>")
            System.exit(1)
        }
        if (params.dbs && !params.run_name) {
            log.info 'Must provide a run_name (--run_name)'
            System.exit(1)
        }
        if (params.input && params.primer_set) {
            log.info 'Provide either --input or --primer_set, not both'
            System.exit(1)
        }
        if (!params.input && !params.primer_set && !enabled(params.build_references)) {
            log.info 'Pipeline requires a sample sheet / primer FASTA directory (--input) or a named primer set (--primer_set) as input'
            System.exit(1)
        }
        if (params.input && isFastaInput(params.input) && (!params.primer_min || !params.primer_max)) {
            log.info 'When --input is a primer FASTA, or a directory of them, provide global amplicon bounds with --primer_min and --primer_max'
            System.exit(1)
        }
    }

    //
    // Is --input primer FASTA (a directory of them, or a single FASTA file), rather than a samplesheet?
    // Decided by content rather than by file extension, so a samplesheet is never parsed as a FASTA
    // (or the other way round) just because of how it was named.
    //
    public static boolean isFastaInput(input) {
        def path = new File(input.toString())
        if (path.isDirectory()) {
            return true
        }

        def first_line = null
        try {
            path.withReader { reader ->
                def line
                while (first_line == null && (line = reader.readLine()) != null) {
                    if (line.trim()) {
                        first_line = line.trim()
                    }
                }
            }
        } catch (Exception e) {
            // Let the schema validator handle missing files
            return false
        }
        return first_line != null && first_line.startsWith('>')
    }

}
