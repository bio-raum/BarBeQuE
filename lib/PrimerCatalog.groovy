//
// Resolves BarBeQuE's primer-set catalog straight from bio-raum/FooDMe2's live repo.
// Fetched and regex-parsed as plain text - never `includeConfig`'d or otherwise executed
// as code - so using this catalog never hands code-execution trust to that repo.
//

import groovy.json.JsonOutput

class PrimerCatalog {

    // Pin the catalog and referenced FASTAs so rebuilding the same pipeline revision
    // cannot silently pick up changed primers from FooDMe2's moving main branch.
    static final String REVISION = 'ad8c8f854cd9a33d4edf995d42a8915eb73f71e8'
    static final String BASE_URL = "https://raw.githubusercontent.com/bio-raum/FooDMe2/${REVISION}"

    // conf/primers.config: name -> [config, description, doi, target, platform]
    static Map fetchCatalog() {
        def text = new URL("${BASE_URL}/conf/primers.config").text
        def catalog = [:]
        text.eachMatch(~/(?s)'([\w-]+)'\s*\{(.*?)\}/) { m ->
            def name = m[1]
            def body = m[2]
            catalog[name] = [
                config     : firstGroup(body, /config\s*=\s*"([^"]+)"/),
                description: firstGroup(body, /description\s*=\s*"([^"]*)"/),
                doi        : firstGroup(body, /doi\s*=\s*"([^"]*)"/),
                target     : splitList(firstGroup(body, /target\s*=\s*\[([^\]]*)\]/)),
                platform   : firstGroup(body, /platform\s*=\s*"([^"]+)"/),
            ]
        }
        return catalog
    }

    // A set's own sub-config (e.g. conf/primers/16S_ILM_ASU184_meat.config): [fasta, min, max]
    static Map fetchSubConfig(String configPath) {
        def text = new URL("${BASE_URL}/${configPath}").text
        return [
            fasta: firstGroup(text, /assets\/primers\/([\w\-.]+\.fasta)/),
            min  : firstGroup(text, /amplicon_min_length\s*=\s*(\d+)/)?.toInteger(),
            max  : firstGroup(text, /amplicon_max_length\s*=\s*(\d+)/)?.toInteger(),
        ]
    }

    static String fastaUrl(String fastaName) {
        return "${BASE_URL}/assets/primers/${fastaName}"
    }

    // --list_primers: print the catalog as JSON and stop.
    static void printCatalog(log) {
        def catalog
        try {
            catalog = fetchCatalog()
        }
        catch (Exception e) {
            log.error("Could not fetch the primer-set catalog from bio-raum/FooDMe2: ${e.message}")
            System.exit(1)
        }
        println(JsonOutput.prettyPrint(JsonOutput.toJson(catalog)))
        System.exit(1)
    }

    // Return one URL per unique primer FASTA used by the catalog.
    static Map fetchFastas() {
        def catalog = fetchCatalog()
        if (!catalog) {
            throw new IllegalStateException('FooDMe2 primer catalog was empty')
        }

        def fastas = [:]
        catalog.values()*.config.findAll().unique().each { configPath ->
            def primer = fetchSubConfig(configPath)
            if (primer.fasta) {
                fastas[primer.fasta] = fastaUrl(primer.fasta)
            }
        }
        if (!fastas) {
            throw new IllegalStateException('FooDMe2 primer catalog contained no FASTAs')
        }
        return fastas
    }

    private static String firstGroup(String text, String pattern) {
        def m = (text =~ pattern)
        return m.find() ? m.group(1) : null
    }

    private static List splitList(String raw) {
        if (!raw) {
            return []
        }
        return raw.split(',').collect { it.trim().replaceAll('"', '') }
    }
}
