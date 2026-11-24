process BARBEQUE_SETTINGS {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    val(primers)    // every primer meta of the run, one per primer x database

    output:
    path("barbeque_settings_mqc.yml"), emit: mqc
    path("samplesheet_mqc.json"), emit: samplesheet

    script:
    // The settings are the same for the whole run, so there is one settings table per run.
    // Read length is only used when masking is on, so it is left out otherwise.
    def read_length = WorkflowPipeline.enabled(params.mask) ? "\"Read length\": {value: \"${params.read_length} bp\"}" : ""
    // Database filter details are useful only when the filter is active. Keep
    // unset optional limits visible so the report describes the exact command.
    def db_filter_settings = ""
    if (WorkflowPipeline.enabled(params.db_filter)) {
        def filter_values = [
            "Header exclusion pattern": params.db_filter_pattern,
            "Minimum database length": params.db_filter_min_length != null ? "${params.db_filter_min_length} bp" : "Not set",
            "Maximum database length": params.db_filter_max_length != null ? "${params.db_filter_max_length} bp" : "Not set",
            "Maximum ambiguous bases": params.db_filter_max_n != null ? params.db_filter_max_n.toString() : "Not set",
        ]
        db_filter_settings = filter_values.collect { name, value ->
            def quoted_name = groovy.json.JsonOutput.toJson(name)
            def quoted_value = groovy.json.JsonOutput.toJson(value?.toString() ?: "Not set")
            "${quoted_name}: {value: ${quoted_value}}"
        }.join("\n      ")
    }

    // Sample sheet: one row per primer with its amplicon length bounds. A primer run
    // against several databases appears once.
    def samplesheet = [
        id: 'samplesheet',
        section_name: 'Sample sheet',
        plot_type: 'table',
        pconfig: [id: 'samplesheet', title: 'Sample sheet', col1_header: 'Primer', scale: false],
        data: primers.collectEntries { m ->
            [(m.primer): ['Minimum length': m.min.toString(), 'Maximum length': m.max.toString()]]
        },
    ]
    def samplesheet_json = groovy.json.JsonOutput.toJson(samplesheet)
    """
    cat > samplesheet_mqc.json <<'EOF'
    ${samplesheet_json}
    EOF

    cat > barbeque_settings_mqc.yml <<EOF
    id: barbeque_settings
    section_name: BarBeQuE Settings
    plot_type: table

    pconfig:
      id: barbeque_settings
      title: BarBeQuE Settings
      col1_header: Setting
      # No colour bars - these are settings, not results to compare
      scale: false

    headers:
      value:
        title: Value

    # Values are quoted, so MultiQC shows them exactly as written (0.97 stays 0.97).
    data:
      "VSEARCH identity": {value: "${params.cluster_id}"}
      "Consensus fraction": {value: "${params.consensus_fraction}"}
      "OBI-PCR mismatches": {value: "${params.obipcr_mismatches}"}
      "Primer-site not-found threshold": {value: "${params.primer_site_not_found_mismatches}"}
      "Primer-site search window": {value: "${params.primer_site_window}"}
      "Fixed 3' bases": {value: "${params.obipcr_fixed_3prime}"}
      "Dereplication": {value: "${params.dereplicate_amplicons}"}
      "Masking": {value: "${params.mask}"}
      ${read_length}
      ${db_filter_settings}
    EOF
    """
}
