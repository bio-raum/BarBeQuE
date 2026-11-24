process OBIPCR_REPORT {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(tsv)

    output:
    // The recovery chart is always written; other plots require at least one amplicon.
    tuple val(meta), path("*_obipcr_*_mqc.json"), optional: true, emit: mqc
    path "versions.yml", emit: versions

    script:
    // primer and db are passed separately so the report can show the database name
    """
    obipcr_report.py \
        ${tsv} \
        "${meta.primer}" \
        "${meta.db}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
