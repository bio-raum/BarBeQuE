process CLUSTER_CONSENSUS_REPORTING {

    tag "${meta.primer}|${meta.db}"
    label 'medium_serial'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(cluster_consensus)

    output:
    tuple val(meta), path("*.taxonomic_resolution_mqc.tsv"), emit: resolution
    tuple val(meta), path("*.overall_summary_mqc.tsv"), emit: summary
    tuple val(meta), path("*.cluster_size_mqc.tsv"), emit: cluster_size
    tuple val(meta), path("*.taxon_search_mqc.html"), emit: taxon_search
    tuple val(meta), path("*.amplified_taxa_mqc.tsv"), emit: amplified_taxa
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    def args = task.ext.args ?: ''

    """
    set -euo pipefail
    # Generate the current custom Taxon Search report content.

    cluster_consensus_reporting_self.py \
        --input "${cluster_consensus}" \
        --output-prefix "${prefix}" \
        --sample "${meta.primer} (${meta.db})" \
        ${args}

    mkdir taxon_search
    (
        cd taxon_search
        python3 \$(command -v cluster_consensus_reporting.py) "../${cluster_consensus}"
    )
    mv taxon_search/taxon_search_mqc.html "${prefix}.taxon_search_mqc.html"
    mv taxon_search/amplified_taxa_mqc.tsv "${prefix}.amplified_taxa_mqc.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
