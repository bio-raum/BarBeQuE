process CLUSTER_CONSENSUS {

    tag "${meta.primer}|${meta.db}"
    label 'medium_serial'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_taxonomy)
    path taxdump

    output:
    tuple val(meta), path("*.cluster_consensus.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    cluster_consensus.py \\
        --input ${cluster_taxonomy} \\
        --taxdump ${taxdump} \\
        --output ${prefix}.cluster_consensus.tsv ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        taxidtools: \$(python3 -c "import taxidTools; print(taxidTools.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
