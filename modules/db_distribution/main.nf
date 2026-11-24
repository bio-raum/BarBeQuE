process DB_DISTRIBUTION {

    tag "${meta.id}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"



    input:
    tuple val(meta), path(taxids_counts)
    path taxdump

    output:
    tuple val(meta), path("*.db_distribution.tsv"), emit: distribution
    tuple val(meta), path("*.db_distribution.html"), emit: plot
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_distribution.py \
        --input ${taxids_counts} \
        --output ${prefix}.db_distribution.tsv \
        --taxdump ${taxdump}

    db_distribution_plot.py \
        --input  ${prefix}.db_distribution.tsv \
        --output ${prefix}.db_distribution.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        taxidTools: \$(python -c "import taxidTools; print(taxidTools.__version__)")
    END_VERSIONS
    """
}
