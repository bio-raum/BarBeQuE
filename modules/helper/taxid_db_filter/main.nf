process TAXID_DB_FILTER {

    tag "${meta.id}|${taxid}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)
    path taxdump
    path accession_taxonomy
    val taxid

    output:
    tuple val(meta), path("*.taxid_filtered.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # Expand the requested taxid to itself + every descendant taxid
    taxonkit list --ids ${taxid} --data-dir ${taxdump} \\
        | sed 's/^ *//' \\
        | awk '{print \$1}' \\
        | grep -E '^[0-9]+\$' \\
        > taxids.txt

    taxid_db_filter.py \\
        --fasta ${db} \\
        --accession-taxid ${accession_taxonomy} \\
        --taxids taxids.txt \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version | sed 's/taxonkit v//')
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
