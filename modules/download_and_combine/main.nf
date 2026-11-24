process DOWNLOAD_AND_COMBINE_DB {
    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(urls)

    output:
    tuple val(meta), path("${meta.filename}"), emit: db, optional: true
    path 'versions.yml', emit: versions
    stdout emit: warnings

    script:
    if (meta.format != 'fasta' || !urls) {
        error "${meta.db} requires at least one FASTA URL"
    }
    def links = urls.collect { "'${it}'" }.join(' ')

    """
    set -euo pipefail

    # Download all parts before creating the final database.
    : > "${meta.filename}.partial"
    failed=0

    for url in ${links}; do
        if ! wget --tries=3 --timeout=60 -O downloaded "\$url"; then
            echo "${meta.db}: Could not download \$url. Skipping database."
            failed=1
            break
        fi

        # Decompress gzip files or copy plain FASTA files.
        if gzip -t downloaded 2>/dev/null; then
            gzip -cd downloaded > part.fasta
        else
            cp downloaded part.fasta
        fi
        # Check whether the file contains FASTA sequences.
        if ! grep -q '^>' part.fasta; then
            echo "${meta.db}: Invalid FASTA from \$url. Skipping database."
            failed=1
            break
        fi
        cat part.fasta >> "${meta.filename}.partial"
    done

    # Publish the database only when all parts succeeded.
    if [ "\$failed" -eq 0 ]; then
        mv "${meta.filename}.partial" "${meta.filename}"
    else
        rm -f "${meta.filename}.partial"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        source_release: '${meta.release}'
        wget: \$(wget --version | head -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
