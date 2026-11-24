process DOWNLOAD_MITOFISH {
    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path('mitofishdb.fa.gz'), optional: true, emit: gz
    stdout emit: warnings

    script:
    """
    set -euo pipefail

    # Download MitoFish; an unavailable source leaves no archive to convert.
    if ! wget --tries=3 --timeout=60 -O mitofishdb.fa.gz '${url}' || [[ ! -s mitofishdb.fa.gz ]]; then
        echo "WARNING: ${meta.db}: MitoFish download failed."
        rm -f mitofishdb.fa.gz
    fi
    """
}

process CONVERT_MITOFISH {
    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gz_file)

    output:
    tuple val(meta), path("${meta.filename}"), emit: db

    script:
    """
    set -euo pipefail

    # Decompress and validate the FASTA before creating its final filename.
    mitofish_fasta.py '${gz_file}' '${meta.filename}'
    """
}
