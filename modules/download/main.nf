process DOWNLOAD_MIDORI {
    tag 'MIDORI'
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    val url

    output:
    path 'midori.zip', optional: true, emit: archive
    stdout emit: warnings

    script:
    """
    set -euo pipefail

    if ! wget --tries=3 --timeout=60 -O midori.zip '${url}'; then
        echo "WARNING: MIDORI archive download failed: ${url}"
        rm -f midori.zip
    fi
    """
}

process EXTRACT_MIDORI {
    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${meta.db}.fasta"), optional: true, emit: db
    stdout emit: warnings

    script:
    """
    set -euo pipefail

    # The ZIP member is itself a gzip-compressed FASTA.
    if unzip -p '${archive}' '${meta.member}' | gzip -cd > '${meta.db}.fasta.partial' && grep -q '^>' '${meta.db}.fasta.partial'; then
        mv '${meta.db}.fasta.partial' '${meta.db}.fasta'
    else
        echo "WARNING: MIDORI ${meta.db}: ${meta.member} could not be extracted as FASTA. Skipping."
        rm -f '${meta.db}.fasta.partial'
    fi
    """
}
