process FILTER_OBIPCR {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_fasta), path(primer_input)

    output:
    tuple val(meta), path('*_filtered_raw.fasta'), emit: raw_fasta
    tuple val(meta), path('*_filtered_insilico.fasta'), emit: fasta
    tuple val(meta), path('*_filter_stats.tsv'), emit: stats
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    set -euo pipefail

    filter_obipcr_simple.py \\
        ${args} \\
        --amplicons ${raw_fasta} \\
        --primer-input ${primer_input} \\
        --primer "${meta.primer}" \\
        --db "${meta.db}" \\
        --mismatches ${params.obipcr_mismatches} \\
        --out "${prefix}_filtered_raw.fasta" \\
        --stats "${prefix}_filter_stats.tsv"

    # The same transformation OBIPCR_INSILICOPCR applies to its own raw output,
    # so filtered amplicons reach the rest of the pipeline in the usual shape.
    seqkit replace -p '_sub\\[.*' -r '' "${prefix}_filtered_raw.fasta" \\
        | seqkit seq -w 0 -u \\
        > "${prefix}_filtered_insilico.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
