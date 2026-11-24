process TAXONOMIC_COVERAGE_REPORT_BLAST {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(with_region)

    output:
    tuple val(meta), path("*.coverage_blast.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.primer}_${meta.db}"
    """
    set -euo pipefail

    # Count distinct core_nt accessions returned by the restricted BLAST search.
    hits=\$(sort -u ${with_region} | wc -l)

    printf 'Sample\\tcore_nt taxid\\tDistinct core_nt accessions with region\\n' > ${prefix}.coverage_blast.tsv
    printf '%s\\t%s\\t%s\\n' \\
        '${meta.primer} (${meta.db})' \\
        '${params.taxonomic_coverage_taxid}' \\
        "\$hits" >> ${prefix}.coverage_blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n 1 | sed 's/.* //')
    END_VERSIONS
    """
}
