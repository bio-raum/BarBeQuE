process VSEARCH_DEREPLICATION {

    tag "${meta.primer}|${meta.db}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.30.6--h0bb26bb_0' :
        'quay.io/biocontainers/vsearch:2.30.6--h0bb26bb_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.derep.fasta"), emit: fasta
    tuple val(meta), path("*.derep.uc"), emit: uc
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    vsearch \\
    --fastx_uniques ${fasta} \\
    --sizeout \\
    --uc ${prefix}.derep.uc \\
    --fastaout ${prefix}.derep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//;s/,.*//')
    END_VERSIONS
    """
}
