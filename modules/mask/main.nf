process MASK {

    tag "${meta.primer}|${meta.db}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*_masked.fasta'), emit: fasta
    path ('versions.yml'), emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
mask.py \
    --input ${fasta} \
    --output "${prefix}_masked.fasta" \
    --read-length ${params.read_length} \
    ${WorkflowPipeline.enabled(params.single_end) ? '--single-end' : ''}

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
