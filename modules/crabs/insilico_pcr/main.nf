process CRABS_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.14.0--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.14.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*insilico.txt'), emit: txt
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.primer + "_" + meta.db

    """
    crabs $args --in-silico-pcr \\
    --input $db \\
    --output ${prefix}_insilico.txt \\
    --forward ${meta.fwd} \\
    --reverse ${meta.rev} \\
    --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
