process CRABS_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swordfish/crabs:1.7.7.0"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*insilico.txt'), emit: txt
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.args ?: meta.primer + "_" + meta.db

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
