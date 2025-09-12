process CRABS_DIVERSITY_FIGURE {

    tag "${meta.primer}|${meta.db}"
    
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.9.0--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.9.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path('*.png'), emit: png
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: txt.getBaseName() +"_diversity"

    """
    crabs --diversity-figure $args \\
    --input $txt \\
    --output ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
