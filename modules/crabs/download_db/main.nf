process CRABS_DOWNLOADDB {

    maxForks 3

    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.14.0--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.14.0--pyhdfd78af_0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path('*.fasta'), emit: db
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'crabs'
    
    """
    crabs $args --output ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
