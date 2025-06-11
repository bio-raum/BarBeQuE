process CRABS_DOWNLOADDB {

    maxForks 3

    tag "${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swordfish/crabs:1.7.7.0"

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
