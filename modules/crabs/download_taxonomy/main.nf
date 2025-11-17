process CRABS_DOWNLOADTAXONOMY {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.12.0--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.12.0--pyhdfd78af_0' }"

    input:
    val(meta)

    output:
    tuple path('names.dmp'), path('nodes.dmp'), path('nucl_gb.accession2taxid'), emit: taxonomy
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''

    """
    
    crabs $args --download-taxonomy --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
