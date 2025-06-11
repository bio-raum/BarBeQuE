process CRABS_DOWNLOADTAXONOMY {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swordfish/crabs:1.7.7.0"

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
