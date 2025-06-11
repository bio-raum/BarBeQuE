process CRABS_AMPLIFICATION_EFFICENCY_FIGURE {

    tag "${meta.primer}|${meta.db}"
    
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swordfish/crabs:1.7.7.0"

    input:
    tuple val(meta), path(txt), path(db)
    val(taxon)

    output:
    tuple val(meta), path('*.png'), emit: png
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: txt.getBaseName() + "_primer_efficency"

    """
    crabs --amplification-efficiency-figure $args \\
    --input $db \\
    --amplicons $txt \\
    --tax-group $taxon \\
    --forward ${meta.fwd} \\
    --reverse ${meta.rev} \\
    --output ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
