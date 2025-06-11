process CRABS_IMPORT {

    tag "${meta.db}"
    
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swordfish/crabs:1.7.7.0"

    input:
    tuple val(meta), path(fasta)
    tuple path(names), path(nodes), path(accs)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.db

    """
    crabs --import $args \\
    --input $fasta \\
    --names $names \\
    --nodes $nodes \\
    --acc2tax $accs \\
    --ranks 'superkingdom;phylum;class;order;family;genus;species' \\
    --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
