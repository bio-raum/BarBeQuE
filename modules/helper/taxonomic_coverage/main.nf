process HELPER_TAXONOMIC_COVERAGE {

    maxForks 1

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_1' :
        'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(clusters), path(db)
    val(taxonomy)

    output:
    tuple val(meta), path('*.tsv') , emit: tsv
    tuple val(meta), path('*.nwk') , emit: nwk
    path 'versions.yml'            , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    cut -f3 $db | sort -n -u > ids.txt

    ete.py --taxon $taxonomy \\
    --reference ids.txt \\
    --report $clusters \\
    --output ${prefix}.tax_coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
