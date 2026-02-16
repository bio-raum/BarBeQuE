process HELPER_CONSENSUS_DISTRIBUTION {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tsv), path(db)

    output:
    tuple val(meta), path('*.json')             , emit: json
    path 'versions.yml'                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    barbeque_size_distro.py \
    --input $tsv \
    --database $db \
    --output ${prefix}_distribution_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
