process VSEARCH_CLUSTER_FAST {
    tag "${meta.primer}|${meta.db}"


    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vsearch:2.30.6--h0bb26bb_0'
        : 'quay.io/biocontainers/vsearch:2.30.6--h0bb26bb_0'}"

    input:
    tuple val(meta), path(fa)

    output:
    tuple val(meta), path("*.cluster.uc"), emit: uc
    tuple val(meta), path("*.centroids.fasta"), emit: centroids
    path ("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    vsearch --cluster_fast ${fa} \
    --threads ${task.cpus} \
    ${args} \
    --centroids ${prefix}.centroids.fasta \
    --uc ${prefix}.cluster.uc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
