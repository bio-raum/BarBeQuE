process CRABS_FILTER {

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.14.0--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.14.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dereplicate)

    output:
    tuple val(meta), path('*filtered.txt')  , emit: txt
    tuple val(meta), path('*filtered.fa')   , emit: fasta
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.args ?: meta.primer + "_" + meta.db

    """
    crabs $args --filter \\
    --input $dereplicate \\
    --output ${prefix}_filtered.txt \\
    --minimum-length $meta.min \\
    --maximum-length $meta.max \\
    --maximum-n 1 \\
    --environmental \\
    --no-species-id \\
    --rank-na 2 \\
    --threads ${task.cpus}

    awk -F "\t" '{print ">"\$1"\\n"\$11}' ${prefix}_filtered.txt > ${prefix}_filtered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs -version 2>&1 | grep "This is CRABS" | sed -e "s/This is CRABS //g")
    END_VERSIONS

    """
}
