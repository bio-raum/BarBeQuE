process DB_FILTER {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0' :
        'quay.io/biocontainers/seqkit:2.13.0--he881be0_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.cleaned.fasta"), emit: fasta
    path("versions.yml"),                       emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pattern = params.db_filter_pattern

    // The two length bounds are independent - seqkit takes --min-len and --max-len
    // separately, so leaving db_filter_max_length unset means "no upper limit"
    // rather than disabling length filtering altogether.
    def length_bounds = []
    if (params.db_filter_min_length != null) {
        length_bounds << "--min-len ${params.db_filter_min_length}"
    }
    if (params.db_filter_max_length != null) {
        length_bounds << "--max-len ${params.db_filter_max_length}"
    }
    def length_filter = length_bounds
        ? "| seqkit seq ${length_bounds.join(' ')} --threads ${task.cpus}"
        : ""
    // SeqKit has no --max-ambig option. Reject records containing more than
    // the allowed number of IUPAC ambiguity symbols with a sequence regex.
    // Independent of the length bounds above.
    def ambiguity_filter = params.db_filter_max_n != null
        ? "| seqkit grep --by-seq --use-regexp --ignore-case --only-positive-strand --invert-match" +
          " --pattern \"[RYKMSWBDHVN]([^RYKMSWBDHVN]*[RYKMSWBDHVN]){${params.db_filter_max_n}}\"" +
          " --threads ${task.cpus}"
        : ""

    """
    seqkit grep -n -v -r -i \
        --threads ${task.cpus} \
        -p "${pattern}" \
        ${db} \
        ${length_filter} \
        ${ambiguity_filter} \
        -o ${prefix}.cleaned.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
