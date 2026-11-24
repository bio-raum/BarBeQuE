// obipcr forbids a mismatch at a base when that base is followed by '#', so clamping
// the last few bases reproduces a real PCR: a mismatch at the primer's 3' end stops the
// polymerase from extending, while mismatches further 5' are still tolerated.
// GGGCAATCCTGAGCCAA -> GGGCAATCCTGAGCC#A#A#
def clamp_3prime(primer, positions) {
    def bases = primer.replace('#', '')
    if (positions <= 0 || bases.length() <= positions) {
        return bases
    }
    return bases[0..<(bases.length() - positions)] + bases[-positions..-1].collect { base -> base + '#' }.join()
}

process OBIPCR_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*_insilico.fasta'), emit: fasta
    tuple val(meta), path('*_raw.fasta'), emit: raw_fasta, optional: true
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    def mismatches = params.obipcr_mismatches as int

    def fixed_3prime = params.containsKey('obipcr_fixed_3prime')
        ? params.obipcr_fixed_3prime as int
        : 0

    def fwd_primer = clamp_3prime(meta.fwd, fixed_3prime)
    def rev_primer = clamp_3prime(meta.rev, fixed_3prime)
    def circular_arg = WorkflowPipeline.enabled(params.obipcr_circular) ? '--circular' : ''

    """
    obipcr ${args} \\
        --forward "${fwd_primer}" \\
        --reverse "${rev_primer}" \\
        --allowed-mismatches ${mismatches} \\
        --min-length ${meta.min} \\
        --max-length ${meta.max} \\
        --max-cpu ${task.cpus} \\
        --fasta-output \\
        --no-progressbar \\
        --skip-empty \\
        ${circular_arg} \\
        "${db}" \\
        | tee "${prefix}_raw.fasta" \\
        | seqkit replace -p '_sub\\[.*' -r '' \\
        | seqkit seq -w 0 -u \\
        > "${prefix}_insilico.fasta"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    obipcr: "\$(obipcr --version 2>&1 | head -n 1)"
END_VERSIONS
    """
}
