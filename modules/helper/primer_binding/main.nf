process PRIMER_BINDING {

    tag "${meta.primer}|${fasta.baseName}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(coverage), path(obipcr_tsv)

    output:
    path "*.seqkit_locate.tsv", emit: locate
    path "*.coverage_binding.tsv", emit: coverage
    path "versions.yml", emit: versions

    script:
    // Example: trnL_4442_Camellia_sinensis_species.with_barcode_not_amplified.fasta -> trnL_4442_Camellia_sinensis_species
    def prefix = fasta.name.replace('.with_barcode_not_amplified.fasta', '')
    """
    # Use the actual original primers selected by FILTER_OBIPCR and recorded by
    # PARSE_OBIPCR, rather than the merged degenerate primers in meta.fwd/meta.rev.

    tail -n +2 ${obipcr_tsv} | cut -f3 | sort -u \\
        | awk '{print ">forward_" NR "\\n" toupper(\$0)}' > primers.fasta
    tail -n +2 ${obipcr_tsv} | cut -f9 | sort -u \\
        | awk '{print ">reverse_" NR "\\n" toupper(\$0)}' >> primers.fasta

    # Find both primers in the records with the barcode but no amplicon,
    # with as many mismatches as obipcr allows
    # Each line: record, primer, strand, start, end, matched bases

    seqkit locate \\
        --max-mismatch ${params.obipcr_mismatches} \\
        --threads ${task.cpus} \\
        --pattern-file primers.fasta \\
        ${fasta} \\
        --out-file ${prefix}.seqkit_locate.tsv

    # Count the records where each primer bound, each record once
    # awk keeps the lines of one primer, cut -f1 keeps the record column

    not_amplified=\$(grep -c '>' ${fasta} || true)
    forward_bound=\$(awk -F'\\t' '\$2 ~ /^forward_/' ${prefix}.seqkit_locate.tsv | cut -f1 | sort -u | wc -l)
    reverse_bound=\$(awk -F'\\t' '\$2 ~ /^reverse_/' ${prefix}.seqkit_locate.tsv | cut -f1 | sort -u | wc -l)

    # Records where a primer did NOT bind = not amplified records - records where it bound

    forward_not_bound=\$((not_amplified - forward_bound))
    reverse_not_bound=\$((not_amplified - reverse_bound))

    # The coverage row with the two new columns added at the end

    head -n 1 ${coverage} | sed 's/\$/\\tforward_not_bound\\treverse_not_bound/' > ${prefix}.coverage_binding.tsv
    tail -n +2 ${coverage} | sed "s/\$/\\t\$forward_not_bound\\t\$reverse_not_bound/" >> ${prefix}.coverage_binding.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
