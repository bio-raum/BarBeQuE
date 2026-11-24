process TAXONOMIC_COVERAGE_BLAST {

    tag "${meta.primer}|${meta.db}"
    label 'process_blast'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(centroids), val(core_nt_db), val(taxid)

    output:
    tuple val(meta), path("*.with_region.txt"), emit: accessions
    tuple val(meta), path("*.blast.tsv"), emit: blast
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.primer}_${meta.db}"
    """
    set -euo pipefail

    # Start the BLAST table with a header line
    # One row per hit: which centroid hit which record, how well, and where in the record
    # record_from > record_to means the region is on the reverse strand

    echo -e "centroid\\trecord\\tidentity\\tcentroid_cov\\taligned\\tmismatches\\tgap_opens\\tcentroid_len\\trecord_len\\trecord_from\\trecord_to\\tstrand\\tevalue\\tbitscore" > ${prefix}.blast.tsv

    # The taxid limits the core_nt search to that taxon and its descendants.
    blastdbcmd -db "${core_nt_db}" -info > /dev/null

    blastn -query ${centroids} -db "${core_nt_db}" \\
        -taxids ${taxid} \\
        -perc_identity ${params.taxonomic_coverage_identity} \\
        -qcov_hsp_perc ${params.taxonomic_coverage_query_cov} \\
        -max_target_seqs 1000000 \\
        -dust no \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen qlen slen sstart send sstrand evalue bitscore" >> ${prefix}.blast.tsv

    # Keep each hit accession once, even if several centroids match it.

    tail -n +2 ${prefix}.blast.tsv | cut -f2 | cut -d';' -f1 | sort -u > ${prefix}.with_region.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n 1 | sed 's/blastn: //')
    END_VERSIONS
    """
}
