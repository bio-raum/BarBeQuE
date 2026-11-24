process COVERAGE_ANALYSIS {

    tag "${meta.primer}|${name}"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(centroids), path(obipcr_tsv), val(name), path(fasta), path(db)

    output:
    path "*.coverage.tsv", emit: tsv
    path "*.blast.tsv", emit: blast
    tuple val(meta), path("*.with_region.fasta"), path("*.no_region.fasta"), path("*.with_barcode_not_amplified.fasta"), emit: split
    tuple val(meta), path("*.with_barcode_not_amplified.fasta"), path("*.coverage.tsv"), path(obipcr_tsv), emit: not_amplified
    tuple val(meta), path("*.with_barcode_not_amplified.fasta"), path("*.blast.tsv"), path(obipcr_tsv), emit: site_check
    path "versions.yml", emit: versions

    script:
    // Split the file name into its parts for the report
    // Example: 155109_Cyclopia_subternata_species -> 155109, Cyclopia subternata, species
    def parts = name.tokenize('_')
    def taxid = parts.first()
    def rank = parts.last()
    def taxon_name = parts[1..-2].join(' ')
    """
    # Record the tool versions first, because the script can stop early below

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n 1 | sed 's/blastn: //')
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS

    # Start the report with a header line

    echo -e "primer\\ttaxid\\tname\\trank\\trecords_db\\twith_region_blast\\tamplified_obipcr\\twith_barcode_not_amplified" > ${meta.primer}_${name}.coverage.tsv

    # Start the BLAST table with a header line (it stays like this for an empty taxid file)

    echo -e "centroid\\trecord\\tidentity\\tcentroid_cov\\taligned\\tcentroid_len\\trecord_from\\trecord_to\\tcentroid_from\\tcentroid_to\\tbitscore" > ${meta.primer}_${name}.blast.tsv

    # Count the records in the taxid file

    records=\$(grep -c '>' ${fasta} || true)

    # Empty taxid file: nothing to search, so the report row is all 0
    # and all FASTAs stay empty

    if [ "\$records" -eq 0 ]; then
        touch ${meta.primer}_${name}.with_region.fasta ${meta.primer}_${name}.no_region.fasta ${meta.primer}_${name}.with_barcode_not_amplified.fasta
        echo -e "${meta.primer}\\t${taxid}\\t${taxon_name}\\t${rank}\\t0\\t0\\t0\\t0" >> ${meta.primer}_${name}.coverage.tsv
        exit 0
    fi

    # Search the centroids against the taxid records
    # A record has the region if one hit covers most of a centroid at high identity
    # -max_target_seqs: BLAST stops at 500 hits per centroid by default, so allow many more
    # -dust no: do not mask simple repeats like ATATAT
    # Each line: which centroid hit which record, how well, and where in the record

    blastn -query ${centroids} -db db/${name} \\
        -perc_identity ${params.coverage_identity} \\
        -qcov_hsp_perc ${params.coverage_query_cov} \\
        -max_target_seqs 1000000 \\
        -dust no \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident qcovhsp length qlen sstart send qstart qend bitscore" >> ${meta.primer}_${name}.blast.tsv

    # List the records that have the region, each record once (even if several centroids hit it)
    # tail skips the header line, cut -f2 keeps the record column
    # BLAST drops a trailing ';' from the record name (e.g. "MZ190997;tax=...;"), so the names
    # are compared by the part before the first ';' only: cut -d';' -f1 here and
    # --id-regexp in seqkit grep below. NCBI names like "MW812458.1" stay unchanged.

    tail -n +2 ${meta.primer}_${name}.blast.tsv | cut -f2 | cut -d';' -f1 | sort -u > with_region.txt

    # Split the taxid FASTA into two files, with the full records unchanged
    # -f with_region.txt: keep the records in the list
    # -v: keep the records NOT in the list
    # --id-regexp: the record name up to the first ';' or space

    seqkit grep --id-regexp '^([^;\\s]+)' -f with_region.txt ${fasta} > ${meta.primer}_${name}.with_region.fasta
    seqkit grep --id-regexp '^([^;\\s]+)' -v -f with_region.txt ${fasta} > ${meta.primer}_${name}.no_region.fasta

    # Count the records with the region

    with_region=\$(wc -l < with_region.txt)

    # Count the taxid records that gave an amplicon
    # Records, not amplicons: a record can give more than one amplicon, so each accession once
    # tail skips the header line, cut -f1 keeps the Sequence_ID column

    tail -n +2 ${obipcr_tsv} | cut -f1 | cut -d';' -f1 | sort -u > amplified.txt
    amplified=\$(seqkit grep --id-regexp '^([^;\\s]+)' -f amplified.txt ${fasta} | grep -c '>' || true)

    # Records with the region that gave no amplicon (with_region minus amplified)
    # -v: keep the records NOT in the amplified list

    seqkit grep --id-regexp '^([^;\\s]+)' -v -f amplified.txt ${meta.primer}_${name}.with_region.fasta > ${meta.primer}_${name}.with_barcode_not_amplified.fasta
    not_amplified=\$(grep -c '>' ${meta.primer}_${name}.with_barcode_not_amplified.fasta || true)

    echo -e "${meta.primer}\\t${taxid}\\t${taxon_name}\\t${rank}\\t\$records\\t\$with_region\\t\$amplified\\t\$not_amplified" >> ${meta.primer}_${name}.coverage.tsv
    """
}
