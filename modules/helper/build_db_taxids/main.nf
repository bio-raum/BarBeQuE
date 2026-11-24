process BUILD_DB_TAXIDS {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"


    input:
    tuple val(meta), path(fasta)
    path genbank2taxid

    output:
    tuple val(meta), path("*.db_taxids.tsv"), emit: taxids
    tuple val(meta), path("*.db_taxids_counts.tsv"), emit: taxids_counts
    tuple val(meta), path("*.accession_taxid.tsv"), emit: accession_taxid
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    tuple val(meta), path("*.omitted_accessions_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # Extract accession IDs from the FASTA headers.
    # SINTAX annotations after the first semicolon are removed.
    grep "^>" "${fasta}" | cut -d' ' -f1 | cut -d'>' -f2 | cut -d';' -f1 > db_accessions.txt

    # Match database accessions against the accession-to-taxid mapping.
    awk -F'\\t' -v acc_taxid_out="${prefix}.accession_taxid.tsv" '
        # Load database accessions and remove accession version suffixes.
        NR == FNR {
            sub(/(\\.[0-9]+)+\$/, "", \$1)
            db_accessions[\$1] = 1
            next
        }

        # Write accession-to-taxid matches and output taxids for counting.
        {
            accession = \$1
            sub(/(\\.[0-9]+)+\$/, "", accession)
            if (accession in db_accessions) {
                taxid = (NF >= 3) ? \$3 : \$2
                print accession"\\t"taxid > acc_taxid_out
                print taxid
                delete db_accessions[accession]
            }
        }

        # Record accessions that could not be assigned to a taxid.
        END {
            for (unmatched in db_accessions) {
                print unmatched > "${prefix}.missing_accessions.tsv"
            }
        }
    ' db_accessions.txt "${genbank2taxid}" \\
    | sort -n \\
    | tee >(sort -u > "${prefix}.db_taxids.tsv") \\
    | uniq -c \\
    | awk '{print \$2"\\t"\$1}' \\
    > "${prefix}.db_taxids_counts.tsv"

    rm db_accessions.txt

    # Ensure optional output files exist even when empty.
    touch "${prefix}.missing_accessions.tsv"
    touch "${prefix}.accession_taxid.tsv"

    # Create a MultiQC summary of accessions without taxid assignments.
    printf "# id: 'omitted_accessions'\\n# parent_id: 'build_db_taxids'\\n# parent_name: 'Database Taxonomy Lookup'\\n# section_name: 'Omitted accessions'\\n# plot_type: 'table'\\nDatabase\\tOmitted accessions\\n%s\\t%s\\n" "${meta.id}" "\$(wc -l < ${prefix}.missing_accessions.tsv)" > "${prefix}.omitted_accessions_mqc.tsv"

    # Record the awk version used.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
