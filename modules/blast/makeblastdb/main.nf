process MAKEBLASTDB {

    tag "${fasta.baseName}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta

    output:
    tuple val(fasta.baseName), path(fasta), path("db"), emit: db
    path "versions.yml", emit: versions

    script:
    """
    # Make a BLAST database from the taxid records
    # Example: 4442_Camellia_sinensis_species.fasta -> db/4442_Camellia_sinensis_species
    # An empty taxid file is skipped (BLAST cannot use it), db/ then stays empty

    mkdir db

    if [ -s ${fasta} ]; then
        makeblastdb -in ${fasta} -dbtype nucl -out db/${fasta.baseName}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$(makeblastdb -version | head -n 1 | sed 's/makeblastdb: //')
    END_VERSIONS
    """
}
