process PARSE_UC {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    tuple val(meta), path(uc)

    output:
    tuple val(meta), path("*.cluster_accessions.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    set -euo pipefail

    # Parse VSEARCH .uc output and keep seed (S) and hit (H) records.
    # Column 9 contains the FASTA sequence label.
    # Accessions are normalised to match BUILD_DB_TAXIDS:
    #   AB189069.1                       -> AB189069
    #   FM163243.1;tax=k:Fungi;          -> FM163243
    #   AY846379.1.1791                  -> AY846379
    awk -F'\\t' '
        BEGIN { OFS = "\\t" }

        \$1 == "S" || \$1 == "H" {

            # Extract the sequence accession.
            acc = \$9

            # Remove annotations after the first semicolon.
            sub(/;.*/, "", acc)

            # Remove accession version/range suffixes.
            sub(/(\\.[0-9]+)+\$/, "", acc)

            # Output cluster ID and cleaned accession.
            if (acc != "*" && acc != "") {
                print \$2, acc
            }
        }
    ' "${uc}" \\
    | sort -k1,1n -k2,2 \\
    > "${prefix}.cluster_accessions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
