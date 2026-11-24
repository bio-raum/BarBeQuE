process JOIN_ACCESSION_TAXONOMY {

    //joins the acessions extracted from vsearch cluster .uc files using PARSE_UC module
    //with the per-db accession->taxid table produced once by BUILD_DB_TAXIDS
    //(rather than re-scanning the full genbank2taxid/accession2taxid file for every primer)
    //

    tag "${meta.primer}|${meta.db}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    tuple val(meta), path(cluster_accessions), path(master_lookup)

    output:
    tuple val(meta), path("*.cluster_taxonomy.tsv"), emit: tsv
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
set -euo pipefail

awk -F'\\t' -v prefix="${prefix}" '
    BEGIN { OFS = "\\t" }

    # First file: load master_lookup (accession -> taxid) into hash
    NR == FNR {
        taxid = (NF >= 3) ? \$3 : \$2
        lookup[\$1] = taxid
        next
    }

    # Second file: cluster_accessions (cluster_id, accession)
    {
        cluster_id = \$1
        acc = \$2
        
        if (acc in lookup) {
            print cluster_id, acc, lookup[acc]
        } else {
            print cluster_id "\\t" acc > (prefix ".missing_accessions.tsv")
        }
    }
' ${master_lookup} ${cluster_accessions} \\
> "${prefix}.cluster_taxonomy.tsv"

# Ensure missing file exists even if empty
touch "${prefix}.missing_accessions.tsv"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
END_VERSIONS
"""
}
