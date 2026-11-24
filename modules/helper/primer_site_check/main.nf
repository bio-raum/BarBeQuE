// EXPERIMENTAL MODULE
// Classifies why a primer did not bind in records that have the barcode but gave no amplicon:
// PRIMER_CUTOFF / PRIMER_SITE_NOT_FOUND / UNSEQUENCED / MISMATCH / BOUND.
// Tested on test_0.9 (trnL, psbA-trnH, 83 taxids) only. Not handled yet: barcodes that wrap around
// the end of a circular genome, insertions/deletions in the primer site (counted as several
// mismatches), the 3' clamp (obipcr_fixed_3prime). Check the results by hand before relying on them.

process PRIMER_SITE_CHECK {

    tag "${meta.primer}|${fasta.baseName}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(blast), path(obipcr_tsv)

    output:
    path "*.primer_sites.tsv", emit: tsv
    tuple val(meta), path("*_mismatch_sites.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    // Example: trnL_4442_Camellia_sinensis_species.with_barcode_not_amplified.fasta -> trnL_4442_Camellia_sinensis_species
    def prefix = fasta.name.replace('.with_barcode_not_amplified.fasta', '')
    """
    # For each record with the barcode but no amplicon: why did each primer not bind?

    primer_site_check.py \\
        --fasta ${fasta} \\
        --blast ${blast} \\
        --obipcr-tsv ${obipcr_tsv} \\
        --mismatches ${params.obipcr_mismatches} \\
        --site-not-found-mismatches ${params.primer_site_not_found_mismatches} \\
        --window ${params.primer_site_window} \\
        --name ${prefix} \\
        --output ${prefix}.primer_sites.tsv \\
        --fasta-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
