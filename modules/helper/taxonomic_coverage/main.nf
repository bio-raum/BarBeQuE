process TAXONOMIC_COVERAGE {

    //provides taxonomic coverage using ETE toolkit
    //input files: vsearch clusters, db taxids file
    //taxonomy is optional 
    cache false

    tag "${meta.primer}|${meta.db}"
    label 'medium_parallel'


    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(clusters), path(db_taxids)
    val taxonomy

    output:
    tuple val(meta), path('*.tax_coverage.tsv'), emit: tsv
    tuple val(meta), path('*.tax_coverage.nwk'), emit: nwk
    tuple val(meta), path('*.tax_coverage_summary_mqc.tsv'), emit: mqc_summary
    tuple val(meta), path('*.tax_coverage_details_mqc.tsv'), emit: mqc_details
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    ete.py --taxon "${taxonomy}" \\
        --reference "${db_taxids}" \\
        --report "${clusters}" \\
        --output "${prefix}--${taxonomy}--.tax_coverage" \\
        --sample "${meta.primer} (${meta.db})"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
