process HELPER_CLUSTER_CONSENSUS {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container 'gregdenay/taxidtools:3.1.0'

    input:
    tuple val(meta), path(clusters), path(filtered)
    path(taxdump)

    output:
    tuple val(meta), path('*.consensus.json')   , emit: json
    path 'versions.yml'                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    cluster_consensus.py \
    --clusters $clusters \
    --table $filtered \
    --taxdump $taxdump \
    --output ${prefix}.consensus.txt

    sed -i '1i SeqID\tSeq_name\tTaxid\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tlca_name\tlca_taxid\tlca_rank\tcluster_members\tamplicon' ${prefix}.consensus.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
