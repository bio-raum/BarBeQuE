#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/**
 * Run the normal BarBeQuE analysis independently for every target taxid in
 * params.teeliste. Each target includes the requested taxid and all of its
 * descendants. Empty target/database combinations are recorded in an audit
 * table but are not passed to the main workflow.
 */

params.version = workflow.manifest.version

include { BARBEQUE } from './workflows/barbeque'
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_completion'
include { DATABASE } from './subworkflows/database'
include { INTERACTIVE_RESULTS } from './workflows/interactive_results'
include { paramsHelp ; paramsSummaryLog ; validateParameters } from 'plugin/nf-schema'


def pathToken(value) {
    def token = value
        .toString()
        .trim()
        .replaceAll(/[^A-Za-z0-9]+/, '_')
        .replaceAll(/^_+|_+$/, '')
    return token ?: 'target'
}


def shellQuote(value) {
    return "'${value.toString().replace("'", "'\"'\"'")}'"
}


process TARGET_TAXID_DB_FILTER {

    tag "${meta.source_db}|${target.taxid}"
    label 'medium_parallel'

    conda "${moduleDir}/modules/helper/taxid_db_filter/environment.yml"

    input:
    tuple val(meta), path(db), val(target)
    path taxdump
    path accession_taxonomy

    output:
    tuple val(meta), path("${meta.id}.taxid_filtered.fasta"), optional: true, emit: fasta
    tuple val(meta), path("${meta.id}.target_status.tsv"), emit: audit
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    # Expand the target taxid to itself and every descendant taxid.
    { \
        printf '%s\\n' ${target.taxid}; \
        taxonkit list --ids ${target.taxid} --data-dir ${taxdump}; \
    } \
        | sed 's/^ *//' \
        | awk '\$1 ~ /^[0-9]+\$/ {print \$1}' \
        | sort -u \
        > taxids.txt

    taxid_db_filter.py \
        --fasta ${db} \
        --accession-taxid ${accession_taxonomy} \
        --taxids taxids.txt \
        --prefix ${prefix}

    filtered=${prefix}.taxid_filtered.fasta
    matched=\$(grep -c '^>' "\$filtered" || true)
    status=matched
    if [[ "\$matched" -eq 0 ]]; then
        status=no_matches
        rm "\$filtered"
    fi

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        source_database target_taxid latin_name german_name target_rank matched_sequence_count status \
        > ${meta.id}.target_status.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        ${shellQuote(meta.source_db)} \
        ${shellQuote(target.taxid)} \
        ${shellQuote(target.latin_name)} \
        ${shellQuote(target.german_name)} \
        ${shellQuote(target.rank)} \
        "\$matched" \
        "\$status" \
        >> ${meta.id}.target_status.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version | sed 's/taxonkit v//')
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}


process TARGET_FILTER_AUDIT {

    tag 'all target/database combinations'
    label 'short_serial'

    publishDir "${params.outdir}/target_analysis", mode: params.publish_dir_mode, overwrite: true

    input:
    path audit_files

    output:
    path 'target_filter_audit.tsv'

    script:
    """
    printf 'source_database\\ttarget_taxid\\tlatin_name\\tgerman_name\\ttarget_rank\\tmatched_sequence_count\\tstatus\\n' \
        > target_filter_audit.tsv

    for audit in ${audit_files}; do
        tail -n +2 "\$audit"
    done \
        | sort -t \$'\\t' -k1,1 -k5,5 -k3,3 -k2,2n \
        >> target_filter_audit.tsv
    """
}


workflow {

    if (!workflow.containerEngine) {
        log.info("\033[1;31mRunning with Conda is not recommended in production!\033[0m\n\033[0;31mConda environments are not guaranteed to be reproducible - for a discussion, see https://pubmed.ncbi.nlm.nih.gov/29953862/.\033[0m")
    }

    if (WorkflowPipeline.enabled(params.help)) {
        log.info(paramsHelp(command: 'nextflow run analysis.nf'))
        System.exit(0)
    }
    if (WorkflowPipeline.enabled(params.helpFull)) {
        System.exit(0)
    }

    validateParameters()
    WorkflowMain.initialise(workflow, params, log)

    if (params.taxid) {
        log.error('analysis.nf reads all target taxids from --teeliste; do not also provide --taxid')
        System.exit(1)
    }
    if (WorkflowPipeline.enabled(params.build_references)) {
        log.error('analysis.nf analyses the targets in --teeliste; use main.nf to build references')
        System.exit(1)
    }

    WorkflowPipeline.initialise(params, log)
    log.info(paramsSummaryLog(workflow))

    ch_targets = channel.fromPath(params.teeliste, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            def taxid = row.taxid?.toString()?.trim()
            def latinName = row['lat.']?.toString()?.trim()
            def germanName = row['German name']?.toString()?.trim()
            def rank = row.taxonomic_rank?.toString()?.trim()

            if (!(taxid ==~ /[0-9]+/)) {
                throw new IllegalArgumentException("Invalid taxid '${taxid}' in ${params.teeliste}")
            }
            if (!latinName) {
                throw new IllegalArgumentException("Missing Latin name for taxid ${taxid} in ${params.teeliste}")
            }

            [
                taxid: taxid,
                latin_name: latinName,
                german_name: germanName ?: '',
                rank: rank ?: '',
                slug: pathToken(latinName),
            ]
        }
        .distinct { target -> target.taxid }

    DATABASE()

    DATABASE.out.db
        .combine(ch_targets)
        .map { dbMeta, db, target ->
            def sourceDb = dbMeta.id.toString()
            def targetMeta = new LinkedHashMap(dbMeta)
            targetMeta.id = "${pathToken(sourceDb)}__taxid_${target.taxid}_${target.slug}"
            targetMeta.source_db = sourceDb
            targetMeta.target_taxid = target.taxid
            targetMeta.target_latin_name = target.latin_name
            targetMeta.target_german_name = target.german_name
            targetMeta.target_rank = target.rank
            tuple(targetMeta, db, target)
        }
        .set { ch_target_database }

    TARGET_TAXID_DB_FILTER(
        ch_target_database,
        DATABASE.out.taxdump,
        DATABASE.out.accession_taxonomy,
    )

    TARGET_FILTER_AUDIT(
        TARGET_TAXID_DB_FILTER.out.audit.map { _meta, audit -> audit }.collect()
    )

    ch_versions = DATABASE.out.versions.mix(TARGET_TAXID_DB_FILTER.out.versions)

    BARBEQUE(
        TARGET_TAXID_DB_FILTER.out.fasta,
        ch_versions,
        DATABASE.out.taxdump,
        DATABASE.out.accession_taxonomy,
    )

    if (WorkflowPipeline.enabled(params.interactive)) {
        BARBEQUE.out.consensus.collect() | map { "${params.outdir}" } | INTERACTIVE_RESULTS
    }

    PIPELINE_COMPLETION()
}
