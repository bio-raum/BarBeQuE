//
// Check input samplesheet and get read channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv

    main:
    samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> fastq_channel(row) }
        .map { meta -> [meta.primer, meta] }
        .groupTuple()
        .map { primer, metas ->
            if (metas.size() > 1) {
                error("Duplicate primer name detected after sanitisation: '${primer}'. Primer names must be unique.")
            }
            return metas[0]
        }
        .set { primers }

    emit:
    primers
}

// Function to get meta hash
def fastq_channel(LinkedHashMap row) {
    def meta = [:]
    
    if (!row.primer) {
        error("No primer name defined (primer), cannot proceed!")
    }
    
    meta.primer = row.primer.replaceAll(/[\s\/\\:*?"<>|]/, '_')
    
    if (row.fwd) {
        meta.fwd = row.fwd
    }
    else {
        error("No forward primer defined (fwd), cannot proceed!")
    }
    if (row.rev) {
        meta.rev = row.rev
    }
    else {
        error("No reverse primer defined (rev), cannot proceed!")
    }
    if (row.min) {
        meta.min = row.min
    }
    else {
        error("Must provide minimum amplicon length for filtering (min)")
    }
    if (row.max) {
        meta.max = row.max
    }
    else {
        error("Must provide maximum amplicon length for filtering (max)")
    }

    return meta
}

