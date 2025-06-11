//
// Check input samplesheet and get read channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row -> fastq_channel(row) }
        .set { primers }

    emit:
    primers 
}

// Function to get meta hash
def fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.primer    = row.primer

    if (row.FWD) {
        meta.fwd       = row.FWD
    } else {
        log.info "No forward primer defined, cannot proceed!\n"
        System.exit(1)
    }
    if (row.REV) {
        meta.rev       = row.REV
    } else {
        log.info "No reverse primer defined, cannot proceed!\n"
        System.exit(1)
    }
    
    meta.min  = row.min
    meta.max  = row.max

    return meta
}
