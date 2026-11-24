# Pipeline Workflow

`main.nf` validates the parameters and then runs one of two modes.

```groovy
if (build_references) {
    BUILD_REFERENCES()
} else {
    DATABASE()
    BARBEQUE(DATABASE.out.db, DATABASE.out.versions, DATABASE.out.taxdump, DATABASE.out.accession_taxonomy)
    if (interactive) INTERACTIVE_RESULTS
}
PIPELINE_COMPLETION()
```

## Benchmarking mode (default)

This mode needs `--input` or `--primer_set`, plus `--dbs` or `--custom_db`.

1. **Database:** `DATABASE` resolves the reference FASTAs and taxonomy. It then applies the
   optional `--taxid` restriction, followed by the optional `--db_filter` cleaning.
2. **Primers:** `--input` FASTAs and `--primer_set` downloads go through `PARSE_PRIMERS`. A
   samplesheet is used as-is. `INPUT_CHECK` then validates the result. See
   [primer_input.md](primer_input.md).
3. **In-silico PCR:** `obipcr` runs once for each primer/database pair.
4. **Filter:** for primer FASTA input, `FILTER_OBIPCR` drops amplicons that match only the merged
   degenerate primer. Disable it with `--filter_collapsed_primers false`.
5. **Parse:** amplicons are parsed into a TSV, with a MultiQC section. Pairs that produced no
   amplicons are dropped with a warning.
6. **Optional processing:** dereplication (`--dereplicate_amplicons`) and masking (`--mask`).
7. **Clustering:** amplicons are clustered with `vsearch --cluster_fast` (`--cluster_id`), and each
   cluster member is joined to its taxid.
8. **Consensus:** each cluster gets a consensus taxon (`--consensus_fraction`).
9. **Reports:**
   - database taxonomic distribution
   - optional target-taxon coverage (`--taxon`)
   - a MultiQC report per primer/database pair, plus one combined report

See [barbeque.md](barbeque.md) for the analysis steps in detail.

## Reference installation mode

`--build_references` skips the analysis. It installs the configured databases, the FooDMe2 primer
FASTAs and, by default, the NCBI taxdump and accession-to-taxid mapping. Everything goes under:

```text
<reference_base>/barbeque/<reference_version>/
```

See [build_references.md](build_references.md).

## Completion

`PIPELINE_COMPLETION` always runs last and handles final bookkeeping.
