# Main Analysis Workflow

The main analysis lives in `workflows/barbeque.nf`. It runs after `DATABASE` has produced the selected database FASTA channel.

## Running it

Always pick one software profile with `-profile` (one dash): `conda`, `docker` or `apptainer`.

```bash
nextflow run main.nf -profile conda --input primers.tsv --dbs midori_cytb --run_name my_test --reference_base /path/to/references
```

BarBeQuE works best with `-profile conda`. Some modules don't have a ready-made Docker container. With `docker` or `apptainer`, Nextflow (through Wave) builds a container for those modules from their conda packages, so the first run needs internet and takes longer.

`--run_name`: a name for your run, e.g. `--run_name my_test`.

## Inputs

The primers can be submitted in three ways. Use only one of them. With `--input`, give the location of a samplesheet with the columns `primer`, `fwd`, `rev`, `min` and `max`. If you have several variants of a primer, please submit them as a `.fasta` file instead, with `fwd` or `rev` in each header (e.g. `>ITS2_fwd1`, `>ITS2_rev1`). If a FASTA is passed, `--primer_min` and `--primer_max` must also be specified.

If you prefer to test any of the primers which come by default with FooDMe2, please use the `--primer_set` option followed by the name of the primer set. These primer sets are resolved from the FooDMe2 catalog.

- `--primer_set`: named primer sets resolved from the FooDMe2 catalog.
- `--input <samplesheet.tsv>`: a TSV with `primer`, `fwd`, `rev`, `min`, and `max` columns.
- `--input <primers.fasta>` or `--input <directory>`: one `.fa`/`.fasta`/`.fna` primer FASTA, or a directory of them. Requires `--primer_min` and `--primer_max`.

A file input is read as a FASTA or as a samplesheet based on its first non-blank line, not its
extension (`WorkflowPipeline.isFastaInput`). FASTA input - and every `--primer_set` download - is
converted to a samplesheet by `modules/parse_primers` first, so all three routes reach `INPUT_CHECK`
in the same shape and share one validation and staging step.

See [Primer Input](primer_input.md) for the parsing rules and set naming.

Database input is one of:

- `--dbs`: comma-separated installed databases, i.e. the preinstalled ones you got when you built the reference base. e.g. `--dbs 'midori_cytb'`
- `--custom_db`: one or more user-provided nucleotide FASTAs (`.fasta`, `.fa`, `.fna`). The accessions in the headers are looked up in the NCBI accession-to-taxid mapping to find each record's species.

Pass just one kind of database: either a custom FASTA via `--custom_db` or a preinstalled one via `--dbs`. If both are given, only `--custom_db` is used.

If you for some reason need to restrict the passed database to a smaller subset, you can pass `--taxid`, e.g. `--taxid 4442` (Camellia sinensis) keeps only the records of that taxon and the taxa below it. If you further need to filter the FASTA, switch on `--db_filter`. It removes records by the information in the header using `--db_filter_pattern 'environmental|recombinant'` (no spaces around the `|`). The reference sequences can also be length filtered using `--db_filter_min_length 100 --db_filter_max_length 200`.

Taxonomy input is resolved from:

- `--taxdump`, otherwise `<reference_base>/barbeque/<reference_version>/taxonomy/new_taxdump`
- `--accession_taxonomy`, otherwise `<reference_base>/barbeque/<reference_version>/taxonomy/nucl_gb.accession2taxid`

If you have given the reference base location, then the taxdump and accession taxonomy are auto-detected from there. Else you can give them directly with `--taxdump` and `--accession_taxonomy`.

`--mask`: cuts the amplicons down to what a sequencer would actually read. By default (paired-end) it keeps `--read_length` bases (default 150) at each end, joined by 20 `N`s. With `--single_end` it keeps only the first `--read_length` bases (default 400).

## Workflow Steps

1. Resolve primers into a common structure: primer id, forward primer, reverse primer, minimum amplicon length, maximum amplicon length.
2. Combine every primer with every selected database.
3. Run in-silico PCR with OBI-PCR.
4. For primer FASTA input, re-check every amplicon against the original, uncollapsed primers and discard those only the combined degenerate primer could explain (`--filter_collapsed_primers`).
5. Parse the raw OBI-PCR output.
6. Write amplicon length summaries.
7. Drop primer/database pairs with no retained amplicons from downstream analysis.
8. Optionally apply `--mask` to mimic single-end or paired-end read coverage.
9. Build accession-to-taxid tables once per database.
10. Cluster amplicons with `vsearch --cluster_fast` using `--cluster_id`.
11. Parse `.uc` cluster assignments.
12. Join clustered accessions to taxids.
13. Calculate consensus taxonomy per cluster.
14. Summarize database taxonomic distribution.
15. Optionally run target-taxon coverage with `--taxon`.
16. Build one MultiQC report per primer/database combination.

## Interactive Results

`--interactive` launches the Streamlit workflow after consensus files are produced. It reads the finished output directory and serves the dashboard on the configured Streamlit port.
