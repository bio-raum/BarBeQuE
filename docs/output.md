# Outputs

Outputs are written to `--outdir` (default `results`).

## Directories

| Directory | Contents |
| --- | --- |
| `reports/individual/` | One MultiQC report per primer/database pair. |
| `reports/combined/` | One MultiQC report covering all primers and databases (`all_primers`). |
| `pipeline_info/` | Staged samplesheet and software versions. |
| `primers/` | Samplesheet generated from primer FASTA or `--primer_set` input, plus parser warnings. |
| `raw/obipcr/` | Raw obipcr amplicons. |
| `filtered_obipcr/` | Amplicons that match the original (unmerged) primers, plus `*_filter_stats.tsv`. Primer FASTA input only. |
| `parsed_obipcr/` | One amplicon table per primer/database pair (see below). |
| `build_db_taxids/` | Accession-to-taxid and taxid-count tables per database. |
| `db_distribution/` | Taxonomic composition of each database. |
| `cluster_fast/` | VSEARCH cluster files (`.uc`). |
| `join_accession_taxonomy/` | Cluster members joined to their taxonomy. |
| `consensus/` | Consensus taxonomy per cluster. This is the main result. |
| `taxid_filtered/` | Databases restricted to `--taxid`. Only with `--taxid`. |
| `tax_coverage/` | Coverage tables and Newick trees for the target taxon. Only with `--taxon`. |
| `taxonomic_coverage_blast/` | Per primer/database: `blast/*.blast.tsv` (centroid hits in core_nt for the configured taxid) and `with_region/*.with_region.txt` (distinct hit accessions). The MultiQC report counts these accessions. Only with `--taxonomic_coverage_blast`. |
## Consensus table

`consensus/` has one TSV per primer/database pair, with no header row. Each row is one accession:

1. cluster id
2. accession
3. accession taxid
4. accession taxon name
5. assigned name
6. assigned taxid
7. assigned rank
8. disambiguation (`;`-separated taxa found in the cluster)

A cluster whose accessions span several taxa can't tell those taxa apart at the chosen
`--cluster_id`. The assigned taxon must be supported by at least `--consensus_fraction` (default
`1.0`) of the accessions that have a valid taxid. Clusters without a consensus are `Unclassified`.

## Parsed obipcr table

`parsed_obipcr/` has one TSV per primer/database pair, with these columns:

- `Sequence_ID`, `Amplicon_Length`, `Amplicon_GC`, `Amplicon_Sequence`
- for each of `Forward_` and `Reverse_`: `Primer`, `Match`, `Errors`,
  `Mismatch_Positions_Primer`, `Mismatch_From_3Prime`, `Primer_GC`, `Match_GC`
