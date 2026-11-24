# Usage

## Quick start

```bash
nextflow run main.nf \
  -profile conda \
  --input primers.tsv \
  --dbs refseq_mito \
  --reference_base /path/to/references \
  --run_name primer_benchmark \
  --outdir results
```

A run needs:

- **Primers:** exactly one of `--input` or `--primer_set`.
- **A database:** `--dbs` (installed) or `--custom_db` (your own FASTA). `--run_name` is required with `--dbs`.
- **Taxonomy:** `--reference_base`, or `--taxdump` plus `--accession_taxonomy`.

## Reference data

Install the databases and taxonomy once:

```bash
nextflow run main.nf -profile conda \
  --build_references --reference_base /path/to/references
```

Later runs then only need `--reference_base`. See [Reference Installation](build_references.md).

To use your own database without a reference base:

```bash
--custom_db /path/to/custom.fasta \
--taxdump /path/to/new_taxdump \
--accession_taxonomy /path/to/nucl_gb.accession2taxid
```

If both `--custom_db` and `--dbs` are given, only `--custom_db` is used.

`--list_dbs` prints the installed database ids and exits.

## Primer input

### Samplesheet

A tab-separated file with the columns `primer`, `fwd`, `rev`, `min` and `max`:

```text
primer  fwd                         rev                         min  max
COI     GGWACWGGWTGAACWGTWTAYCCYCC  TAIACYTCIGGRTGICCRAARAAYCA  100  500
```

### Primer FASTA

`--input` also accepts a primer FASTA file or a directory of them. Each file is one
primer set, named after the file. In a directory, only files ending in `.fa`, `.fasta`
or `.fna` are read. Subdirectories are ignored.

FASTA input requires amplicon bounds:

```bash
--input primer_fastas/ --primer_min 100 --primer_max 500
```

See [Primer Input](primer_input.md) for the header and file rules.

### Named primer sets

```bash
--primer_set COI,Fish16S
```

Named sets come from the [FooDMe2](https://github.com/bio-raum/FooDMe2) catalog. Each
set brings its own amplicon bounds. `--list_primers` prints the catalog and exits.
Both options need internet access.

## Options

| Option | Default | Description |
|---|---|---|
| `--taxid` | – | Restrict the database to this taxon and its descendants. |
| `--db_filter` | `false` | Clean the database. Tune it with `--db_filter_pattern`, `--db_filter_min_length`, `--db_filter_max_length` and `--db_filter_max_n`. |
| `--obipcr_mismatches` | `2` | Maximum mismatches per primer. |
| `--primer_site_not_found_mismatches` | `4` | In primer-site diagnostics, classify a candidate as `PRIMER_SITE_NOT_FOUND` only above this mismatch count. Must be at least `--obipcr_mismatches`; intermediate counts remain `MISMATCH` and are written to `coverage_analysis/primer_mismatch_sites/<sample>.{forward,reverse}_mismatch_sites.fasta`. |
| `--primer_site_window` | `100` | In primer-site diagnostics, search this many bases on each side of the primer site predicted from the BLAST hit. |
| `--obipcr_fixed_3prime` | `3` | Number of 3' bases where no mismatch is allowed (`0` = off). |
| `--obipcr_circular` | `true` | Treat database records as circular, so amplicons crossing the start/end of a complete organelle genome are found. Applies to all records, including linear ones. |
| `--filter_collapsed_primers` | `true` | Drop amplicons that match only the merged degenerate primer. Applies to FASTA input only. |
| `--cluster_id` | `0.98` | `vsearch --cluster_fast` identity threshold (0–1). |
| `--consensus_fraction` | `1.0` | Fraction of a cluster that must support a taxon (>0.5–1.0). `1.0` is a strict LCA. |
| `--mask` | `false` | Mask amplicons to simulate reads. Paired-end mode keeps `--read_length` bases at each end, joined by 20 `N`s. |
| `--single_end` | `false` | With `--mask`, truncate each amplicon to `--read_length` instead. |
| `--read_length` | `150` (`400` single-end) | Read length used by `--mask`. |
| `--taxon` | – | Add coverage reports for this target taxon. |
| `--taxonomic_coverage_blast` | `true` | BLAST the centroids against core_nt, restricted to `--taxonomic_coverage_taxid`. |
| `--taxonomic_coverage_taxid` | `3398` | Taxid used to restrict the core_nt search, including descendants. |
| `--core_nt_db` | Installed database under `--reference_base` | Override with the prefix of an existing core_nt BLAST database, such as `/refs/core_nt/core_nt`. |
| `--taxonomic_coverage_identity` | `75` | Minimum % identity of a centroid hit for a record to count as having the region. |
| `--taxonomic_coverage_query_cov` | `75` | Minimum % of the centroid that one hit must cover. |
| `--interactive` | `false` | Start the experimental Streamlit dashboard on port 8501 after the run. |

See [Outputs](output.md) for results and [Troubleshooting](troubleshooting.md) for common failures.
