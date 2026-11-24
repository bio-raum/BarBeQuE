# Troubleshooting

## No Primer Input Was Accepted

Normal benchmarking requires exactly one primer source:

- `--input <samplesheet.tsv>`
- `--input <primer_fasta_directory>`
- `--primer_set <name[,name...]>`

Do not combine `--input` and `--primer_set`.

## FASTA Input Fails

When `--input` is a primer FASTA, or a directory of them, provide global amplicon bounds:

```bash
--primer_min 100 --primer_max 500
```

A directory must contain `.fa`, `.fasta`, or `.fna` files; everything else in it is ignored.

If a FASTA is rejected, the error names every unusable file at once and nothing is run. The
causes are: no FASTA records in the file, a record with an empty sequence, a non-nucleotide
character in a sequence, a prefix missing its forward or reverse primer, or a record with no
`fwd`/`rev` token in a file that is not a plain two-record pair.

A warning that a prefix "could not be collapsed" is not a failure - it means that prefix's
variants have different lengths, so each forward/reverse combination is being benchmarked as its
own primer set. Check `primers/` to see exactly what was run.

## Accession Taxonomy Mapping Is Missing

BarBeQuE needs accession-to-taxid data for taxonomy-aware steps. Provide either:

```bash
--reference_base /path/to/references
```

or:

```bash
--taxdump /path/to/new_taxdump \
--accession_taxonomy /path/to/nucl_gb.accession2taxid
```

## Database Id Is Rejected

Run:

```bash
nextflow run main.nf -profile conda --list_dbs --reference_base /path/to/references
```

Use ids exactly as shown. Prebuilt BLAST databases such as `core_nt` are not valid `--dbs` inputs for in-silico PCR.

## No Amplicons Are Produced

Common causes:

- primers target a region absent from the selected database
- `min` and `max` amplicon bounds are too strict
- mismatch settings are too strict
- gene-specific databases do not include primer binding sites outside the gene

Try a broader database such as `refseq_mito` and inspect raw in-silico PCR output.

## `--taxon` Runs Slowly

Target-taxon coverage can be expensive for large groups because the workflow traverses taxonomy and checks many species. Use a narrower taxon when possible.

## Full Test Runner Fails On Conda

`run_all_tests.sh` executes Python tests first, then Nextflow module tests. The module tests require the selected process backend. If the run fails with `conda: command not found`, install Conda/Mamba or run with a container backend/profile where the module environments can be resolved.

## Remote Config Include Fails

`nextflow.config` includes a shared remote config by default. If your environment blocks network access, provide a local config with `-c` or adjust `params.custom_config_base` as appropriate for your deployment.
