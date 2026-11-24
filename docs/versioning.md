# Versioning

BarBeQuE reproducibility depends on three things:

1. pipeline code version
2. reference data version
3. external database release dates

## Pipeline Version

Pin production runs by checking out a release tag in your clone:

```bash
git checkout <tag>
nextflow run main.nf ...
```

Use release tags for published analyses. Use branches such as `main` only for development.

## Reference Version

Reference installs are scoped by:

```bash
--reference_version
```

Default:

```text
1.1
```

Installed paths use:

```text
<reference_base>/barbeque/<reference_version>/
```

Install a new reference version beside the old one instead of overwriting an existing analysis reference set.

## Database Versions

Some upstream databases are versioned, and some are effectively snapshots of whatever was current on installation day.

Record at least:

- BarBeQuE git tag or commit
- `--reference_version`
- `--midori_version`
- date references were installed
- selected `--dbs` or custom database filename/checksum
- NCBI taxdump and accession mapping source/date

## Recommended Run Metadata

Keep the following with each analysis:

```text
nextflow run command
pipeline version/tag
reference_base
reference_version
midori_version
database ids
primer input file or primer_set names
outdir
run date
```
