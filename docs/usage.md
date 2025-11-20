# Usage information

[Basic execution](#basic-execution)

[Pipeline version](#specifying-pipeline-version)

[Basic options](#basic-options)

[Expert options](#expert-options)

## Basic execution

Please see our [installation guide](installation.md) to learn how to set up this pipeline first. 

A basic execution of the pipeline looks as follows:

a) Without a site-specific config file

```bash
nextflow run bio-raum/BarBeQuE -profile singularity --input samples.tsv \\
--reference_base /path/to/references \\
--run_name pipeline-test
```

where `path_to_references` corresponds to the location in which you have [installed](installation.md) the pipeline references (this can be omitted to trigger an on-the-fly temporary installation, but is not recommended in production). 

In this example, the pipeline will assume it runs on a single computer with the singularity container engine available. Available options to provision software are:

`-profile singularity`

`-profile apptainer`

`-profile docker` 

`-profile podman` 

`-profile conda` 

Additional software provisioning tools as described [here](https://www.nextflow.io/docs/latest/container.html) may also work, but have not been tested by us. Please note that conda may not work for all packages on all platforms. If this turns out to be the case for you, please consider switching to one of the supported container engines. 

**IMPORTANT** We do not recommend you use Conda for production purposes due to issues with reproducibility of environments across platforms and time. For a discussion, see [here](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30140-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471218301406%3Fshowall%3Dtrue)

b) with a site-specific config file

```bash
nextflow run bio-raum/BarBeQuE -profile my_profile --input samples.tsv \\
--run_name pipeline-test \\
--dbs midoria_lrna
```

In this example, both `--reference_base` and the choice of software provisioning are already set in the local configuration `lsh` and don't have to be provided as command line argument. 

## Specifying pipeline version

If you are running this pipeline in a production setting, you will want to lock the pipeline to a specific version. This is natively supported through nextflow with the `-r` argument:

```bash
nextflow run bio-raum/BarBeQuE -profile my_profile -r 1.0 <other options here>
```

The `-r` option specifies a github [release tag](https://github.com/marchoeppner/THIS_PIPELINE/releases) or branch, so could also point to `main` for the very latest code release. Please note that every major release of this pipeline (1.0, 2.0 etc) comes with a new reference data set, which has the be [installed](installation.md) separately.

## Basic options

### `--input` [default = null ]

This pipeline expects a samplesheet with information on one or several primer sets for benchmarking. The format is a simple tab-delimited text file (.tsv):

```TSV
primer  fwd rev min max
16SMeat  GACGAGAAGACCCTRTGGAGC   TCCRAGRTCGCCCCAAYC  50  100
```

| Column | Description |
| ------ | ----------- |
| primer | The name of the primer set, determines naming of output files |
| fwd    | The sequence of the forward primer; can include ambigious (IUPAC) bases |
| rev    | The sequence of the reverse primer; can include ambigious (IUPAC) bases |
| min    | minimum length of expected amplicons |
| max    | maximum length of expected amplicons |

Please make sure that your values for min/max are somewhat realistic for your primer set; else the results may be very noisy and unreliable. 

### `--dbs` [default = null]

A list of one or several pre-installed databases to benchmark against. If multiple databases are requested, they have to be separated by a comma

```bash
nextflow run bio-raum/BarBeQuE --input primers.tsv --dbs midori_lrrna,midori_cytb ...
```

### `--run_name` [default = null]

A descriptive name for this analysis run. Should be a single workd without special characters or white spaces  (i.e. my_analysis_run_01). 

### `--list_dbs` [default=null]

List all pre-installed databases and exit. Some examples include:

| Name | Description |
| ---- | ----------- |
| refseq_mito | The RefSeq dataset of mitochondrial genomes (v1.1) |
| midori_lrna | The lrRNA database from MIDORI (aka 16S) |
| midori_cytb | The CYTB database from MIDORI |
| mitofish | The Mitofish database |

## Expert options

### `--taxon` [default=null]

The default mode of this analysis is to run against the entire target database; use this option to focus on a specific [taxonomic sub group](https://www.ncbi.nlm.nih.gov/taxonomy) and get additional information/visualization. The argument must be a valid taxon identifier (such as: 'Chordata' or 'Mammalia'). For the moment, the pipeline cannot validate your taxon argument and we found that some combinations of databases and perfectly valid taxon arguments will nevertheless crash CRABS.

## Expert options

### `--crabs_insilicopcr_options` [ default = null ]

Pass custom options to the insilico PCR stage of the CRABS analysis. This may be used to e.g. increase to decrease the number of allowed mismatched bases (`--mismatch 2`). 