# Usage information

[Basic execution](#basic-execution)

[Pipeline version](#specifying-pipeline-version)

[Resources](#resources)

## Basic execution

Please see our [installation guide](installation.md) to learn how to set up this pipeline first. 

A basic execution of the pipeline looks as follows:

a) Without a site-specific config file

```bash
nextflow run marchoeppner/THIS_PIPELINE -profile singularity --input samples.tsv \\
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
nextflow run marchoeppner/THIS_PIPELINE -profile my_profile --input samples.tsv \\
--run_name pipeline-test 
```

In this example, both `--reference_base` and the choice of software provisioning are already set in the local configuration `lsh` and don't have to be provided as command line argument. 

## Specifying pipeline version

If you are running this pipeline in a production setting, you will want to lock the pipeline to a specific version. This is natively supported through nextflow with the `-r` argument:

```bash
nextflow run marchoeppner/THIS_PIPELINE -profile my_profile -r 1.0 <other options here>
```

The `-r` option specifies a github [release tag](https://github.com/marchoeppner/THIS_PIPELINE/releases) or branch, so could also point to `main` for the very latest code release. Please note that every major release of this pipeline (1.0, 2.0 etc) comes with a new reference data set, which has the be [installed](installation.md) separately.
