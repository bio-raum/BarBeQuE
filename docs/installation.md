# Installation

[Installating Nextflow](#installing-nextflow)

[Installing references](#installing-the-references)

[Site-specific config](#site-specific-config-file)

[Testing](#testing-the-installation)

## Installing nextflow

Nextflow is a highly portable pipeline engine. Please see the official [installation guide](https://www.nextflow.io/docs/latest/getstarted.html#installation) to learn how to set it up.

This pipeline expects Nextflow version 24.10.8 or later, available [here](https://github.com/nextflow-io/nextflow/releases/tag/v24.10.8).

## Software provisioning

This pipeline is set up to work with a range of software provisioning technologies - no need to manually install packages. 

You can choose one of the following options:

[Docker](https://docs.docker.com/engine/install/)

[Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/)

[Apptainer](https://apptainer.org/docs/admin/main/installation.html)

[Podman](https://podman.io/docs/installation)

[Conda](https://github.com/conda-forge/miniforge)

The pipeline comes with simple pre-set profiles for all of these as described [here](usage.md); if you plan to use this pipeline regularly, consider adding your own custom profile to our [central repository](https://github.com/bio-raum/nf-configs) to better leverage your available resources. 

**IMPORTANT** We do not recommend you use Conda for production purposes due to issues with reproducibility of environments across platforms and time. For a discussion, see [here](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30140-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471218301406%3Fshowall%3Dtrue)

## Installing the references

This pipeline requires locally stored databases in [CRABS-compliant](https://github.com/gjeunen/reference_database_creator) format. To build these, do:

```
nextflow run bio-raum/BarBeQuE -profile singularity \\
-r <DESIRED_VERSION> \\
--build_references \\
--run_name build_refs \\
--reference_base /path/to/references
```

where `/path/to/references` could be something like `/data/pipelines/references` or whatever is most appropriate on your system. The <DESIRED_VERSION> should be set to whatever version of BarBeQue you plan on using. Databases will typically only be updated with every major release, so installation is only required once per major release (i.e. 1.0). Also see [versioning.md](versioning.md)

If you do not have singularity on your system, you can also specify docker, apptainer, podman or conda for software provisioning - see the [usage information](usage.md).

Please note that the build process will create a pipeline-specific subfolder that must not be given as part of the `--reference_base` argument. This pipeline is part of a collection of pipelines that use a shared reference directory and it will choose the appropriate subfolder by itself. 

## Site-specific config file

If you run on anything other than a local system, this pipeline requires a site-specific configuration file to be able to talk to your cluster or compute infrastructure. Nextflow supports a wide range of such infrastructures, including Slurm, LSF and SGE - but also Kubernetes and AWS. For more information, see [here](https://www.nextflow.io/docs/latest/executor.html).

Site-specific config-files for our pipeline ecosystem are stored centrally on [github](https://github.com/bio-raum/nf-configs). Please talk to us if you want to add your system.

### Custom config

If you absolutely do not want to add your system to our central config repository, you can manually pass a compatible configuration to nextflow using the `-c`  command line option:

```bash
nextflow -c my.config run bio-raum/BarBeQuE -profile myprofile -r 1.0.1 --input samples.csv --run_name my_run_name --reference_base /path/to/references
```

A basic example using Singularity may look as follows:

```GROOVY
params {
  reference_base = "/path/to/references"
}

process {
  resourceLimits = [ cpus: 16, memory: 64.GB, time: 72.h ]
}

singularity {
  enabled = true
  cacheDir = "/path/to/singularity_cache"
}
``` 
This would be for a single computer, with 16 cores and 64GB Ram, using Singularity. Containers are cached to the specified location to be re-used on subsequent pipeline runs.  

Or with the Conda/Mamba package manager:

```GROOVY
params {
  reference_base = "/path/to/references"
}

process {
  resourceLimits = [ cpus: 16, memory: 64.GB, time: 72.h ]
}

conda {
  enabled = true
  useMamba = true
  cacheDir = "/path/to/conda_cache"
}
```

Finally, if you are planning to run this pipeline on a cluster, your config my look as follows:

```GROOVY
params {
  reference_base = "/path/to/references"
}

process {
  executor = 'slurm'
  queue = 'all' 
  resourceLimits = [ cpus: 20, memory: 250.GB, time: 240.h ]
}

singularity {
  enabled = true
  cacheDir = "/path/to/singularity_cache"
}
``` 

Here, we specify that the pipeline should use the SLURM resource manager for job submission to a partition called "all", and that the configuration of individual nodes is 20 cores with 250GB of RAM. Again, we use Singularity (which could be switched to Conda, Apptainer etc - whatever fits your situation).

## Testing the installation

BarBeQue comes with a simple test data sets, which you can run to verify that everything works. To run the test data, combine your profile of choice with the built-ind test profile:

```bash
nextflow run bio-raum/BarBeQue -r DESIRED_VERSION -profile singularity,test
```

where `DESIRED_VERSION` should be one of the releases of the pipeline (or 'main' for the latest version of the code).