# Installation

## Requirements

- Nextflow DSL2, version 24.10.5 or newer
- Java compatible with your Nextflow version
- One software provisioning backend:
  - Singularity
  - Apptainer
  - Docker
  - Podman
  - Conda

Container backends are recommended for production. Conda is useful for development and small local tests, but it is less reproducible across platforms and time.

## Install Nextflow

Follow the official Nextflow installation guide:

```text
https://www.nextflow.io/docs/latest/getstarted.html#installation
```

Check the installed version:

```bash
nextflow -version
```

## Get The Pipeline

BarBeQuE is not yet merged into bio-raum/BarBeQuE, so clone it from the fork and run `main.nf`
from the checkout:

```bash
git clone https://github.com/jissksaji/BarBeQuE.git
cd BarBeQuE
```

## Choose A Profile

| Profile | How software is provided |
| --- | --- |
| `conda` | A conda environment per module. Works offline once the environments are built. |
| `docker`, `podman`, `singularity`, `apptainer` | BioContainers images for most modules. Modules with no image yet (`obipcr`, `filter_obipcr`, `taxid_db_filter`, `taxonomic_coverage`, `cluster_consensus`, `db_distribution`, `download`, `streamlit`) get their image built by [Wave](https://seqera.io/wave/), which needs internet access. |

For cluster or cloud execution, provide a site-specific Nextflow config with `-c` or use a shared profile from the bio-raum config repository.

## Install References

```bash
nextflow run main.nf \
  -profile conda \
  --build_references \
  --reference_base /path/to/references
```

Do not include `barbeque/<reference_version>` in `--reference_base`; the pipeline adds that structure itself.

## Verify The Installation

List installed/known databases:

```bash
nextflow run main.nf \
  -profile conda \
  --list_dbs \
  --reference_base /path/to/references
```

List named primer sets from the upstream catalog:

```bash
nextflow run main.nf \
  -profile conda \
  --list_primers
```

Run a small local test when the required software backend is available:

```bash
bash run_all_tests.sh
```

The module tests require the selected backend and external bioinformatics packages. The Python unit tests at the start of `run_all_tests.sh` only require `python3`.
