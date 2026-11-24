# Software

BarBeQuE is a Nextflow DSL2 workflow. Process-level software is declared in the individual module `environment.yml` files and module definitions.

## Main Runtime

| Tool | Use | Citation / DOI |
| --- | --- | --- |
| Nextflow | Workflow execution | Di Tommaso et al. 2017, Nature Biotechnology, DOI: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| Java | Nextflow runtime | No single project DOI used here. |
| Conda / Docker / Singularity / Apptainer / Podman | Software provisioning and containers | Cite the selected platform separately if required by your publication or infrastructure policy. |

## Bioinformatics Tools

| Tool | Use | Citation / DOI |
| --- | --- | --- |
| OBITools / OBITools4 | Default in-silico PCR engine (`obipcr`) | Boyer et al. 2016, Molecular Ecology Resources, DOI: [10.1111/1755-0998.12428](https://doi.org/10.1111/1755-0998.12428) |
| VSEARCH | Dereplication and clustering | Rognes et al. 2016, PeerJ, DOI: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584) |
| SeqKit | Sequence summaries and FASTA filtering | Shen et al. 2016, PLOS ONE, DOI: [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) |
| TaxonKit | NCBI taxonomy lineage utilities | Shen and Xiong 2021, Journal of Genetics and Genomics, DOI: [10.1016/j.jgg.2021.03.006](https://doi.org/10.1016/j.jgg.2021.03.006) |
| ETE Toolkit / ETE3 | Taxonomy and tree handling | Huerta-Cepas et al. 2016, Molecular Biology and Evolution, DOI: [10.1093/molbev/msw046](https://doi.org/10.1093/molbev/msw046) |
| GrapeTree | Tree visualization / clustering helper dependency | Zhou et al. 2018, Genome Research, DOI: [10.1101/gr.232397.117](https://doi.org/10.1101/gr.232397.117) |
| MultiQC | Report aggregation | Ewels et al. 2016, Bioinformatics, DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |

## Python And Analysis Libraries

| Tool | Use | Citation / DOI |
| --- | --- | --- |
| Python | Helper scripts and dashboard runtime | No single project DOI used here. |
| Biopython | FASTA parsing/writing in helper scripts | Cock et al. 2009, Bioinformatics, DOI: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) |
| pandas | Tabular data handling | The pandas development team 2020, DOI: [10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134); McKinney 2010, DOI: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a) |
| NumPy | Numerical arrays | Harris et al. 2020, Nature, DOI: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| Plotly | Dashboard/report plotting | No single formal DOI used here. |
| Streamlit | Interactive dashboard | No single formal DOI used here. |
| taxidTools | Python taxonomy helper | No formal DOI identified; cite the package repository if needed. |

## Command-Line Utilities Without Project DOIs

Some modules also use standard command-line utilities such as `awk`/`gawk`, `grep`, `sed`, `tar`, `gzip`/`gunzip`, `wget`, and `coreutils`. These are runtime dependencies rather than primary scientific methods in this pipeline; no DOI is listed here.

## Reproducibility

Use a container profile for production runs when possible. Conda environments are convenient but may resolve differently over time or across platforms.
