# sc_preprocess

A Nextflow pipeline for preprocessing of 10X scRNA-seq data using Alevin.

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/basic_test.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>


## Workflow

- download genome, transcriptome and GTF reference files from [GENCODE version M25](https://www.gencodegenes.org/mouse/release_M25.html).

- create an index of spliced- and unspliced transcript sequences based on the reference annotations, following the [Alevin Velocity tutorial(https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

- create a merged of the exon+intron transcriptome with the mouse genome into a `gentrome` and index it with `salmon`

- perform cell barcode detection, UMI deduplication and read quantification with [Alevin](https://salmon.readthedocs.io/en/latest/)

- read quantifications into R with `tximeta`, split into spliced- and unspliced counts and save as MatrixMarket (mtx) to disk for easy distribution of data

## Usage

As the workflow is hosted on GitHub it can be pulled directly via nextflow. On the HPC we use:

```bash

NXF_VER=21.04.1 nextflow run atpoint/sc_preprocess -r <commit_sha> -profile singularity,slurm --fastq 'path/to/*_{1,2}.fastq.gz'  

```

**Details:**

-  `--fastq`: path to fastq file pairs, e.g. `--fastq path/to/*_{1,2}.fastq.gz`

-  `-profile`: choose a container profile, either `docker`, `singularity` or `conda`.
This repository contains an `environment.yml` that was used to build a [Docker image](https://hub.docker.com/r/atpoint/sc_preprocess). Choosing `-profile singularity` will pull it and convert to a Singulartiy image.
We also added a platform-agnostic environment (`environment_conda.yml`) which sould solve on both Linux and macOS.
The `-profile conda` uses the latter environment. For testing one can add `-profile test`. With `-profile slurm` nextflow will submit jobs via the SLURM scheduler.
The `params.queue/clusteroptions` can be used to specify a queue and additional arguments for the scheduler.

-  each process has options to change the threads and memory, intuitively named in the config file. The defaults are tailored for a HPC/cluster node.
The `-profile test` sets threads/memory to a minimum so it runs on any machine and is suitable for the GitHub Actions CI tests.

We use a naming convention of the fastq files like:

-  `basename_1.fastq.gz` for R1 (CB/UMI)

-  `basename_2.fastq.gz` for R2 (cDNA)

-  by default the pipeline will create folders for the output file in `$launchDir` so the directory from which the Nextflow run was launched. This can be changed with the `--*_outdir` params,
see the `nextflow.config` file. params are named intuitively.

## Software

A Docker container is available at the [Docker Hub](https://hub.docker.com/r/atpoint/sc_preprocess) which can be used using `-profile docker/singularity`. If using `-profile conda` then the `environment_conda.yml` will be used to create a conda environment with all required software.

## Citations

-  [nf-core project](https://nf-co.re/)

-  [Ewels et al (2020) Nature Biotechnology volume 38, pages 276–278](https://www.nature.com/articles/s41587-020-0439-x)

-  [Nextflow Docs](https://www.nextflow.io/docs/latest/index.html#)

-  [Seqera Training](https://seqera.io/training/)

-  [https://github.com/nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf)

-  [Merkel, D (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal](https://dl.acm.org/doi/10.5555/2600239.2600241)

-  [Kurtzer et al (2017) Singularity: Scientific containers for mobility of compute. PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)

-  [Srivastava et al (2019) Alevin efficiently estimates accurate gene abundances from dscRNA-seq data. Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y)

-  [Soneson et al (2021) Preprocessing choices affect RNA velocity results for droplet scRNA-seq data. PLoS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008585)

-  [Love et al (2020) Tximeta: Reference sequence checksums for provenance identification in RNA-seq](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007664)
