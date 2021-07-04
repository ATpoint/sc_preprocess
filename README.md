# sc_preprocess

A Nextflow pipeline for preprocessing of 10X scRNA-seq data using Alevin.


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

-  `-profile`: choose a container profile, `-profile singularity` or `-profile docker` (see #software)
This repository contains an `environment.yml` to create a conda environment. This is not available directly via the workflow so one would need to manually create an environment, e.g. `mamba create --name sc_preprocess --file environment.yml` and then run the workflow locally.  
For testing one can add `-profile test`. With `-profile slurm` nextflow will submit jobs via the SLURM scheduler.

-  `--idx_label` & `--quant_label` allocate resources:  
=> for the indexing there is `--idx_label idx_small/idx_big` which allocates 4 or 36 cores and 15 or 40GB of RAM.  
=> for the quantification there is `--quant_label quant_small/quant_big` which allocates 3 or 36 cores with 6 or 50GB RAM.  
=> the defaults are the big labels intended for HPC use.  

We use a naming convention of the fastq files like:

-  `basename_1.fastq.gz` for R1 (CB/UMI)

-  `basename_2.fastq.gz` for R2 (cDNA)

## Software

A Docker container is available at [Docker Hub](https://hub.docker.com/r/atpoint/sc_preprocess).  
When using the Docker or Singularity profile Nextflow will take care of pulling it.

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

- [Love et al (2020) Tximeta: Reference sequence checksums for provenance identification in RNA-seq](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007664)