# sc_preprocess

A Nextflow pipeline for preprocessing of 10X scRNA-seq data using [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) based on the [COMBINE lab tutorial](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/).

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>


## Workflow

The workflow first creates a mapping index using both the spliced- and unspliced transcripts as well as the entire genome (as a mapping decoy) using Salmon. The spliced- and unspliced transcripts are parsed from a reference GTF using [eisaR](https://bioconductor.org/packages/release/bioc/html/eisaR.html). It then uses Alevin for barcode detection, UMI deduplication and Salmon for read quantification. The obtained quantifications are then split into spliced- and unspliced counts at the gene level using [tximeta](https://bioconductor.org/packages/release/bioc/html/tximeta.html) and then saved in as `mtx.gz` files together with the column- and row annotations.

## Usage / Parameters

On HPC using SLURM and Singularity use:

```bash

NXF_VER=21.04.1 nextflow run atpoint/sc_preprocess -r <commit_sha> -profile singularity,slurm --fastq 'path/to/*_{1,2}.fastq.gz'  

```

### Params

- pubmode: publishing mode, default relllink
- fastq: path to fastq file pairs, naming conventions must be `basename_1/2.fastq.gz`
- ref_genome: path or remote link to `fa.gz` genome fasta
- ref_gtf: path or remote link to `gtf.gz` GTF file
- ref_gene_name: name of gene name column in GTF, default `gene_name`
- ref_gene_id: name of gene id column in GTF, default `gene_id`
- ref_gene_type: name of gene type column in GTF, default `gene_type`
- ref_chrM: name of mitochondrial chromosome, default `chrM`
- cdna_readlength: length of cDNA (R2) read, default 91bp
- parse_intron_threads: threads for exon/intron parsing, leave at 1
- parse_intron_mem: memory for exin/intron parsing, default 16.GB, enough for mouse, not tested for human
- idx_outdir: name for index file directiry, default `./alevinIndex`
- idx_name: name of the folder containing the salmon/alevin idx files, default `idx`
- idx_salmonargs: further arguments beyond [ -t -d -i -p ] for salmon/alevin indexing, default `--sparse --gencode'`, we generally expect GENCODE reference files!
- idx_threads: threads for indexing, default 32
- idx_mem: memory for indexing, default 50.GB
- skip_quant: logical, whether to stop after indexing and skip quantification
- quant_outdir: output dir for quantifications, default `./alevinQuant`
- quant_libtype: library type, leave at default `ISR` for 10X data
- quant_additional: additional quant args, including the chemistry, beyond [ -i --tgMap --mrna --rrna -l -p -o -1 -2 ], default `--chromiumV3` is appropriate for 10X v3 experiments
- quant_threads: threads for quantification, default 8
- quant_mem: memory for quantification, default 40.GB
- skip_mtx: logical, whether to save counts as `mtx.gz` files
- mtx_outdir: outdir for `mtx.gz` files, default `./mtx`
- mtx_threads: threads for mtx generation, leave at 1
- mtx_mem: memory for mtx generation, default 8.GB
- queue: name of cluster queue, default `normal` to submit to normal queue
- clusteroptions: options for cluster scheduler, default `--time=12:00:00` to submit a 12h job

## Software

A Docker container is available at the [Docker Hub](https://hub.docker.com/r/atpoint/sc_preprocess) which can be used using `-profile docker/singularity`. If using `-profile conda` then the `environment.yml` will be used to create a conda environment with all required software.

## Citations

-  [nf-core website](https://nf-co.re/)

-  [nf-core paper => Ewels et al (2020) Nature Biotechnology volume 38, pages 276–278](https://www.nature.com/articles/s41587-020-0439-x)

-  [Nextflow Docs](https://www.nextflow.io/docs/latest/index.html#)

-  [Seqera Training](https://seqera.io/training/)

-  [https://github.com/nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf)

-  [Merkel, D (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal](https://dl.acm.org/doi/10.5555/2600239.2600241)

-  [Kurtzer et al (2017) Singularity: Scientific containers for mobility of compute. PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)

-  [Srivastava et al (2019) Alevin efficiently estimates accurate gene abundances from dscRNA-seq data. Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y)

-  [Soneson et al (2021) Preprocessing choices affect RNA velocity results for droplet scRNA-seq data. PLoS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008585)

-  [Love et al (2020) Tximeta: Reference sequence checksums for provenance identification in RNA-seq](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007664)
