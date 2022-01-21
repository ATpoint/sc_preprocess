# sc_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

## Introduction

**sc_preprocess** is an automated preprocessing pipeline for 10X scRNA-seq data implemented using [Nextflow](https://www.nextflow.io/) which is fully containerized to take care of all required software and ensure reproducibility. It supports feature barcode experiments such as CITE-Seq and cell hashing (HTO).

## Details

When running with default parameters the following steps are being executed:

1. Generate an expanded reference transcriptome containing all spliced (=exonic) and unspliced (intronic) transcript sequences. This expanded transcriptome is required for quantification to allow output of both **spliced- and unspliced counts**, e.g. for [velocity](https://www.embopress.org/doi/full/10.15252/msb.202110282) analysis. Also, this expanded reference is decoyed by the entire reference genome to perform [selective alignments](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) in order to improve mapping accuracy by capturing reads that better align to the genome than the transcriptome, removing potential gDNA contaminations and other spurious mapppings.

2. Quantification of the RNA reads against the expanded reference using [Alevin](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y) performing cell barcode detection, read mapping, UMI deduplication and gene count estimation, resulting in per-cell gene-level abundance estimations.

3. Split the obtained quantifications into spliced and unspliced count tables and save these in [mtx](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) format for easy distribution and loading into downstream analysis environments such as R.

4. Create a per-sample summary report using [alevinQC](https://csoneson.github.io/alevinQC/articles/alevinqc.html) which includes relevant QC metrics such as the average number of reads per cell, number of detected genes per cell [and others](https://bioconductor.org/packages/3.15/bioc/vignettes/alevinQC/inst/doc/alevinqc.html#generate-individual-plots).

5. Optionally, quantify matched reads from a [feature barcoding](https://www.biolegend.com/en-us/blog/cite-seq-and-totalseq-reagents) experiment such as [CITE-seq](https://cite-seq.com/) or [Cell Hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) against a provided set of barcode sequences. The obtained counts are filtered for the cellular barcodes detected in the RNA quantification. By default the pipeline performs [barcode translation](https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-) as the cellular barcodes for RNA and feature barcodes on the same gel bead are slightly different. This option is activated by default and makes sense when using TotalSeqB or TotalSeqC feature barcodes. This can be turned off with `--translate_barcodes false`. The translation list is taken from CellRanger the [CellRanger GitHub](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz), see also this thread at [biostars.org](https://www.biostars.org/p/9506747/). In any case, the feature barcode quantifications are written as `mtx.gz` in the same order as the cells detected in the RNA experiment.

## Usage

The typical CL on a HPC with Singularity and SLURM (see sections below) would be:

```bash

NXF_VER=21.10.6 \
    nextflow run main.nf -profile singularity,slurm \
    --genome path/to/genome.fa.gz --gtf path/to/gtf.gz \
    --samplesheet path/to/samplesheet.csv --features_file /path/to/features_file.tsv \
    -with-trace -with-report report.log -bg > pipeline.log

```    

### Samplesheet

The pipeline reads the fastq file pairs from a [samplesheet](https://github.com/ATpoint/sc_preprocess/blob/main/test/samplesheet.csv) which is a four-column CSV with a header.

```bash
sample_id,R1,R2,is_sf
sample1,/full/path/to/test/sample1_1.fastq.gz,/full/path/to/test/sample1_2.fastq.gz,false
sample1,/full/path/to/test/sample1a_1.fastq.gz,/full/path/to/test/sample1a_2.fastq.gz,false
sample1,/full/path/to/test/sample1SF_1.fastq.gz,/full/path/to/test/sample1SF_2.fastq.gz,true
sample2,/full/path/to/test/sample2_1.fastq.gz,/full/path/to/test/sample2_2.fastq.gz,false
```

The first line is the mandatory header followed by the sample to quantify:

- column1 (`sample_id`) is the per-sample name to be used in the output. This can be any Unix-compatible name, and there is no need to match this name with the fastq file names as other software may require you to do. 
- column2/3 (`R1/2`) are the paths to the fastq files for that sample. It either must be the full **absolute path** (don't use `~`) or alternative a path **relative** to the directory that this pipeline is started from using either of the three [implicit Nextflow variables](https://www.nextflow.io/docs/latest/script.html?highlight=basedir#implicit-variables) `$baseDir`, `$projectDir` and `$launchDir` (see example below).
- column4 (`is_sf`) is a logical column with either `true` or `false` indicating whether this fastq file pair is a feature barcode experiment (sf=surface features). If so then these files will be quantified against the barcode library provided by `--features-file` (see below). This column can be empty and then defaults to `false`.

An example using relative paths (same samples as in the example above) could be:

```bash
sample_id,R1,R2,is_sf
sample1,$baseDir/test/sample1_1.fastq.gz,$baseDir/test/sample1_2.fastq.gz,false
sample1,$baseDir/test/sample1a_1.fastq.gz,$baseDir/test/sample1a_2.fastq.gz,false
sample1,$baseDir/test/sample1SF_1.fastq.gz,$baseDir/test/sample1SF_2.fastq.gz,true
sample2,$baseDir/test/sample2_1.fastq.gz,$baseDir/test/sample2_2.fastq.gz,false
```

**Technical replicates** from sequencing the same library on multiple lanes or in individual sequencing runs can be specified in the samplesheet by using the same `sample_id` for these files. In the example above `sample1` has two technical replicates (`sample_1/2.fastq.gz` and `sample1a_1/2.fastq.gz`) which will be merged by the pipeline prior to quantification.

### Reference files for indexing

The quantification requires building an index against reference files. This requires two files, first a **reference genome fasta file** and second a **reference annotation in GTF format**. This is controlled by the options `--genome` and `--gtf`. The defaults are the mouse GENCODE references from version [M25](https://www.gencodegenes.org/mouse/release_M25.html). The user can provide either a download link here (Nextflow will then pull the files automatically) or provide paths to local files. We recommend to download these files manually though as we found the automated staging of files to be unreliable at times.<br>
**Note** that when using non-GENCODE files the user should adjust the options `--gene_name`, `--gene_id` and `--gene_type` which are the names of the columns in the GTF storing these information and probably also `--chrM` (the name of chrM in that annotation) and `-rrna` which is the gene type containing the rRNA information. The chrM and rRNA information are used during the cell barcode whitelisting by Alevin. We currently by default pass the `--gencode` flag to the Alevin indexing process. If non-GENCODE references are used consider removing this by passing `--idx_args ''` to the Nextflow command line.

Optionally, one can build the index and then exit by using the `--idx_only` flag.<br>
Also, one could provide an existing index, e.g. the one build with `--idx_only`. This works via the following params:<br>
- `--idx`: Path to the index folder with the Alevin index files. This is what is outputted in `sc_preprocess_results/alevinIndex/` as folder named `idx_gentrome`
- `--tgmap`: the tx2gene map. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.tx2gene.tsv`
- `--rrna`: gene names of rRNA genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.tx2gene.tsv`
- `--mtrna`: gene names of mitochondrial genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.mtRNA.txt`
- `--expanded_features`: gene names of rRNA genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.features.tsv.gz`
- `--gene2type`: gene names of rRNA genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.gene2type.txt`

In the future we will add an option to simply provide the path to a folder with all these files (so the output of this pipeline when running `--idx_only`).

### Quantification

The defaults assume 10X Chromium V3 libraries which correspond to the `--chromiumV3` flag in Alevin. If using other types of libraries that are supported by Alevin the user should provide the respective Alevin flag via `--quant_args`, e.g. `--quant_args '--dropseq'`, see [Alevin manual](https://salmon.readthedocs.io/en/latest/alevin.html#using-alevin) or the help via `salmon alevin -h` for details on supported library types.

### Feature Barcode Files

If feature barcode experiments shall be quantified then the user must provide a `--features_file` which is a tab-separated file that in column1 stores the name of the feature barcode and in column2 the sequence, e.g.:

```bash
hto_1	ACCCACCAGTAAGAC
hto_2	GGTCGAGAGCATTCA
hto_3	CTTGCCGCATGTCAT
```

The defaults assume 15bp features in R2 with 10X Chromium V3 structure in R1 (16bp CB, 12bp UMI) which corresponds to the Alevin options:<br>

```bash
--bc-geometry 1[1-16] --umi-geometry 1[17-28] --read-geometry 2[11-25]
```

See the [Alevin manual](https://salmon.readthedocs.io/en/latest/alevin.html) for details. 

As mentioned above, the default settings assume 16bp CBs with 12bp UMIs in R1 and 15bp feature barcodes starting from position 11 in R2 as in totalSeqB/C cell hashing kits. Also, by default feature barcode translation (see above, point 5 of "Details") is performed. If not quantifying totalseqB/C experiments but other types of feature barcoding/CITE-seq use `--translate_barcode false`.

### Resources

As we quantify against the expanded and genome-decoyed reference transcriptome the pipeline requires quite some memory. Currently, 30GB of RAM are hardcoded as requirements for indexing and RNA quantification. This works well for mouse samples. We did not test it for human samples. If resources are not sufficient modify the process labels on top of the `nextflow.config` file or use a custom config file via `nextflow -c custom.config`. Indexing and quantification by default use up to 6 CPUs.

### Schedulers

If the pipeline is supposed to run via HPC Schedulers such as `SLURM` we recommend adding the desired configuration into `configs/schedulers.config` (or any other custom Nextflow config file) and then importing that file with the `-c` option, see the [Nextflow docs](https://www.nextflow.io/docs/latest/config.html?highlight=config) for details on config files. The default `-profile slurm` submits a < 8h job to queue `normal`.

### Software

The pipeline is fully containerized via both Docker and Singularity. Use `-profile docker` or `profile singularity` to enable them. Alternatively, use `-profile conda` if none of these container engines are available, but note that containerization ensures fully identical software between runs while conda environments may differ between platforms, e.g. macOS vs Linux distros and upon release of new software versions. As a last option the user can compile all requires software locally, see the `environment.yml` and the module definitions in the `modules` folder for required software, but the latter is not recommended. Use the container engines if possible!

### Output

The pipeline will produce an output folder `sc_preprocess_results` in the location from which the pipeline was launched that contains:

- `alevinIdx`: folder with the expanded transcriptome index(`idx_gentrome`) and the feature barcode index (`idx_features`)
- `alevinQuant`: folder with the alevin outputs (one folder per sample)
- `mtx`: folder with the expression matrices as `mtx.gz` and the column and row annotations as `tsv.gz`. If feature barcode libraries were provided then these files will have a suffix after the basename, by default that is `_SF_`, e.g. `sample1_SF.mtx.gz`. The feature barcode libraries contain all detected/unfiltered cellular barcodes whereas the transcriptomic libraries are already filtered for proper/non-noisy/whitelisted barcodes. That means the user has to remove barcodes (=columns) from these libraries that have no match in the transcriptomic (spliced/unspliced) matrices.
- `alevinQC`: folder with the alevinQC html reports

Depending on `--publishmode` the files in these folders are either real files or soft/hard/relative links to the Nextflow `work` directory.