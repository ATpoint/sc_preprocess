# sc_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

<br>

## Introduction

**sc_preprocess** is an automated preprocessing pipeline for 10x scRNA-seq data implemented using [Nextflow](https://www.nextflow.io/) which is fully containerized to take care of all required software and ensure reproducibility. It optionally supports feature barcode experiments such as CITE-Seq and cell hashing (HTO). The workhorse of this pipeline is the quantification software [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) from [Rob Patro's lab](https://combine-lab.github.io/) and its scRNA-seq module [alevin](https://salmon.readthedocs.io/en/latest/alevin.html). See also the publications for [salmon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/) and [alevin](https://pubmed.ncbi.nlm.nih.gov/30917859/) in Nature Methods and Genome Biology.

## Details

When running with default parameters the following steps will be executed:

1. Generate an expanded reference transcriptome containing all spliced (=exonic) and unspliced (intronic) transcript sequences. This expanded transcriptome is required for quantification to allow output of both **spliced- and unspliced counts**, e.g. for [velocity](https://www.embopress.org/doi/full/10.15252/msb.202110282) analysis. Also, this expanded reference is decoyed by the entire reference genome to perform [selective alignments](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) in order to improve mapping accuracy by capturing reads that better align to the genome than the transcriptome, removing potential gDNA contaminations and other spurious mappings. The indexing of this expanded reference is done by `salmon`.

2. Processing of the reads (fastq files) against the expanded reference using `alevin` performing cellular barcode (CB) detection, read mapping, Unique Molecular Identifier (UMI) deduplication and gene count estimation, resulting in per-cell gene-level abundance estimations.

3. Split the obtained quantifications into spliced and unspliced count tables and save these in [mtx](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) format for easy distribution and loading into downstream analysis environments such as R.

4. Create a per-sample summary report using [alevinQC](https://csoneson.github.io/alevinQC/articles/alevinqc.html) which includes relevant QC metrics such as the average number of reads per cell, number of detected genes per cell [and others](https://bioconductor.org/packages/3.15/bioc/vignettes/alevinQC/inst/doc/alevinqc.html#generate-individual-plots).

5. Optionally, quantify matched reads from a [feature barcoding](https://www.biolegend.com/en-us/blog/cite-seq-and-totalseq-reagents) experiment such as [CITE-seq](https://cite-seq.com/) or [HTO cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) against a provided set of reference barcode sequences. The obtained counts are filtered for the cellular barcodes detected in the RNA quantification. In feature barcode mode the pipeline will eventually return a count matrix that only contains those cellular barcodes detected in both the RNA and feature barcode experiment. By default the pipeline performs [barcode translation](https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-) assuming totalSeqB/C HTOs which can be turned off, see below. We set this as default as in our lab this is the most common type of feature barcoding experiment. The translation tabke is taken from CellRanger the [CellRanger GitHub](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz), see also this thread at [biostars.org](https://www.biostars.org/p/9506747/). It can be turned off with `--translate_barcodes false`, in this case no translation happens and the pipeline assumes that the cellular barcodes between the RNA and feature barcode experiments are the same.

## Usage

The typical command line, e.g. on a HPC with Singularity and SLURM (see sections below) would be:

```bash

NXF_VER=21.10.6 \
    nextflow run main.nf -profile singularity,slurm \
    --genome path/to/genome.fa.gz \
    --gtf path/to/gtf.gz \
    --samplesheet path/to/samplesheet.csv \
    -with-trace -with-report report.html -bg > report.log

```    

For a feature barcoding experiment one would use:

```bash

NXF_VER=21.10.6 \
    nextflow run main.nf -profile singularity,slurm \
    --genome path/to/genome.fa.gz \
    --gtf path/to/gtf.gz \
    --samplesheet path/to/samplesheet.csv \
    --features_file path/to/feature_barcode_file.txt \
    -with-trace -with-report report.html -bg > report.log

```
For details on the individual steps and the input samplesheet see below. The default output folder for all results is called `sc_preprocess_results` generated in the location from which the pipeline was launched.

### Indexing

The quantification requires building an index against reference files. This requires two files, first a **reference genome fasta file** and second a **reference annotation in GTF format**. This is controlled by the options `--genome` and `--gtf`. The defaults are the mouse GENCODE references from version [M25](https://www.gencodegenes.org/mouse/release_M25.html). The user can provide either a download link here (Nextflow will then pull the files automatically) or provide paths to local files. We recommend to download these files manually though as we found the automated staging of files by Nextflow was unreliable at times.

The pipeline then parses the spliced- and unspliced transcript sequences from these files and combines them into an expanded transcriptome. This requires the three flags `--gene_name`, `--gene_id` and `--gene_type` which have defaults to be compatible with GENCODE GTF files. They're used to find the columns in the GTF that store gene_name, gene_id and the gene_type. We also require the flags `--chrM` and `-rRNA`. The chrM flag takes the name of the mitochondrial chromosome, and the rRNA flag takes the gene_type in the GTF that indicates ribosomal genes. Both these information are later used during the cellular barcode whitelisting procedure and can be used in downstream analysis as QC metrics. The defaults for these flags are "chrM" and "rRNA" which is compatible with GENCODE annotations.

 The entire genome from `--genome` is used as a decoy for [selective alignment](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8), meaning that any reads better aligning to the genome rather than transcriptome will be "decoyed" aka removed. This information will be directly incorporated into the index so the user does not need to worry about that once the index is built.

Finally, there is an option `--idx_args` which allows to pass any custom options to the indexing process. Available options can be 
**Note** that when using non-GENCODE files the user should adjust the options `--gene_name`, `--gene_id` and `--gene_type` which are the names of the columns in the GTF storing these information and probably also `--chrM` (the name of chrM in that annotation) and `-rrna` which is the gene type containing the rRNA information. The chrM and rRNA information are used during the cell barcode whitelisting by Alevin. We currently by default pass the `--gencode` flag to the Alevin indexing process. If non-GENCODE references are used consider removing this by passing `--idx_args ''` to the Nextflow command line.

The user can either build an index for every new run or build the index first using the  `--idx_only` flag. In the latter case the pipeline stops after indexing so one could move the index to a permanent storage location and then run the rest of the pipeline separately.
Running with a premade index would require:<br>
- `--idx`: Path to the index folder with the Alevin index files. This is what is outputted in `sc_preprocess_results/alevinIndex/` as folder named `idx_gentrome`
- `--tgmap`: the tx2gene map. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.tx2gene.tsv`
- `--rrna`: gene names of rRNA genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.tx2gene.tsv`
- `--mtrna`: gene names of mitochondrial genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.mtRNA.txt`
- `--expanded_features`: gene names of rRNA genes. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.expanded.features.tsv.gz`
- `--gene2type`: a map connecting gene_name/id to gene_type. This is what is outputted in `sc_preprocess_results/alevinIndex/` as `annotation.gene2type.txt`

Finally, there is `--idx_args` which allows to pass further arguments to `salmon index`. See `salmon index -h` or check the salmon docs for available arguments. By default we pass `--gencode --sparse` telling salmon that we use GENCODE-formatted files and to build a sparse index. The latter saves some memory but makes analysis slightly slower. Use `--idx_args ''` to remove these arguments or pass custom ones. Usually there is no need to modify the defaults.

### Samplesheet

The pipeline reads the fastq file pairs from a [samplesheet](https://github.com/ATpoint/sc_preprocess/blob/main/test/samplesheet.csv) which is a four-column CSV with a header.

```bash
sample_id,R1,R2,is_fb
sample1,/full/path/to/test/sample1_1.fastq.gz,/full/path/to/test/sample1_2.fastq.gz,false
sample1,/full/path/to/test/sample1a_1.fastq.gz,/full/path/to/test/sample1a_2.fastq.gz,false
sample1,/full/path/to/test/sample1SF_1.fastq.gz,/full/path/to/test/sample1SF_2.fastq.gz,true
sample2,/full/path/to/test/sample2_1.fastq.gz,/full/path/to/test/sample2_2.fastq.gz,false
```

The first line is the mandatory header followed by the sample to quantify:

- column1 (`sample_id`) is the per-sample name to be used in the output. This can be any Unix-compatible name, and there is no need to match this name with the fastq file names as other software may require you to do. 
- column2/3 (`R1/2`) are the paths to the fastq files for that sample. It either must be the full **absolute path** (don't use `~`) or alternative a path **relative** to the directory that this pipeline is started from using either of the three [implicit Nextflow variables](https://www.nextflow.io/docs/latest/script.html?highlight=basedir#implicit-variables) `$baseDir`, `$projectDir` and `$launchDir` (see example below).
- column4 (`is_fb`) is a logical column with either `true` or `false` indicating whether this fastq file pair is a feature barcode experiment. If so then these files will be quantified against the barcode library provided by `--features-file` (see below). This column can be empty and then defaults to `false`.

An example using relative paths (same samples as in the example above) could be:

```bash
sample_id,R1,R2,is_fb
sample1,$baseDir/test/sample1_1.fastq.gz,$baseDir/test/sample1_2.fastq.gz,false
sample1,$baseDir/test/sample1a_1.fastq.gz,$baseDir/test/sample1a_2.fastq.gz,false
sample1,$baseDir/test/sample1SF_1.fastq.gz,$baseDir/test/sample1SF_2.fastq.gz,true
sample2,$baseDir/test/sample2_1.fastq.gz,$baseDir/test/sample2_2.fastq.gz,false
```

**Technical replicates** from sequencing the same library on multiple lanes or in individual sequencing runs can be specified in the samplesheet by using the same `sample_id` for these files. In the example above `sample1` has two technical replicates (`sample_1/2.fastq.gz` and `sample1a_1/2.fastq.gz`) which will be merged by the pipeline prior to quantification.

If fastq files are present in multiple directories and the user still wants a single folder with all input fastq files we recommend to create a new folder and use `ln -s` to make symlinks of these files into that new folder. This will collect all files in one place without that the fastq actually need to be moved. **As a remark:** The same strategy makes sense if renaming of files is desired. Be smart, don't touch or manually rename your precious raw data. Do it programatically and save the script or even better, make symlinks to the files and then rename those. That way the original files stay save and untouched.

The pipeline will validate the integrity of the samplesheet and detect if a fastq files was provided in the sheet more than once or if the fastq does not exist. If validation fails the pipeline exists with an error and a message that helps debugging the problem. 

### Feature Barcode Files

If feature barcode experiments are to be quantified then the user must provide a `--features_file` which is a tab-separated list that in column1 stores the name of the feature barcode and in column2 the sequence, e.g.:

```bash
hto_1	ACCCACCAGTAAGAC
hto_2	GGTCGAGAGCATTCA
hto_3	CTTGCCGCATGTCAT
```

As mentioned above, the pipeline by default assumes 10x Chromium V3 libraries with totalSeqB/C feature barcodes which require barcode translation (see point 5 of "Details"). If using any other feature barcoding type then this should be turned off by providing the `--translate_barcode false` flag to the Nextflow command line.
### Read processing

The read processing (=quantification, CB detection, UMI deduplication) process requires a couple of flags. The defaults assume 10x Chromium V3 libraries with CBs/UMIs in read1 and cDNA/features barcodes in read2. <br>

- `--r1_type`: this flag defines structure of read1 in the experiment so the position and length of the CB and UMI. This must consist of the two options `--bc-geometry --umi-geometry` from `alevin`. The defaults are:<br>
`--bc-geometry 1[1-16] --umi-geometry 1[17-28]` and indicate that the BC is in read1 from position 1-16 and the UMI from position 17-28, so 16bp CB and 12bp UMI as in Chromium V3. For Chromium V2 it would be 10bp UMIs so change the UMI flag to `--umi-geometry 1[17-26]`. 
- `--r2_type`: same as above but defining the structure of read2. For Chromium V3 (the default) it would be:<br>
`--read-geometry 2[1-91]` meaning that the first 91bp of R2 should be used for quantification, as recommended by 10x. One could also use `1-end` here so `alevin` would only the entire read, e.g. in case of 150bp reads. We stick with the default of 91bp here fir consistency between runs as not every sequencing run may produce 150bp reads, depending on the machine and run mode.
- `--R2_type_fb`: same as `--r2_type` but the read2 structure for feature barcode experiments. The defaults here assume totalSeqB feature barcoding:<br>
`--read-geometry 2[11-25]` meaning that the feature barcodes are 15bp long and starting from at position 11. The original CITE-seq protocol would be `[1-15]`. 
- `--libtype`: this flag defines the [fragment library type](https://salmon.readthedocs.io/en/latest/library_type.html). For 10x libraries that is `ISR` meaning a stranded paired-end library.
- `--quants_args`: this flag allows to pass further arguments to the `alevin` processing. This can be any of the allowed `alevin` arguments (see its manual) such as generating inferential replicates. Is must **not** include any of `-o -i -p` as these options are already defined internally. See the [Alevin manual](https://salmon.readthedocs.io/en/latest/alevin.html#using-alevin) or the help via `salmon alevin -h` for details on further options.

As said above, the default assume 10x Chromium V3. If using feature barcodes the defaults go with totalSeqB/C sequences.

### Resources

Currently, 30GB of RAM are hardcoded as requirements for indexing and RNA quantification. This is necessary because we index and quantify against the expanded (spliced+unspliced) transcriptome with the entire genome as decoy. This works well for mouse samples. We did not test it for human samples yet. If resources are not sufficient then modify the process labels on top of the `nextflow.config` file or use a custom config file via `nextflow -c custom.config`. Indexing and quantification by default uses up to 6 CPUs. On our HPC the indexing usually takes 60-90' and the processing of a typical sample is usually completed after 2-3h. Note that Nextflow will use all available resources on the machine it runs this pipeline on.

### Job Schedulers

Any [scheduler/executor](https://www.nextflow.io/docs/latest/executor.html) that Nextflow supports can be used here. If the pipeline is supposed to run via `SLURM` we recommend adding the desired configuration into `configs/schedulers.config` (or any other custom Nextflow config file) and then importing that file with the `-c` option, see the [Nextflow docs](https://www.nextflow.io/docs/latest/config.html?highlight=config) for details on config files. The default `-profile slurm` submits a < 8h job to queue `normal`.

### Software

The pipeline is fully containerized via both Docker and Singularity. Use `-profile docker` or `profile singularity` to enable them. Alternatively, use `-profile conda` if none of these container engines are available, but note that containerization ensures fully identical software between runs while conda environments may differ between platforms, e.g. macOS vs Linux distros and upon release of new software versions. As a last option the user can compile all requires software locally, see the `environment.yml` and the module definitions in the `modules` folder for required software, but the latter is not recommended. Use the container engines if possible!

### Output

The pipeline will produce an output folder `sc_preprocess_results` in the location from which the pipeline was launched that contains:

- `alevinIdx`: folder with the expanded transcriptome index(`idx_gentrome`) and the feature barcode index (`idx_features`)
- `alevinQuant`: folder with the alevin outputs (one folder per sample)
- `mtx`: folder with the expression matrices as `mtx.gz` and the column and row annotations as `tsv.gz`. If feature barcode (FBs) libraries were present for a particular sample then the returned files will already be filtered for CBs found in both the RNA and FB experiment. The FB counts will be appended to the `mtx.gz` so will be the last entries in these files. The same goes for the `sample_feature.tsv.gz` file, where the feature barcode names will be the last entries.<br>
- `alevinQC`: folder with the alevinQC html reports for each quantified library -- RNA and FB (if present). Note that in case of FBs this QC report will be based on the full RNA/FB experiment that has not been filtered for CBs present in both experiments. It is useful to judge the quality of both experiments independently.<br>
In this folder there will also be a file `summary_cellnumbers.txt` which summarizes the number of detected cells per sample in the RNA experiment and the number of intersecting cells with the FB experiment. If no FBs were present for that sample NAs are returned. The numbers are also presented as a barplot in `summary_cellnumbers.pdf`.

The typical files required for downstream analysis are the `sample.mtx.gz` and associated barcode/features files. In case of feature barcoding that file will contain the FB quantifications at the bottom of that file. "Gene" names are 

Depending on `--publishmode` the files in these folders are either real files or soft/hard/relative links to the Nextflow `work` directory.
See the Nextflow docs on [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) for details. The default is mode `copy` so the output folder contains real files rather than links. 
