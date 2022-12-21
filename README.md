# sc_preprocess

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.4-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000&logo=data%3Aimage%2Fjpeg%3Bbase64%2C%2F9j%2F4AAQSkZJRgABAQAAkACQAAD%2F2wBDABwcHBwcHDAcHDBEMDAwRFxEREREXHRcXFxcXHSMdHR0dHR0jIyMjIyMjIyoqKioqKjExMTExNzc3Nzc3Nzc3Nz%2F2wBDASIkJDg0OGA0NGDmnICc5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ub%2FwAARCACAAHgDASIAAhEBAxEB%2F8QAGgAAAgMBAQAAAAAAAAAAAAAABAUAAgMBBv%2FEADMQAAIBAwIDBgYBAwUAAAAAAAECAAMEESExEpHRBRVBUVJxEyIyYYGxQhQzwUNTgpLh%2F8QAGAEAAwEBAAAAAAAAAAAAAAAAAAIDAQT%2FxAAfEQEBAQACAgMBAQAAAAAAAAAAAQIDESExEkFRIjL%2F2gAMAwEAAhEDEQA%2FADe5rX1PzHSTua19T8x0jYbTsAUdzWvqfmOknc1r6n5jpG8yq1kpDLn8eMAW9zWvqfmOko3ZNmmrOw9yOk0qXdR9F%2BUfbeCHU5OspOP9SvJ%2BIbHs4f6jn2x0lk7OsnVmDVMLvt0lYZb%2FANqr7CG8yS0Z3bQw7OsD%2FNx746TVeybNvpdz%2BR0kk%2B4nLOSn7W7mtfU%2FMdJO5rX1PzHSbJcOujfMIalRagyplM6lb2WdzWvqfmOknc1r6n5jpG8kZpR3Na%2Bp%2BY6SRvJAODadnBtMq9UUU4vHwEILelLi4FEYGrGKWZnPExyTOMxYlmOSZtRoNWbTQDcy8kzHPbdVkqM54UGTDqdid6h%2FA6w6nTSmvCgxLxLu%2FSk459h1taA%2Fjn31mop01BCqADvLznEvmIndp%2BoyNCkf449pg9r4ofwYZkHYzsS4lHRQyspwwwZFYqeJTgxo6K4wwzF1Wk1I66jwMjrHXllg2jWFQYOjTeJwSpBG4jOlUFRc%2BPjKY334rZWskkko1wbRPc1fi1DjYaCMq7%2FDolhvsPzEsrxz7S5L9NKdM1XCL4x0iLTUIuwgtlT4U%2BId2%2FUNi7vd6bjPU7Ud0poXc4UakmILjtaq5K244F8zqeXhO9r1yai242UcR9ztEsMw9rZ69aocu7H8n%2FGJjgeIEk0p06lVxTpqWY%2BAjMZj5dV09tP1N0ubin9FRh%2Bc%2FvM2fs%2B8QZamSPsQYHAHFn2jd1KyUWw%2FEcZIwceJ0noGUOpVtjEHY9HiqPXP8Rwj3OpnoZPTYUuhRipl6L8DjyOhhNymU4xuP1AJzWfGsOZJlRbjpgyS8vZgt63yIvmcxbvp5w69yXT2MFpqTUUHzE6M%2BnPvzo7VQqhR4CWkkkHQ8l2ln%2BtqZ%2B2PbEBjzti3IIul20Vv8RJKy%2BC1yF2d0bSr8Th4gRgjx%2FEEkgHsaF7bXGlNhnyOh5TC9sEuV40wKg2Pn9jPKw%2B37QuaGnFxr5N13i%2FH8b29BYUDb2qIww27e5hkGtbqndoXp5GDgg%2BBhMStcYcSlT4xRtpHEU1Bh2H3MlyRlF2p0ZfIySlr9Te0kbHoRtVdkxjxEzFZiQDiXrj5QfKCx2mckqp4lB85aAUcIylKgBDaYPjEdz2QQS1qdPS3%2BDL9sViPh0VODniOPttNLPtRKgFO5PC%2B3F4HoY079sIKlKpRbhqqUP36zOe3dadVOFwGU%2BeonjrhaaV6iUjlAcCNL2yxjJJJGAuzuGt7hXB0JCsPsek9hPCHY4nuV%2Bke0TTYtF7XNQMQMYz5Q52CqWPgIonNy6666V4537H29V6hIbGkk5aD5S3mZI%2BP8%2BS79%2BBDrxIRF8ZDaB1k4WyNjHK1oNkcB8IRFqsVIYeEPVhUXIgHkr6r8a6dxsDwj2H%2FALmCR1V7GqDWjUDfZtDzHSBP2fepvTJ9iDKSwoMEgcIJA8gTicmzW9wv1UnH%2FEyvwq3%2B2%2F8A1MYM5IStpdP9NJvyMfuH0eyKznNdgg8hqekzsA7G3NzcKP4qQzH22H5nr5jRoUrdPh0hgfv3l6jimvEZPWvs0ga6fAFMeOpgMszF2LNuZtb0%2BN8nZZyW%2FLTon8wdSTgphZJpJOmTrw564NpV1Drgyw2nZoLmUqcGdRyhyIZUphxrv5wJlKHDQA5KiuNOUvFoJGomy12H1awAySYCuh3yJf4qFS2dBvANJIMbqmNsmYPdO2i%2FL%2B4l5Mw8xaMqVVpjXfyi2pUao2W5ShJJydZZEZzwqJDW7rwrnMy4ql2CruY1poKahRK0qS0h5k7may3HjrzUt67SSSSUI4Np2J%2B%2Bbb0vyHWd75tvS%2FIdYA3lWUMMMMxV3zbel%2BQ6yd823pfkOsAMegRqmswKlfqGJl3zbel%2BQ6znfNqf4PyHWAazVf7VT2gR7VsjvTbkOsnetngrwPg77dZlnhsWnVVm%2BkZmY7TsRtTbkOs075tRsj8h1kZw%2FtVvJ%2BCUtWOtQ4%2BwhqoqDCjAirvm29L8h1k75tvS%2FIdZXOJPSd1abyRR3zbel%2BQ6yd823pfkOsYpvJE%2FfNt6X5DrJAP%2F2Q%3D%3D)](https://sylabs.io/docs/) 
  
==> [Introduction](#Introduction)  
==> [Details](#Details)  
==> [Usage](#Usage)  
====> [Indexing](#Indexing)  
====> [Quantification](#Quantification)    
==> [Output](#Output)  

## Introduction

**sc_preprocess** is a containerized preprocessing pipeline for 10x scRNA-seq data written in [Nextflow](https://www.nextflow.io/).
The pipeline basically implements the code suggestions of the salmon/alevin developers in [their tutorial](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/).
This covers generation of a genome-decoyed and expanded (spliced+unspliced) transcriptome index directly from reference annotations, quantification of reads against this index, generation of spliced- and unspliced count tables and a basic QC summary report. It optionally supports feature barcoding experiments such as CITE-Seq or cell hashtag oligos (HTO). The indexing procedure relies on extraction of spliced and unspliced counts from a genome/GTF using [eisaR](https://bioconductor.org/packages/release/bioc/html/eisaR.html) followed by index building with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html). The quantification happens via [alevin](https://salmon.readthedocs.io/en/latest/alevin.html). Generation of count tables from the quantification results is achieved via [tximeta](https://bioconductor.org/packages/release/bioc/html/tximeta.html). A QC summary report is provided by [alevinQC](https://www.bioconductor.org/packages/release/bioc/html/alevinQC.html).

## Workflow

The pipeline covers these steps:

1. Generate a genome-decoyed and expanded reference transcriptome containing all spliced (=exonic) and unspliced (intronic) transcript sequences. This expanded transcriptome is required for quantification to allow output of both **spliced- and unspliced counts**, e.g. for [velocity](https://www.embopress.org/doi/full/10.15252/msb.202110282) analysis. Also, this expanded reference is decoyed by the entire reference genome to perform [selective alignments](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) in order to improve mapping accuracy by capturing reads that better align to the genome than the transcriptome, removing potential gDNA contaminations and other spurious mappings. The indexing of this expanded reference is done by `salmon`.

2. Processing of the reads (fastq files read from a CSV samplesheet) against the expanded reference using `alevin` performing cellular barcode (CB) detection, read mapping, Unique Molecular Identifier (UMI) deduplication and gene count estimation, resulting in per-cell gene-level abundance estimations.

3. Split the obtained quantifications into spliced and unspliced count tables and save these in [mtx](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) format for easy distribution and loading into downstream analysis environments such as R.

4. Create a per-sample summary report using [alevinQC](https://csoneson.github.io/alevinQC/articles/alevinqc.html) which includes relevant QC metrics such as the average number of reads per cell, number of detected genes per cell [and others](https://bioconductor.org/packages/3.15/bioc/vignettes/alevinQC/inst/doc/alevinqc.html#generate-individual-plots). Also, a table and plot summarizing number of cells per sample is generated using custom scripts.

5. Optionally, quantify matched reads from a [feature barcoding](https://www.biolegend.com/en-us/blog/cite-seq-and-totalseq-reagents) experiment such as [CITE-seq](https://cite-seq.com/) or [HTO cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) against a provided set of reference barcode sequences. The obtained counts are filtered for the cellular barcodes detected in the RNA quantification. In feature barcode mode the pipeline will eventually return a count matrix that only contains those cellular barcodes detected in both the RNA and feature barcode experiment. Some feature barcode sets such as totalSeq-B/C require [barcode translation](https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-). We provide the default translation table from the [CellRanger GitHub](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz) in the `assets/` folder to perform this translation if `--translate_barcodes` is set and `--translate_table` points to that table. See also this thread at [biostars.org](https://www.biostars.org/p/9506747/) for details.

## Usage

Typically one would build an index first and then run the quantification separately as shown below. Still, the entire workflow
including indexing and quantification can also be run as a single command.

### Indexing

Typical command line for the indexing on a cluster with SLURM and Singularity, using GENCODE annotations. If not using GENCODE then the `--idx_args` argument can be omitted.

```bash
NXF_VER=22.10.4 nextflow run main.nf -profile singularity,slurm --idx_only --genome path/to/genome.fa.gz --gtf path/to/gtf.gz --idx_args '\--gencode'
```    

This will parse both spliced and unspliced transcript sequences directly from the genome using the GTF as guide and then run the indexing
process on these files, using the provided genome as decoy. Six CPUs and 30GB of RAM are hardcoded in the default profile for this step,
and this is sufficient for mouse data.

**Indexing parameters**

These parameters can be used to customize the indexing:  

`--only_idx`: logical, whether to only run the indexing process without any quantification and then exit  
`--genome`: path to the genome fasta file (`.fa.gz`)   
`--gtf`: path to the GTF reference file (`.gtf.gz`)  
`--gene_id`: name of the column storing the gene ID in the GTF, default `gene_id`   
`--gene_name`: name of the column storing the gene name in the GTF, default `gene_name`  
`--gene_type`: name of the column storing the gene biotype in the GTF, default `gene_type`    
`--chrM`: name of the mitochondrial chromosome, default `chrM`  
`--rrna`: gene biotype for ribosomal RNAs, default `rRNA`  
`--read_length`: an integer indicating the read length of R2 (the non-CB/non-UMI read), default 150
`--idx_outdir`: name of output folder inside the `sc_preprocess_results` directly to store the output  
`--idx_name`: name of the index in which the actual index is stored inside `--idx_outdir`, default `idx`  
`--idx_args`: additional arguments passed to `salmon index` step, default `--sparse`  

### Quantification

Typical command line for the quantification on a HPC with SLURM and Singularity, assuming presence of feature barcodes that require barcode translation:  

```bash
NXF_VER=22.10.4 nextflow run main.nf -profile singularity,slurm \
    --samplesheet path/to/samplesheet.csv \
    --idx path_to_idx_folder \
    --tgmap path_to_tgmap --rrnagenes path_to_rrnagenes_file --mtgenes path_to_mtrnagenes_file \
    --features path_to_expanded_features_file --gene2type path_to_gene2type_file \
    --features_file path/to/hto.tsv --translate_table path/to/assets/translate_table.txt.gz --translate_barcodes
```    

This will quantify the reads against the index, create spliced- and unspliced count tables and create the alevinQC report. As resources we hardcoded 6 CPUs and 30GB of RAM which works well for mouse samples of typical size. See details below.  

**Samplesheet**

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
- column4 is the [library type](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype) so basically the strandedness and layout of the library. That would be `ISR` for 10x Chromium data.
- column5 (`is_fb`) is a logical column with either `true` or `false` indicating whether this fastq file pair is a feature barcode experiment. If so then these files will be quantified against the barcode library provided by `--features-file` (see below). This column can be empty and then defaults to `false`. See below for details on feature barcoding experiments.

An example using relative paths (same samples as in the example above) could be:

```bash
sample,r1,r2,libtype,is_fb
sample1,$baseDir/test/sample1_1.fastq.gz,$baseDir/test/sample1_2.fastq.gz,ISR,false
sample1,$baseDir/test/sample1a_1.fastq.gz,$baseDir/test/sample1a_2.fastq.gz,ISR,false
sample1,$baseDir/test/sample1SF_1.fastq.gz,$baseDir/test/sample1SF_2.fastq.gz,ISR,true
sample2,$baseDir/test/sample2_1.fastq.gz,$baseDir/test/sample2_2.fastq.gz,ISR,false
```

**Technical replicates**   

...from sequencing the same library on multiple lanes or in individual sequencing runs can be specified in the samplesheet by using the same `sample_id` for these files. In the example above `sample1` has two technical replicates (`sample_1/2.fastq.gz` and `sample1a_1/2.fastq.gz`) which will be merged by the pipeline prior to quantification. 

**Feature barcoding experiments**

If feature barcode experiments are to be quantified then the user must provide a `--features_file` which is a tab-separated list that in column1 stores the name of the feature barcode and in column2 the sequence, for example:  

```bash
hto_1	ACCCACCAGTAAGAC
hto_2	GGTCGAGAGCATTCA
hto_3	CTTGCCGCATGTCAT
```

If this file is provided then an index will be made for these files and the feature barcode fastq files will be quantified against it.
The counts of the feature barcodes will be included to the bottom of the gene expression count tables. If barcode translation must be performed (totalSeq-B/C),
then use `--translate_barcodes --translate_table path/to/translate_table.txt.gz`.

**Quantification parameters**

These parameters can be used to customize the quantification:  

- `--r1_type`: this flag defines structure of read1 in the experiment so the position and length of the CB and UMI. This must consist of the two options `--bc-geometry --umi-geometry` from `alevin`. The defaults are:  
`--bc-geometry 1[1-16] --umi-geometry 1[17-28]` and indicate that the BC is in read1 from position 1-16 and the UMI from position 17-28, so 16bp CB and 12bp UMI as in Chromium V3. For Chromium V2 it would be 10bp UMIs so one would use `--r1_type '\--bc-geometry 1[1-16] --umi-geometry 1[17-26]'`.  
- `--r2_type`: same as above but defining the structure of read2. For Chromium V3 (the default) it would be:<br>
`--r2_type '\--read-geometry 2[1-91]'` meaning that the first 91bp of R2 should be used for quantification, as recommended by 10x. One could also use `1-end` here so `alevin` would only the entire read, e.g. in case of 150bp reads. We stick with the default of 91bp here fir consistency between runs as not every sequencing run may produce 150bp reads, depending on the machine and run mode.  
- `--R2_type_fb`: same as `--r2_type` but the read2 structure for feature barcode experiments. The defaults here assume totalSeqB feature barcoding:<br>
`--read-geometry 2[11-25]` meaning that the feature barcodes are 15bp long and starting from at position 11. The original CITE-seq protocol would be `[1-15]`.  
- `--quants_args`: this flag allows to pass further arguments to the `alevin` processing. This can be any of the allowed `alevin` arguments (see its manual) such as generating inferential replicates. Is must **not** include any of `-o -i -p` as these options are already defined internally. See the [Alevin manual](https://salmon.readthedocs.io/en/latest/alevin.html#using-alevin) or the help via `salmon alevin -h` for details on further options.  

## Output

The pipeline will produce an output folder `sc_preprocess_results` in the location from which the pipeline was launched that contains:

- `alevin_idx`: folder with the expanded transcriptome index (`idx_gentrome`) and the feature barcode index (`idx_features`)
- `alevin_quant`: folder with the alevin outputs (one folder per sample)
- `mtx`: folder with the expression matrices as `mtx.gz` and the column and row annotations as `tsv.gz`. If feature barcode (FBs) libraries were present for a particular sample then the returned files will already be filtered for CBs found in both the RNA and FB experiment. The FB counts will be appended to the `mtx.gz` so will be the last entries in these files. The same goes for the `sample_feature.tsv.gz` file, where the feature barcode names will be the last entries.<br>
- `qc_dir`: folder with the alevinQC html reports for each quantified library -- RNA and FB (if present). Note that in case of FBs this QC report will be based on the full RNA/FB experiment that has not been filtered for CBs present in both experiments. It is useful to judge the quality of both experiments independently.<br>
In this folder there will also be a file `summary_cellnumbers.txt` which summarizes the number of detected cells per sample in the RNA experiment and the number of intersecting cells with the FB experiment. If no FBs were present for that sample NAs are returned. The numbers are also presented as a barplot in `summary_cellnumbers.pdf`.
- `pipeline_info`: folder contains `software_versions.txt` listing all software versions used in the pipeline and `command_lines.txt` listing the exact command lines used in the processes  
