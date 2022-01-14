# sc_preprocess

A Nextflow pipeline for preprocessing of 10X scRNA-seq data using [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) supporting output of spliced- and unspliced counts for velocity analysis as well as feature barcoding experiments (CITE-seq).

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

Documentation will follow...

The pipeline expects a samplesheet as input. It must be a CSV file with four fields and a header as below:

```bash
sample_name,R1,R2,is_sf
```

The first field `sample_name` is the sample name to be used in the output folder. This can be anything the user feels good with. Unlike CellRanger the pipeline does not expect any naming conventions for either sample names or fastq files. Be smart and do not use whitespaces here. The 2n and 3rd fields are `R1` and `R2` which are the paths to the R1 and R2 fastq files. Lastly, the 4th field `is_sf` should be `true` or `false`. It is `true` if that particular fastq pair is a feature barcode sample (e.g. ADT/HTO) and `false` or just empty if it is the "normal" transcriptomic readout. Lane replicates should be specified by using **the same `sample_id`** in different lines, see examples below. Also, the feature barcode reads of a certain sample should have the same `sample_id` as the transcriptomic readout fastq files.

Example:

```bash
sample_id,R1,R2,is_sf
sample1,path/to/sample1_1.fastq.gz,path/to/sample1_2.fastq.gz,false
sample1,path/to/sample1a_1.fastq.gz,path/to/sample1a_2.fastq.gz,false
sample1,path/to/sample1SF_1.fastq.gz,path/to/sample1SF_2.fastq.gz,true
```

Here we have a single sample called `sample1` which has two fastq file pairs for the transcriptomic readout (lane replicates). As these share the same `sample_id` the pipeline will automatically merge them prior to quantification. That sample has one fastq file pair with feature barcode reads indicated by `true` in the 4th field. That pair has the same `sample_id` and in the output folder will be automatically named `sample1_SF`, note the `_SF`, so it can easily be identified as feature barcode (SF=surface feature) quantificaitons. Empty lines in the samplesheet are allowed to visually group samples.


In this `sample_name` is the **group** name of the respective sample. It can be something like `treatment_rep1`. Lane replicates must have the same `sample_name`and then will internally be combined. The `fastq_1/2` is then the path (relative to samplesheet location) to the fastq files. Finally, `is_sf` is a logical (`true`/`false`) indicating whether that sample is a surface feature (HTO/ADT) or not. If true then one needs to provide a tab-separated file with the feature names and sequences, e.g.:

```bash
hto_1	ACCCACCAGTAAGAC
hto_2	GGTCGAGAGCATTCA
hto_3	CTTGCCGCATGTCAT
```

...which will then be used to quantify the SF libraries against. THe fourth field of the samplesheet can be left empty and default then to `false`.



## Input

Samplesheet, it is recommended to just make symlinks.

