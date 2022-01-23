# Changelog

## v2.2
- added error strategy to each module, now submitted processes will finish if any process fails rather than killing the entire pipeline
- added `--keepCBFraction 1.0` flag to the FB quantification process
- the summary table/plot towards cell numbers in the `alevinQC` dir now only lists cells per RNA experiment and intersect with the FB experiment
- fixed minor bug related to barcode translation in the `mtx.R` script that happened when `alevin` identified CBs not listed in the translation list
## v2.1
- extended documentation in `README.md`
- now output one QC report for RNA and feature barcode (FB) experiments
- support barcode translation for totalSeqB/C with the CellRanger 3M whitelist
- changed params for the quantification process, the user can now freely define the CB/UMI/read geometry via the params `--r1_type`, `-r2_type` and `--r2_type_fb` rather than using the in-build flags like `--chromiumV3`, this makes the pipeline more generic
- We now only output one set of mtx files per sample (spliced+unspliced) with the counts of the FBs appended to it. Same goes for barcodes/features in the `tsv.gz` files. These mtx files will represent only cells detected in both the RNA and FB experiment
- summarize number of detected cells in RNA and FB and its intersect as sumamry table and barplot in the alevinQC directory as `summary_detected_cells.txt/pdf`
## v2.0
- updated with [nf_blank](https://github.com/ATpoint/nf_blank) template
  - params validation with params listed in `schema.nf`
  - intro message with summary of all params upon startup of pipeline
  - validate minimal NF version
- add support for feature barcoding experiments
- read fastq file pairs directly from a samplesheet rather than from `--fastq` (retired `--fastq`)
- most modules now have a hardcoded container from quay
- the custom container via `params.container` is now only used for the modules that require multiple Bioc packages simultaneously
- added a CI test for singularity
- added alevinQC module and an overall summary report of all samples based on its output
- allow premade index
- validate that all input files and samplesheet exist
## v1.1.0
- added multiqc
- now use micromamba base image
- make environment platform agnostic
## v1.0.0
- initial versionl
