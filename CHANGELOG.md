# Changelog

## v2.0
- updated with nf_blank template
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
## v1.1.0
- added multiqc
- now use micromamba base image
- make environment platform agnostic

## v1.0.0
- initial versionl
